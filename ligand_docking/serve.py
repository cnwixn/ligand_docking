"""Serve the viewer locally and let it run AutoDock Vina from the page.

A static HTML file cannot start a docking run, so this module starts a small
HTTP server on the loopback interface: the page loads from it, and its Run
button posts the two PDBQT files to ``/api/dock``, which executes Vina and sends
the poses back.
"""

from __future__ import annotations

import base64
import json
import shutil
import threading
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

from .dock import VinaNotFound, run_vina
from .pdbqt import pdbqt_to_pdb
from .prep import PrepError, prepare_ligand, prepare_receptor
from .report import ASSETS, render_html

LOCK = threading.Lock()


def _numbers(values, name: str) -> tuple[float, float, float]:
    try:
        triple = tuple(float(v) for v in values)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be three numbers.") from None
    if len(triple) != 3:
        raise ValueError(f"{name} must be three numbers.")
    return triple  # type: ignore[return-value]


class DockingHandler(BaseHTTPRequestHandler):
    server_version = "LigandDocking/0.1"
    page: dict = {}          # payload rendered into index.html
    work_root: Path = Path(".")
    vina_path: str | None = None

    def log_message(self, fmt, *args):  # keep the console quiet
        pass

    def _send(self, code: int, body: bytes, content_type: str):
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, code: int, payload: dict) -> None:
        self._send(code, json.dumps(payload).encode("utf-8"), "application/json; charset=utf-8")

    def do_GET(self):  # noqa: N802
        route = urlparse(self.path).path
        if route in ("/", "/index.html"):
            self._send(200, render_html(self.page).encode("utf-8"), "text/html; charset=utf-8")
            return
        name = route.lstrip("/")
        if name in ("ngl.js", "parser.js"):
            target = ASSETS / name
            ctype = "application/javascript; charset=utf-8"
            self._send(200, target.read_bytes(), ctype)
            return
        self._json(404, {"ok": False, "error": "not found"})

    def do_POST(self):  # noqa: N802
        if urlparse(self.path).path != "/api/dock":
            self._json(404, {"ok": False, "error": "not found"})
            return
        try:
            length = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            length = 0
        if length <= 0 or length > 64 * 1024 * 1024:
            self._json(400, {"ok": False, "error": "bad request size"})
            return
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8"))
        except Exception as exc:  # noqa: BLE001
            self._json(400, {"ok": False, "error": f"bad JSON: {exc}"})
            return

        try:
            result = self._dock(body)
        except ValueError as exc:
            self._json(400, {"ok": False, "error": str(exc)})
            return
        except VinaNotFound as exc:
            self._json(500, {"ok": False, "error": str(exc)})
            return
        except Exception as exc:  # noqa: BLE001
            self._json(500, {"ok": False, "error": f"{type(exc).__name__}: {exc}"})
            return
        self._json(200, result)

    def _dock(self, body: dict) -> dict:
        receptor_b64 = body.get("receptor_b64")
        ligand_b64 = body.get("ligand_b64")
        if not receptor_b64 or not ligand_b64:
            raise ValueError("Both a receptor and a ligand PDBQT are required.")
        center = _numbers(body.get("center"), "centre")
        size = _numbers(body.get("size"), "size")
        exhaustiveness = int(body.get("exhaustiveness") or 32)
        if exhaustiveness < 1:
            raise ValueError("exhaustiveness must be at least 1.")

        work = self.work_root / f"run-{uuid.uuid4().hex[:10]}"
        work.mkdir(parents=True, exist_ok=True)
        receptor_raw = work / Path(body.get("receptor_name") or "receptor.pdbqt").name
        ligand_raw = work / Path(body.get("ligand_name") or "ligand.pdbqt").name
        receptor_raw.write_bytes(base64.b64decode(receptor_b64))
        ligand_raw.write_bytes(base64.b64decode(ligand_b64))

        try:
            receptor = prepare_receptor(receptor_raw, work / "receptor.pdbqt")["pdbqt"]
            ligand = prepare_ligand(ligand_raw, work / "ligand.pdbqt")["pdbqt"]
        except PrepError as exc:
            raise ValueError(str(exc)) from None

        with LOCK:
            outcome = run_vina(receptor, ligand, work / "out.pdbqt",
                               center=center, size=size, exhaustiveness=exhaustiveness,
                               vina=self.vina_path)

        poses = outcome["poses"]
        receptor_pdb = pdbqt_to_pdb(receptor.read_text(encoding="utf-8", errors="ignore"))
        return {
            "ok": True,
            "ligand_name": ligand.name,
            "receptor_pdb": receptor_pdb,
            "poses": [pose.as_dict() for pose in poses],
            "log": outcome["log"][-4000:],
        }


def run_server(page: dict, *, host: str = "127.0.0.1", port: int = 8766,
               vina: str | None = None, open_browser: bool = False,
               workdir: str | Path | None = None) -> None:
    # Vina is a separate process and must be able to write its output here, so
    # the work folder lives somewhere ordinary (the system temp folder is not
    # writable for child processes in every environment); --workdir overrides it
    work_root = Path(workdir).expanduser().resolve() if workdir else Path.cwd() / ".ligand-docking-tmp"
    work_root.mkdir(parents=True, exist_ok=True)
    DockingHandler.page = page
    DockingHandler.work_root = work_root
    DockingHandler.vina_path = vina

    server = ThreadingHTTPServer((host, port), DockingHandler)
    url = f"http://{host}:{port}/"
    print(f"Docking viewer running at {url}")
    print("Press Ctrl+C to stop.")
    if open_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        server.server_close()
        shutil.rmtree(work_root, ignore_errors=True)
