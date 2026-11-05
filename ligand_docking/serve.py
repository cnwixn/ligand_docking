"""Serve the viewer locally, keep docking projects and let the page run Vina.

A static HTML file cannot start a docking run, so this module serves the page
from a loopback HTTP server: the Run button posts the structures to ``/api/dock``,
which prepares them, runs Vina, saves the result as a project folder and sends
the poses back. ``/api/projects`` lists what was saved so old runs can be reopened.
"""

from __future__ import annotations

import base64
import json
import shutil
import threading
import time
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

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


def _safe_name(name: str) -> str:
    keep = "".join(ch if ch.isalnum() or ch in "-_. " else "_" for ch in name).strip()
    return keep[:60] or "docking"


class DockingHandler(BaseHTTPRequestHandler):
    server_version = "LigandDocking/0.2"
    page: dict = {}
    projects_root: Path = Path("./projects")
    work_root: Path = Path("./.ligand-docking-tmp")
    vina_path: str | None = None

    def log_message(self, fmt, *args):
        pass

    # ---------------------------------------------------------------- helpers
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

    def _read_json(self) -> dict:
        try:
            length = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            length = 0
        if length <= 0 or length > 64 * 1024 * 1024:
            raise ValueError("bad request size")
        return json.loads(self.rfile.read(length).decode("utf-8"))

    # ------------------------------------------------------------------ routes
    def do_GET(self):  # noqa: N802
        parsed = urlparse(self.path)
        route = parsed.path

        if route in ("/", "/index.html"):
            self._send(200, render_html(self.page).encode("utf-8"), "text/html; charset=utf-8")
            return

        if route == "/api/projects":
            self._json(200, {"ok": True, "projects": self._list_projects()})
            return

        if route == "/api/project":
            project_id = (parse_qs(parsed.query).get("id") or [""])[0]
            try:
                self._json(200, {"ok": True, **self._load_project(project_id)})
            except (FileNotFoundError, ValueError) as exc:
                self._json(404, {"ok": False, "error": str(exc)})
            return

        name = route.lstrip("/")
        if name in ("ngl.js", "parser.js"):
            self._send(200, (ASSETS / name).read_bytes(), "application/javascript; charset=utf-8")
            return

        self._json(404, {"ok": False, "error": "not found"})

    def do_POST(self):  # noqa: N802
        if urlparse(self.path).path != "/api/dock":
            self._json(404, {"ok": False, "error": "not found"})
            return
        try:
            body = self._read_json()
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

    # ------------------------------------------------------------------- logic
    def _dock(self, body: dict) -> dict:
        receptor_b64 = body.get("receptor_b64")
        ligand_b64 = body.get("ligand_b64")
        if not receptor_b64 or not ligand_b64:
            raise ValueError("Both a receptor and a ligand are required.")
        center = _numbers(body.get("center"), "centre")
        size = _numbers(body.get("size"), "size")
        exhaustiveness = int(body.get("exhaustiveness") or 8)
        if exhaustiveness < 1:
            raise ValueError("effort must be at least 1.")
        name = _safe_name(body.get("project_name") or "")

        run_dir = self.work_root / f"run-{uuid.uuid4().hex[:10]}"
        run_dir.mkdir(parents=True, exist_ok=True)
        receptor_raw = run_dir / Path(body.get("receptor_name") or "receptor.pdbqt").name
        ligand_raw = run_dir / Path(body.get("ligand_name") or "ligand.pdbqt").name
        receptor_raw.write_bytes(base64.b64decode(receptor_b64))
        ligand_raw.write_bytes(base64.b64decode(ligand_b64))

        try:
            receptor = prepare_receptor(receptor_raw, run_dir / "receptor.pdbqt")["pdbqt"]
            ligand = prepare_ligand(ligand_raw, run_dir / "ligand.pdbqt")["pdbqt"]
        except PrepError as exc:
            raise ValueError(str(exc)) from None

        with LOCK:
            outcome = run_vina(receptor, ligand, run_dir / "out.pdbqt",
                               center=center, size=size, exhaustiveness=exhaustiveness,
                               vina=self.vina_path)

        poses = outcome["poses"]
        if not poses:
            raise ValueError("Vina produced no poses; check the box and the inputs.")

        project = self._save_project(name, receptor, ligand, run_dir, poses, center, size,
                                     exhaustiveness)
        return {
            "ok": True,
            "project": project,
            "receptor_pdb": pdbqt_to_pdb(receptor.read_text(encoding="utf-8", errors="ignore")),
            "poses": [pose.as_dict() for pose in poses],
            "log": outcome["log"][-4000:],
        }

    def _save_project(self, name: str, receptor: Path, ligand: Path, run_dir: Path,
                      poses, center, size, exhaustiveness) -> dict:
        stamp = time.strftime("%Y%m%d-%H%M%S")
        project_id = f"{stamp}-{uuid.uuid4().hex[:6]}"
        folder = self.projects_root / f"{project_id}-{name}"
        folder.mkdir(parents=True, exist_ok=True)
        for source, target in ((receptor, "receptor.pdbqt"), (ligand, "ligand.pdbqt"),
                               (run_dir / "out.pdbqt", "poses.pdbqt")):
            if source.is_file():
                shutil.copy2(source, folder / target)
        meta = {
            "id": project_id,
            "name": name,
            "created": time.strftime("%Y-%m-%d %H:%M:%S"),
            "center": list(center),
            "size": list(size),
            "exhaustiveness": exhaustiveness,
            "poses": [{"mode": pose.mode, "affinity": pose.affinity,
                       "rmsd_lb": pose.rmsd_lb, "rmsd_ub": pose.rmsd_ub} for pose in poses],
        }
        (folder / "project.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
        return meta

    def _list_projects(self) -> list[dict]:
        items = []
        if self.projects_root.is_dir():
            for folder in sorted(self.projects_root.iterdir(), reverse=True):
                meta_file = folder / "project.json"
                if meta_file.is_file():
                    try:
                        items.append(json.loads(meta_file.read_text(encoding="utf-8")))
                    except json.JSONDecodeError:
                        continue
        return items

    def _load_project(self, project_id: str) -> dict:
        if not project_id:
            raise ValueError("missing project id")
        for folder in self.projects_root.iterdir() if self.projects_root.is_dir() else []:
            if folder.name.startswith(project_id):
                meta = json.loads((folder / "project.json").read_text(encoding="utf-8"))
                receptor_pdbqt = folder / "receptor.pdbqt"
                poses_file = folder / "poses.pdbqt"
                from .pdbqt import parse_vina_out
                poses = parse_vina_out(poses_file.read_text(encoding="utf-8", errors="ignore"))
                return {
                    "project": meta,
                    "receptor_pdb": pdbqt_to_pdb(
                        receptor_pdbqt.read_text(encoding="utf-8", errors="ignore")),
                    "poses": [pose.as_dict() for pose in poses],
                }
        raise FileNotFoundError(f"project {project_id} not found")


def run_server(page: dict, *, host: str = "127.0.0.1", port: int = 8766,
               vina: str | None = None, open_browser: bool = False,
               workdir: str | Path | None = None,
               projects: str | Path | None = None) -> None:
    # Vina is a separate process and must be able to write here, so the work
    # folder lives somewhere ordinary rather than the system temp folder
    work_root = Path(workdir).expanduser().resolve() if workdir else Path.cwd() / ".ligand-docking-tmp"
    work_root.mkdir(parents=True, exist_ok=True)
    projects_root = Path(projects).expanduser().resolve() if projects else Path.cwd() / "projects"
    projects_root.mkdir(parents=True, exist_ok=True)

    DockingHandler.page = page
    DockingHandler.work_root = work_root
    DockingHandler.projects_root = projects_root
    DockingHandler.vina_path = vina

    server = ThreadingHTTPServer((host, port), DockingHandler)
    url = f"http://{host}:{port}/"
    print(f"Docking viewer running at {url}")
    print(f"Projects are saved in {projects_root}")
    print("Press Ctrl+C to stop.")
    if open_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        server.server_close()
