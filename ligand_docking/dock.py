"""Run AutoDock Vina (or a compatible build) from Python and read its output."""

from __future__ import annotations

import os
import platform
import shutil
import stat
import subprocess
import sys
from pathlib import Path

from .pdbqt import Pose, read_poses

VINA_ENV_VAR = "VINA_EXE"
_PACKAGE_ROOT = Path(__file__).resolve().parent.parent


class VinaNotFound(RuntimeError):
    """Raised when no Vina executable can be located."""


def platform_binary_names() -> list[str]:
    """Names the bundled binary can have on this machine, best first."""
    if sys.platform.startswith("win"):
        return ["vina.exe", "vina_win.exe"]
    machine = platform.machine().lower()
    if sys.platform == "darwin":
        if machine in ("arm64", "aarch64"):
            return ["vina", "vina_mac_aarch64", "vina_mac_x86_64"]
        return ["vina", "vina_mac_x86_64", "vina_mac_aarch64"]
    if machine in ("aarch64", "arm64"):
        return ["vina", "vina_linux_aarch64", "vina_linux_x86_64"]
    return ["vina", "vina_linux_x86_64", "vina_linux_aarch64"]


def _make_executable(path: Path) -> None:
    """Bundled binaries lose their permission bits when checked out on POSIX."""
    if os.name != "posix":
        return
    try:
        mode = path.stat().st_mode
        if not mode & stat.S_IXUSR:
            path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    except OSError:
        pass


def find_vina(explicit: str | Path | None = None) -> Path:
    """Locate the Vina executable.

    Order: the path you pass, ``$VINA_EXE``, a bundled binary under ``bin/``
    (picked for this platform), then whatever is on ``PATH``.
    """
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    if os.environ.get(VINA_ENV_VAR):
        candidates.append(Path(os.environ[VINA_ENV_VAR]).expanduser())
    for base in (Path.cwd(), _PACKAGE_ROOT):
        candidates += [base / "bin" / name for name in platform_binary_names()]
    for on_path in (shutil.which("vina"), shutil.which("vina.exe")):
        if on_path:
            candidates.append(Path(on_path))

    for candidate in candidates:
        if candidate and candidate.is_file():
            _make_executable(candidate)
            return candidate.resolve()
    wanted = ", ".join(platform_binary_names())
    raise VinaNotFound(
        "AutoDock Vina was not found. Put a binary in bin/ (this platform looks "
        f"for {wanted}), download one from "
        "https://github.com/ccsb-scripps/AutoDock-Vina/releases, or set "
        f"{VINA_ENV_VAR} to its path."
    )


def run_vina(
    receptor: str | Path,
    ligand: str | Path,
    out: str | Path | None = None,
    *,
    center: tuple[float, float, float] | None = None,
    size: tuple[float, float, float] | None = None,
    config: str | Path | None = None,
    exhaustiveness: int = 32,
    num_modes: int = 9,
    seed: int | None = 42,
    cpu: int = 0,
    vina: str | Path | None = None,
    extra_args: list[str] | None = None,
) -> dict:
    """Dock ``ligand`` into ``receptor`` and return the parsed poses.

    Either give ``config`` (a Vina config file with the box) or ``center`` and
    ``size``. Both structures must already be PDBQT.
    """
    receptor_path = Path(receptor).expanduser().resolve()
    ligand_path = Path(ligand).expanduser().resolve()
    for path in (receptor_path, ligand_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    out_path = Path(out).expanduser().resolve() if out else ligand_path.with_name(
        f"{ligand_path.stem}_out.pdbqt")
    exe = find_vina(vina)

    command = [str(exe), "--receptor", str(receptor_path), "--ligand", str(ligand_path),
               "--out", str(out_path), "--exhaustiveness", str(exhaustiveness),
               "--num_modes", str(num_modes)]
    if config:
        command += ["--config", str(Path(config).expanduser().resolve())]
    elif center and size:
        command += ["--center_x", str(center[0]), "--center_y", str(center[1]),
                    "--center_z", str(center[2]),
                    "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2])]
    else:
        raise ValueError("Give either config= or both center= and size=.")
    if seed is not None:
        command += ["--seed", str(seed)]
    if cpu is not None:
        command += ["--cpu", str(cpu)]
    if extra_args:
        command += list(extra_args)

    completed = subprocess.run(command, capture_output=True, text=True, timeout=24 * 3600)
    log = (completed.stderr or "") + (completed.stdout or "")
    if completed.returncode != 0:
        raise RuntimeError(f"Vina exited with {completed.returncode}:\n{log.strip()}")

    poses: list[Pose] = read_poses(out_path) if out_path.is_file() else []
    return {"out": out_path, "log": log, "poses": poses, "command": command}
