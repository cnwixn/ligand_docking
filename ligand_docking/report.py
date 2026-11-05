"""Build a standalone, offline NGL.js page for a docking result."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from .pdbqt import Pose, parse_vina_out, pdbqt_to_pdb, rmsd_to_reference

ASSETS = Path(__file__).resolve().parent / "viewer"
ASSET_FILES = ("ngl.js", "parser.js")


def _structure(path: Path) -> tuple[str, str]:
    """Return (text, NGL extension) for a structure file."""
    text = path.read_text(encoding="utf-8", errors="ignore")
    suffix = path.suffix.lower()
    if suffix == ".pdbqt":
        return pdbqt_to_pdb(text), "pdb"
    if suffix in (".cif", ".mmcif"):
        return text, "cif"
    if suffix == ".mol2":
        return text, "mol2"
    if suffix in (".sdf", ".mol"):
        return text, "sdf"
    return text, "pdb"


def _load_poses(poses: str | Path | list[Pose]) -> list[Pose]:
    if isinstance(poses, list):
        return poses
    path = Path(poses).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    parsed = parse_vina_out(path.read_text(encoding="utf-8", errors="ignore"))
    if parsed:
        return parsed
    # a plain ligand file: one pose, no scores
    text, _ = _structure(path)
    return [Pose(mode=1, affinity=None, rmsd_lb=None, rmsd_ub=None, pdb=text)]


def viewer_payload(
    receptor: str | Path,
    poses: str | Path | list[Pose],
    *,
    title: str = "",
    subtitle: str = "",
    crystal: str | Path | None = None,
    receptor_label: str | None = None,
    ligand_label: str | None = None,
    server_mode: bool = False,
) -> dict:
    """Everything the page needs, as plain data (used by the builder and the server)."""
    receptor_path = Path(receptor).expanduser().resolve() if isinstance(receptor, (str, Path)) else None
    if receptor_path is not None:
        if not receptor_path.is_file():
            raise FileNotFoundError(receptor_path)
        receptor_text, receptor_ext = _structure(receptor_path)
        receptor_name = receptor_path.name
    else:
        raise TypeError("receptor must be a path to a structure file")

    pose_list = _load_poses(poses)
    if not pose_list:
        raise ValueError("No poses found in the given file.")

    if not ligand_label:
        if isinstance(poses, (str, Path)):
            ligand_label = Path(poses).name
        else:
            ligand_label = "ligand"

    payload: dict = {
        "title": title or f"{Path(ligand_label).stem} docked into {receptor_path.stem}",
        "subtitle": subtitle,
        "receptor": {"text": receptor_text, "ext": receptor_ext},
        "receptorLabel": receptor_label or receptor_name,
        "ligandLabel": ligand_label,
        "poses": [{**pose.as_dict(), "ext": "pdb"} for pose in pose_list],
        "crystal": None,
        "crystalRmsd": None,
        "serverMode": bool(server_mode),
    }

    if crystal:
        crystal_path = Path(crystal).expanduser().resolve()
        crystal_text, crystal_ext = _structure(crystal_path)
        payload["crystal"] = {"text": crystal_text, "ext": crystal_ext}
        payload["crystalRmsd"] = rmsd_to_reference(pose_list[0].pdb, crystal_text)

    return payload


def render_html(payload: dict) -> str:
    template = (ASSETS / "index.html").read_text(encoding="utf-8")
    return template.replace("/*__DATA__*/{}", json.dumps(payload))


def build_viewer(
    receptor: str | Path,
    poses: str | Path | list[Pose],
    out_html: str | Path,
    *,
    title: str = "",
    subtitle: str = "",
    crystal: str | Path | None = None,
    receptor_label: str | None = None,
) -> Path:
    """Write ``out_html`` plus the two viewer assets next to it.

    ``receptor`` is a structure file (.pdb, .pdbqt, .cif, .mol2).
    ``poses`` is a Vina output file, a single ligand file, or a list of Pose.
    ``crystal`` is an optional reference ligand used for the RMSD badge and the
    grey overlay.
    """
    payload = viewer_payload(
        receptor, poses, title=title, subtitle=subtitle, crystal=crystal,
        receptor_label=receptor_label,
    )

    target = Path(out_html).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    for name in ASSET_FILES:
        shutil.copy2(ASSETS / name, target.parent / name)
    target.write_text(render_html(payload), encoding="utf-8")
    return target
