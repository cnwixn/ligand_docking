"""Small-molecule docking helpers: run AutoDock Vina and view the result in 3D."""

from .dock import VinaNotFound, find_vina, platform_binary_names, run_vina, vina_version
from .pdbqt import Pose, parse_vina_out, pdbqt_to_pdb, read_poses, rmsd_to_reference
from .prep import PrepError, meeko_available, prepare_ligand, prepare_receptor
from .report import build_viewer, render_html, viewer_payload
from .serve import run_server

__all__ = [
    "Pose",
    "PrepError",
    "VinaNotFound",
    "build_viewer",
    "find_vina",
    "meeko_available",
    "parse_vina_out",
    "pdbqt_to_pdb",
    "platform_binary_names",
    "prepare_ligand",
    "prepare_receptor",
    "read_poses",
    "render_html",
    "rmsd_to_reference",
    "run_server",
    "run_vina",
    "viewer_payload",
    "vina_version",
]
__version__ = "0.2.0"
