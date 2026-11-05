"""Parse AutoDock PDBQT files and Vina output files.

Vina and AutoDock-GPU write PDBQT, which the browser viewer cannot read
directly, so the atom lines are converted to plain PDB here. Everything in this
module is pure standard library.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# AutoDock atom types -> element symbols
ELEMENT_BY_TYPE = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O", "O": "O",
    "SA": "S", "S": "S", "HD": "H", "H": "H", "HS": "H", "P": "P", "F": "F",
    "CL": "Cl", "BR": "Br", "I": "I", "SI": "Si", "MG": "Mg", "ZN": "Zn",
    "FE": "Fe", "CA": "Ca", "MN": "Mn", "SE": "Se",
}
_KNOWN_ELEMENTS = ["Cl", "Br", "Si", "Na", "Mg", "Zn", "Fe", "Mn", "Ca", "Se",
                   "F", "I", "P", "S", "O", "N", "C", "H"]


@dataclass
class Pose:
    """One docking mode."""

    mode: int
    affinity: float | None
    rmsd_lb: float | None
    rmsd_ub: float | None
    pdb: str

    def as_dict(self) -> dict:
        return {"mode": self.mode, "affinity": self.affinity,
                "rmsd_lb": self.rmsd_lb, "rmsd_ub": self.rmsd_ub, "pdb": self.pdb}


def element_of(atom_type: str, atom_name: str) -> str:
    key = (atom_type or "").strip().upper()
    if key in ELEMENT_BY_TYPE:
        return ELEMENT_BY_TYPE[key]
    for symbol in _KNOWN_ELEMENTS:
        if key.startswith(symbol.upper()):
            return symbol
    first = (atom_name or "").strip()[:1].upper()
    return first or "C"


def atom_line_to_pdb(line: str) -> str:
    """Rewrite one PDBQT atom line as a PDB atom line with an element column."""
    record = line[:6]
    atom_name = line[12:16]
    atom_type = line[77:79] if len(line) >= 79 else line[77:]
    element = element_of(atom_type, atom_name)
    body = line[:54].ljust(54)
    occupancy = line[54:60] if len(line) >= 60 and line[54:60].strip() else "  1.00"
    bfactor = line[60:66] if len(line) >= 66 and line[60:66].strip() else "  0.00"
    head = record + body[6:54] + occupancy + bfactor
    label = element if len(element) == 1 else element[0].upper() + element[1:].lower()
    return head + " " * 10 + label.rjust(2)


def pdbqt_to_pdb(text: str) -> str:
    """Convert receptor or ligand PDBQT text to PDB (coordinates unchanged)."""
    lines = [atom_line_to_pdb(line) for line in text.splitlines()
             if line.startswith(("ATOM", "HETATM"))]
    return "\n".join(lines) + "\nEND\n"


def parse_vina_out(text: str) -> list[Pose]:
    """Read the poses and scores out of a Vina ``--out`` file."""
    poses: list[Pose] = []
    atoms: list[str] = []
    scores: dict[str, float] = {}
    started = False

    def flush() -> None:
        if atoms:
            poses.append(Pose(
                mode=len(poses) + 1,
                affinity=scores.get("affinity"),
                rmsd_lb=scores.get("rmsd_lb"),
                rmsd_ub=scores.get("rmsd_ub"),
                pdb="\n".join(atoms) + "\nEND\n",
            ))
        atoms.clear()

    for line in text.splitlines():
        if line.startswith("MODEL"):
            flush()
            scores = {}
            started = True
        elif line.startswith("REMARK VINA RESULT:"):
            fields = line.split()
            scores = {"affinity": float(fields[3]), "rmsd_lb": float(fields[4]),
                      "rmsd_ub": float(fields[5])}
        elif line.startswith(("ATOM", "HETATM")):
            atoms.append(atom_line_to_pdb(line))
        elif line.startswith("ENDMDL"):
            flush()
            scores = {}
            started = False
    flush()
    return poses


def read_poses(path: str | Path) -> list[Pose]:
    return parse_vina_out(Path(path).read_text(encoding="utf-8", errors="ignore"))


def coordinates(pdb_text: str) -> list[tuple[float, float, float]]:
    points = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            points.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return points


def rmsd_to_reference(pose_pdb: str, reference_pdb: str) -> float | None:
    """Heavy-atom RMSD after optimal superposition, or None if it cannot be computed."""
    first, second = coordinates(pose_pdb), coordinates(reference_pdb)
    count = min(len(first), len(second))
    if count < 3:
        return None
    try:
        import numpy as np
    except ImportError:  # translation only, still a useful rough number
        from math import sqrt
        ca = [sum(p[i] for p in first[:count]) / count for i in range(3)]
        cb = [sum(p[i] for p in second[:count]) / count for i in range(3)]
        total = sum(sum((first[k][i] - ca[i] - (second[k][i] - cb[i])) ** 2 for i in range(3))
                    for k in range(count))
        return round(sqrt(total / count), 2)

    a = np.asarray(first[:count], dtype=float)
    b = np.asarray(second[:count], dtype=float)
    a -= a.mean(axis=0)
    b -= b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    correction = np.diag([1.0, 1.0, 1.0 if np.linalg.det(u @ vt) > 0 else -1.0])
    diff = a @ (u @ correction @ vt) - b
    return round(float(np.sqrt((diff ** 2).sum() / count)), 2)
