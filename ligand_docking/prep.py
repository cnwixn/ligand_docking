"""Turn ordinary structure files into the PDBQT files Vina needs.

Meeko (the Forli lab's preparation tool) does the chemistry: it assigns atom
types, Gasteiger charges and rotatable bonds. This module only locates it, feeds
it sensible input and reports clearly when something cannot be prepared.

    pip install meeko rdkit gemmi
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

LIGAND_SUFFIXES = {".sdf", ".mol", ".mol2"}
RECEPTOR_SUFFIXES = {".pdb", ".ent", ".cif", ".mmcif"}
BACKBONE = {"N", "CA", "C", "O"}


class PrepError(RuntimeError):
    """Raised when a structure cannot be turned into PDBQT."""


def meeko_available() -> bool:
    try:
        subprocess.run([sys.executable, "-c", "import meeko"], check=True,
                       capture_output=True, timeout=120)
        return True
    except Exception:  # noqa: BLE001
        return False


def _run_meeko(module: str, args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-m", module, *args],
                          capture_output=True, text=True, timeout=3600)


def _meeko_missing(stderr: str) -> PrepError:
    return PrepError(
        "Preparing structures needs Meeko: pip install meeko rdkit gemmi. "
        f"Meeko said: {(stderr or '').strip()[-300:]}"
    )


def clean_pdb(text: str) -> tuple[str, list[str]]:
    """Drop residues with an incomplete backbone, which strict parsers reject.

    Returns the cleaned PDB text and the labels of the residues removed.
    """
    residues: dict[tuple[str, str, str], list[str]] = {}
    order: list[tuple[str, str, str]] = []
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            key = (line[21:22].strip(), line[22:27].strip(), line[26:27].strip())
            if key not in residues:
                residues[key] = []
                order.append(key)
            residues[key].append(line)

    kept, dropped = [], []
    for key in order:
        lines = residues[key]
        names = {line[12:16].strip() for line in lines if line.startswith("ATOM")}
        if not names:
            continue                    # HETATM only: waters, ions, ligands
        if not BACKBONE.issubset(names):
            chain, number, icode = key
            dropped.append(f"{chain}:{number}{icode}")
            continue
        kept.extend(lines)

    body = "\n".join(kept)
    return f"{body}\nEND\n", dropped


def prepare_receptor(source: str | Path, out_pdbqt: str | Path) -> dict:
    """PDB/mmCIF -> receptor PDBQT. Returns the path and any dropped residues."""
    source_path = Path(source).expanduser().resolve()
    target = Path(out_pdbqt).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)

    if source_path.suffix.lower() == ".pdbqt":
        target.write_bytes(source_path.read_bytes())
        return {"pdbqt": target, "dropped": []}

    cleaned_path = target.with_suffix(".cleaned.pdb")
    text = source_path.read_text(encoding="utf-8", errors="ignore")
    cleaned, dropped = clean_pdb(text)
    if not cleaned.strip():
        raise PrepError(f"{source_path.name} has no complete protein residues.")
    cleaned_path.write_text(cleaned, encoding="utf-8")

    base = target.with_suffix("")
    result = _run_meeko("meeko.cli.mk_prepare_receptor",
                        ["--read_pdb", str(cleaned_path), "-o", str(base), "-p"])
    if result.returncode != 0 and "No module named 'meeko'" in (result.stderr or ""):
        raise _meeko_missing(result.stderr or "")
    if result.returncode != 0 or not target.is_file():
        raise PrepError(
            f"Could not prepare the receptor from {source_path.name}: "
            f"{(result.stderr or result.stdout or '').strip()[-400:]}"
        )
    return {"pdbqt": target, "dropped": dropped}


def prepare_ligand(source: str | Path, out_pdbqt: str | Path) -> dict:
    """SDF/MOL/MOL2 -> ligand PDBQT."""
    source_path = Path(source).expanduser().resolve()
    target = Path(out_pdbqt).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)

    if source_path.suffix.lower() == ".pdbqt":
        target.write_bytes(source_path.read_bytes())
        return {"pdbqt": target}

    result = _run_meeko("meeko.cli.mk_prepare_ligand",
                        ["-i", str(source_path), "-o", str(target)])
    if result.returncode != 0 and "No module named 'meeko'" in (result.stderr or ""):
        raise _meeko_missing(result.stderr or "")
    if result.returncode != 0 or not target.is_file():
        raise PrepError(
            f"Could not prepare the ligand from {source_path.name}: "
            f"{(result.stderr or result.stdout or '').strip()[-400:]}"
        )
    return {"pdbqt": target}
