"""Turn ordinary structure files into the PDBQT files Vina needs.

Meeko (the Forli lab's preparation tool) does the chemistry: it assigns atom
types, Gasteiger charges and rotatable bonds. This module only locates it, cleans
the input, converts between formats and reports clearly when something cannot be
prepared.

    pip install meeko rdkit gemmi
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

LIGAND_SUFFIXES = {".sdf", ".mol", ".mol2"}
SMILES_SUFFIXES = {".smi", ".smiles"}
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


def _cif_to_pdb(source: Path, target: Path) -> None:
    """mmCIF -> PDB, so the rest of the pipeline only sees PDB."""
    try:
        import gemmi
    except ImportError:
        raise PrepError(
            "Reading .cif needs gemmi: pip install gemmi (or give a .pdb)."
        ) from None
    structure = gemmi.read_structure(str(source))
    structure.write_pdb(str(target))


def _write_ligand_sdf(source: Path, target: Path) -> None:
    """Normalise any ligand input to an SDF with explicit hydrogens and 3D coordinates."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise PrepError(
            "Reading this ligand format needs RDKit: pip install rdkit."
        ) from None

    suffix = source.suffix.lower()
    molecules = []

    if suffix in SMILES_SUFFIXES:
        for line in source.read_text(encoding="utf-8", errors="ignore").splitlines():
            text = line.split("#", 1)[0].strip()
            if not text:
                continue
            molecule = Chem.MolFromSmiles(text.split()[0])
            if molecule is None:
                raise PrepError(f"Could not read the SMILES: {text[:60]}")
            molecules.append(molecule)
        if not molecules:
            raise PrepError("The SMILES file is empty.")
    elif suffix == ".pdb":
        molecule = Chem.MolFromPDBFile(str(source), removeHs=False, sanitize=True)
        if molecule is None:
            raise PrepError(f"Could not read a molecule from {source.name}.")
        molecules.append(molecule)
    else:
        supplier = Chem.SDMolSupplier(str(source), removeHs=False)
        molecules = [mol for mol in supplier if mol is not None]
        if not molecules:
            raise PrepError(f"Could not read any molecule from {source.name}.")

    writer = Chem.SDWriter(str(target))
    try:
        for molecule in molecules:
            molecule = Chem.AddHs(molecule, addCoords=molecule.GetNumConformers() > 0)
            if molecule.GetNumConformers() == 0:
                if AllChem.EmbedMolecule(molecule, randomSeed=42) != 0:
                    raise PrepError("Could not generate 3D coordinates for the ligand.")
                AllChem.MMFFOptimizeMolecule(molecule, maxIters=500)
            writer.write(molecule)
    finally:
        writer.close()


def prepare_receptor(source: str | Path, out_pdbqt: str | Path) -> dict:
    """PDB / mmCIF -> receptor PDBQT. Returns the path and any dropped residues."""
    source_path = Path(source).expanduser().resolve()
    target = Path(out_pdbqt).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)

    if source_path.suffix.lower() == ".pdbqt":
        target.write_bytes(source_path.read_bytes())
        return {"pdbqt": target, "dropped": []}

    pdb_path = source_path
    if source_path.suffix.lower() in (".cif", ".mmcif"):
        pdb_path = target.with_suffix(".from_cif.pdb")
        _cif_to_pdb(source_path, pdb_path)

    cleaned_path = target.with_suffix(".cleaned.pdb")
    cleaned, dropped = clean_pdb(pdb_path.read_text(encoding="utf-8", errors="ignore"))
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
    """PDBQT, SDF/MOL/MOL2, SMILES or a small-molecule PDB -> ligand PDBQT."""
    source_path = Path(source).expanduser().resolve()
    target = Path(out_pdbqt).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)

    if source_path.suffix.lower() == ".pdbqt":
        target.write_bytes(source_path.read_bytes())
        return {"pdbqt": target}

    if source_path.suffix.lower() in (LIGAND_SUFFIXES | SMILES_SUFFIXES | {".pdb"}):
        sdf_path = target.with_suffix(".prepared.sdf")
        _write_ligand_sdf(source_path, sdf_path)
    else:
        sdf_path = source_path

    result = _run_meeko("meeko.cli.mk_prepare_ligand",
                        ["-i", str(sdf_path), "-o", str(target)])
    if result.returncode != 0 and "No module named 'meeko'" in (result.stderr or ""):
        raise _meeko_missing(result.stderr or "")
    if result.returncode != 0 or not target.is_file():
        raise PrepError(
            f"Could not prepare the ligand from {source_path.name}: "
            f"{(result.stderr or result.stdout or '').strip()[-400:]}"
        )
    return {"pdbqt": target}
