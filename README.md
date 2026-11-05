# ligand_docking

Dock a small molecule into a protein with AutoDock Vina, then look at the result
in 3D: the ligand on its own, the receptor with the search box drawn on it, and
the complex, with the score table beside them. Runs are kept as projects, so an
old docking can be reopened weeks later.

No server has to be installed anywhere. The viewer is one HTML file plus two
JavaScript files, and it opens straight from disk. AutoDock Vina stays the
docking engine — this package prepares the inputs, runs the official binary and
reads what it writes; it does not re-implement scoring.

![Docking result viewer](preview.png)

## Install and prepare

Two halves: the Python package, and the Vina binary it calls. Both are needed
before the first run.

### The Python package

```
pip install -e .
```

The viewer needs nothing beyond the standard library. Two optional pieces:

| Extra | Command | What it adds |
| --- | --- | --- |
| `numpy` | `pip install numpy` | a properly superposed RMSD to a reference pose; without it a translation-only estimate is used |
| `prep` | `pip install -e .[prep]` | Meeko, RDKit and gemmi, which convert `.pdb`, `.sdf`, `.mol` and `.mol2` into PDBQT |

`pip install -e .[prep]` is what a normal user wants: give the page a protein
`.pdb` and a ligand `.sdf` and it prepares both for you.

### The Vina binary

AutoDock Vina 1.2.7 is bundled, so nothing extra is downloaded:

| File | Platform |
| --- | --- |
| `bin/vina_win.exe` | Windows x86_64 |
| `bin/vina_linux_x86_64` | Linux x86_64 |
| `bin/vina_mac_x86_64` | macOS Intel |
| `bin/vina_mac_aarch64` | macOS Apple silicon |

The package picks the right one for the machine and sets the executable bit on
POSIX if a checkout lost it. `ligand-docking vina` prints which binary will be
used. A build of your own wins over the bundled files: pass `--vina /path/to/vina`
or set `VINA_EXE`. On Linux aarch64, drop a `bin/vina_linux_aarch64` from the same
release next to the others.

The binaries are unmodified files from the
[1.2.7 release](https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7),
redistributed under the Apache License 2.0 — see `bin/LICENSE-AutoDock-Vina`. They
are third-party executables: on macOS, run `xattr -d com.apple.quarantine bin/vina_mac_*`
if Gatekeeper refuses to start them.

### Which files you can give it

The Run panel and `ligand-docking dock` accept the usual formats and convert them
when Meeko is installed:

| Input | Receptors | Ligands |
| --- | --- | --- |
| prepared | `.pdbqt` | `.pdbqt`, a Vina `--out` file |
| converted for you | `.pdb`, `.ent`, `.cif`, `.mmcif` | `.sdf`, `.mol`, `.mol2`, `.smi`/`.smiles`, small-molecule `.pdb` |

mmCIF is read with gemmi and turned into PDB first. Ligands given as SMILES get 3D
coordinates from RDKit (ETKDG embedding, MMFF refinement); SDF, MOL and PDB inputs
get explicit hydrogens added, which Meeko requires. Receptor residues whose
backbone is incomplete are dropped before conversion, and the run reports which
ones. Without Meeko the page says so instead of failing silently.

### A first smoke test

Two examples ship with the repository, both taken from the official Vina set:

| Folder | System | Best score |
| --- | --- | --- |
| `examples/1iep/` | imatinib in Abl kinase (PDB 1IEP) | −13.22 kcal/mol |
| `examples/1s63/` | inhibitor in a zinc metalloprotein (PDB 1S63) | −9.36 kcal/mol |

Each folder holds the receptor as `.pdb` and `.pdbqt`, the ligand, the Vina output
and the box the run used:

```
ligand-docking view --receptor examples/1iep/receptor.pdb \
    --poses examples/1iep/vina_out.pdbqt --crystal examples/1iep/ligand.pdbqt \
    -o out/1iep.html --open

ligand-docking view --receptor examples/1s63/receptor.pdb \
    --poses examples/1s63/vina_out.pdbqt --crystal examples/1s63/ligand.pdbqt \
    -o out/1s63.html --open
```

## Run a docking from the browser

A static page cannot start a program, so `serve` starts a small local server on
the loopback interface and the page grows a **Run docking** panel:

```
ligand-docking serve --receptor examples/1iep/receptor.pdb \
    --poses examples/1iep/vina_out.pdbqt --open
```

Then: name the project, choose a receptor and a ligand, press **Run Vina**. The
chosen files appear in the views immediately, PDB/SDF are converted to PDBQT, and
when the run finishes the poses and the score table are replaced. Nothing is
uploaded anywhere: files go from your browser to the local server on 127.0.0.1.

The box is filled in for you from a pose that already sits in the pocket (the
reference pose if there is one, otherwise the first docked pose):

| Field | Meaning |
| --- | --- |
| **centre** | centroid of that pose, in Å — where the box sits |
| **size** | its extent plus 12 Å, never below 18 Å — how far Vina may search |
| **effort** | Vina exhaustiveness: how many independent searches. 8 is a quick look, 32–64 a result worth keeping. It does not set how many poses come back |
| **poses** | the cap on returned poses (Vina's `--num-modes`, 9). Fewer are written when the searches only find a few distinct clusters |

The search box is drawn in the receptor view, so you can see what is being
sampled; change centre or size and the box moves with you. The views spin by
default — untick **Spin all views** in the Display panel to stop them.

Choosing a new receptor or a new ligand clears the complex: the old pose belongs
to the old pair, so it is removed rather than left floating in the new structure.
Run again to fill it.

### Projects

Every press of **Run Vina** becomes a project under `projects/`:

```
projects/20261105-173506-809add-imatinib-test/
├── receptor.pdbqt     prepared receptor
├── ligand.pdbqt       prepared ligand
├── poses.pdbqt        the Vina output
└── project.json       name, time, box, effort and the score of every pose
```

The **Projects** list in the page shows what has been run; click one to reopen its
structures, poses, scores and box. The list is read from disk on every request, so
projects survive restarts and can be copied or archived like any other folder.
They stay local — `projects/` is in `.gitignore`.

## Command line

| Command | What it does |
| --- | --- |
| `ligand-docking dock` | run Vina on prepared inputs, optionally write the viewer |
| `ligand-docking view` | build a viewer page from an existing run |
| `ligand-docking serve` | serve the viewer with the Run panel and projects |
| `ligand-docking vina` | print the Vina binary that will be used |

```
ligand-docking dock --receptor receptor.pdbqt --ligand ligand.pdbqt \
    --config receptor.box.txt --exhaustiveness 32 \
    -o run_out.pdbqt --html run.html --open

ligand-docking view --receptor receptor.pdb --poses run_out.pdbqt \
    --crystal ligand.pdbqt -o run.html --open
```

`dock` options: `--receptor` and `--ligand` (required), `-o/--out`,
`--config` or `--center X Y Z` with `--size X Y Z`, `--exhaustiveness` (32),
`--num-modes` (9), `--seed` (42), `--cpu` (0 = every core), `--vina`,
`--html`, `--crystal`, `--title`, `--open`.

`view` options: `--receptor`, `--poses`, `-o/--out`, `--crystal`, `--title`,
`--subtitle`, `--open`. `--poses` takes a Vina `--out` file (scores and every
pose) or a single ligand file (one pose, empty score table).

`serve` options: `--receptor` and `--poses` (what the page shows first),
`--crystal`, `--title`, `--subtitle`, `--port` (8766), `--host` (127.0.0.1),
`--workdir` (`.ligand-docking-tmp`), `--projects` (`./projects`), `--vina`,
`--open`. The work folder must be writable by the Vina process, which is why it
defaults to the current directory rather than the system temp folder.

## What the page shows

- **Ligand** — the chosen or docked pose on its own, ball and stick by element.
- **Receptor** — cartoon coloured by secondary structure, with the search box drawn as a wireframe.
- **Complex** — the pose inside the pocket.
- **Docking modes** — affinity and RMSD l.b./u.b. per pose; click a row to switch
  every view to that pose.
- **Display** — receptor surface, the reference pose in grey, and spin.

The title is built from the file names unless `--title` says otherwise; the line
under it names the Vina build that produced the poses.

## Using it from Python

```python
from ligand_docking import run_vina, build_viewer, prepare_receptor, prepare_ligand

receptor = prepare_receptor("protein.pdb", "work/receptor.pdbqt")["pdbqt"]
ligand = prepare_ligand("ligand.sdf", "work/ligand.pdbqt")["pdbqt"]
result = run_vina(receptor, ligand, "work/out.pdbqt",
                  center=(15.19, 53.90, 16.92), size=(20, 20, 20),
                  exhaustiveness=32)
print(result["poses"][0].affinity)
build_viewer(receptor, result["poses"], "run.html", crystal="ligand.pdbqt")
```

## Limits worth knowing

- Vina scores rank poses; they are not binding free energies. Do not read a
  −10 kcal/mol result as a measured affinity.
- A low `effort` returns fewer distinct poses, not a cheaper version of the same
  answer; raise it before believing a ranking.
- The viewer needs a WebGL-capable browser. Without the local server it is
  read-only: no Run panel, no projects.
- Nothing here prepares a protein beyond atom typing: missing loops, alternate
  locations and protonation states remain your call.

## Licence

MIT, see [LICENSE](LICENSE). The bundled Vina binaries are Apache 2.0, see
`bin/LICENSE-AutoDock-Vina`.
