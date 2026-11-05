"""Command line interface: ``ligand-docking dock`` and ``ligand-docking view``."""

from __future__ import annotations

import argparse
import sys
import webbrowser
from pathlib import Path

from .dock import VinaNotFound, find_vina, run_vina, vina_version
from .pdbqt import read_poses
from .report import build_viewer, viewer_payload
from .serve import run_server


def default_subtitle() -> str:
    """The line under the title: which Vina build produced the poses."""
    version = vina_version()
    return f"AutoDock Vina {version}" if version else "AutoDock Vina"


def _print_summary(result: dict) -> None:
    poses = result["poses"]
    print(f"Poses : {len(poses)}")
    print(f"Output: {result['out']}")
    if poses:
        print("  mode  affinity  rmsd_lb  rmsd_ub")
        for pose in poses:
            print(f"  {pose.mode:>4}  {pose.affinity:>8.2f}  {pose.rmsd_lb:>7.2f}  {pose.rmsd_ub:>7.2f}")


def cmd_dock(args: argparse.Namespace) -> int:
    try:
        result = run_vina(
            args.receptor,
            args.ligand,
            args.out,
            center=tuple(args.center) if args.center else None,
            size=tuple(args.size) if args.size else None,
            config=args.config,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            seed=args.seed,
            cpu=args.cpu,
            vina=args.vina,
        )
    except VinaNotFound as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    _print_summary(result)

    if args.html:
        viewer = build_viewer(
            args.receptor,
            result["poses"],
            args.html,
            title=args.title or "",
            subtitle="AutoDock Vina docking result",
            crystal=args.crystal,
        )
        print(f"Viewer: {viewer}")
        if args.open:
            webbrowser.open(viewer.as_uri())
    return 0


def cmd_view(args: argparse.Namespace) -> int:
    try:
        viewer = build_viewer(
            args.receptor,
            args.poses,
            args.out,
            title=args.title or "",
            subtitle=args.subtitle or default_subtitle(),
            crystal=args.crystal,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    poses = read_poses(args.poses) if Path(args.poses).is_file() else []
    if poses:
        best = poses[0]
        print(f"Poses : {len(poses)} (best {best.affinity:.2f} kcal/mol)")
    print(f"Viewer: {viewer}")
    print("Open it in a browser; ngl.js and parser.js must sit beside it.")
    if args.open:
        webbrowser.open(viewer.as_uri())
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ligand-docking",
        description="Dock a small molecule into a protein with AutoDock Vina and view the result in 3D.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    dock = sub.add_parser("dock", help="run AutoDock Vina on prepared PDBQT files")
    dock.add_argument("--receptor", required=True, help="receptor PDBQT")
    dock.add_argument("--ligand", required=True, help="ligand PDBQT")
    dock.add_argument("-o", "--out", help="Vina output PDBQT (default: <ligand>_out.pdbqt)")
    box = dock.add_argument_group("search box (config file or centre and size)")
    box.add_argument("--config", help="Vina config file with centre_x/size_x ...")
    box.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"))
    box.add_argument("--size", nargs=3, type=float, metavar=("X", "Y", "Z"))
    dock.add_argument("--exhaustiveness", type=int, default=32)
    dock.add_argument("--num-modes", type=int, default=9)
    dock.add_argument("--seed", type=int, default=42)
    dock.add_argument("--cpu", type=int, default=0, help="0 uses every core")
    dock.add_argument("--vina", help="path to the Vina executable")
    dock.add_argument("--html", help="also write a 3D viewer page here")
    dock.add_argument("--crystal", help="reference ligand for the RMSD badge")
    dock.add_argument("--title", help="title shown in the viewer")
    dock.add_argument("--open", action="store_true", help="open the viewer when it is ready")
    dock.set_defaults(func=cmd_dock)

    view = sub.add_parser("view", help="build a static 3D viewer page from docking results")
    view.add_argument("--receptor", required=True, help="receptor .pdb / .pdbqt / .cif / .mol2")
    view.add_argument("--poses", required=True, help="Vina output PDBQT, or a single ligand file")
    view.add_argument("-o", "--out", default="viewer.html", help="HTML file to write")
    view.add_argument("--crystal", help="reference ligand for the RMSD badge and grey overlay")
    view.add_argument("--title", help="page title")
    view.add_argument("--subtitle", help="line under the title")
    view.add_argument("--open", action="store_true", help="open the page when it is ready")
    view.set_defaults(func=cmd_view)

    serve = sub.add_parser("serve", help="serve the viewer locally so its Run button can dock")
    serve.add_argument("--receptor", required=True, help="receptor .pdb / .pdbqt")
    serve.add_argument("--poses", required=True, help="poses to show first (Vina out or a ligand file)")
    serve.add_argument("--crystal", help="reference ligand for the RMSD badge")
    serve.add_argument("--title", help="page title (default: built from the file names)")
    serve.add_argument("--subtitle", help="line under the title")
    serve.add_argument("--port", type=int, default=8766)
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--vina", help="path to the Vina executable")
    serve.add_argument("--workdir", help="where run folders are written (default: ./.ligand-docking-tmp)")
    serve.add_argument("--open", action="store_true", help="open the page when the server is up")
    serve.set_defaults(func=cmd_serve)

    version = sub.add_parser("vina", help="print the Vina executable that will be used")
    version.add_argument("--vina", help="path to check instead of the automatic search")
    version.set_defaults(func=cmd_vina)

    return parser


def cmd_serve(args: argparse.Namespace) -> int:
    try:
        page = viewer_payload(
            args.receptor,
            args.poses,
            title=args.title or "",
            subtitle=args.subtitle or default_subtitle(),
            crystal=args.crystal,
            server_mode=True,
        )
    except (FileNotFoundError, ValueError, TypeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    run_server(page, host=args.host, port=args.port, vina=args.vina,
               open_browser=args.open, workdir=args.workdir)
    return 0


def cmd_vina(args: argparse.Namespace) -> int:
    try:
        print(find_vina(args.vina))
    except VinaNotFound as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
