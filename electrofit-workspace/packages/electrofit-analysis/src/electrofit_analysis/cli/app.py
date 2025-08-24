"""
Unified CLI entry point for electrofit-analysis.

Provides subcommands that forward to existing CLI modules in this package.

Subcommands
-----------
- distance:     Compute Na–P distances for IP6 microstates.
- count:        Compute Na–P distances, Na+ counts, buried volume, and excess.
- summarize-nap:    Summarize distances/counts/coordination across patterns (requires out dir).
- plot-2d:      Render 2D molecule depictions for microstates.
- coordination: Analyze Na–IP6 coordination (counts, RDFs, optional projections).

Usage examples
--------------
python -m electrofit_analysis.cli.app distance -p /path/to/project
python -m electrofit_analysis.cli.app count -p /path/to/project
python -m electrofit_analysis.cli.app summarize-nap -p /path/to/project -o /path/to/out
python -m electrofit_analysis.cli.app plot-2d -p /path/to/project [--subdir process]
python -m electrofit_analysis.cli.app coordination -p /path/to/project [--subdir process] [--determine-global-y] [--rdf-y-max 1800] [--plot-projection]
"""

from __future__ import annotations

import argparse
import os


def _cmd_distance(args: argparse.Namespace) -> None:
    from .coordination.na_p_distance_ip6 import main as distance_main

    project = os.path.abspath(args.project)
    distance_main(project)


def _cmd_count(args: argparse.Namespace) -> None:
    from .coordination.na_p_dist_count_ip6 import main as count_main

    project = os.path.abspath(args.project)
    count_main(project)


def _cmd_summarize_nap(args: argparse.Namespace) -> None:
    from .coordination.summerize_nap_dist_count_ip6 import main as summarize_main

    project = os.path.abspath(args.project)
    root_dir = os.path.join(project, "process")
    out_dir = os.path.abspath(args.output_dir)
    summarize_main(root_dir, out_dir)


def _cmd_plot_2d(args: argparse.Namespace) -> None:
    from .plot_molecule_2d import main as plot2d_main

    project = os.path.abspath(args.project)
    plot2d_main(project, args.subdir)


def _cmd_coordination(args: argparse.Namespace) -> None:
    from .coordination.Na_IP6_coordination import main as coord_main

    project = os.path.abspath(args.project)
    coord_main(
        project_dir=project,
        subdir=args.subdir,
        determine_global_y=args.determine_global_y,
        rdf_y_max=args.rdf_y_max,
        plot_projection=args.plot_projection,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="electrofit-analysis",
        description="Unified CLI for electrofit-analysis tools.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # distance
    p_dist = sub.add_parser(
        "distance",
        help="Compute Na–P distances for an IP6 project (per microstate).",
    )
    p_dist.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_dist.set_defaults(func=_cmd_distance)

    # count
    p_count = sub.add_parser(
        "count",
        help=(
            "Run Na–P distances, Na+ counts, buried volume, and excess analysis for an IP6 project."
        ),
    )
    p_count.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_count.set_defaults(func=_cmd_count)

    # summarize-nap
    p_sum = sub.add_parser(
        "summarize-nap",
        help=(
            "Summarize Na–P distances/counts/coordination across patterns. Requires an output directory."
        ),
    )
    p_sum.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_sum.add_argument(
        "-o", "--out", "--output", dest="output_dir", required=True, help="Output directory."
    )
    p_sum.set_defaults(func=_cmd_summarize_nap)

    # plot-2d
    p_plot = sub.add_parser(
        "plot-2d", help="Render 2D molecule depictions (SVG) for microstates."
    )
    p_plot.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_plot.add_argument(
        "--subdir",
        default="process",
        help="Relative subdirectory under the project root to traverse (default: process).",
    )
    p_plot.set_defaults(func=_cmd_plot_2d)

    # coordination
    p_coord = sub.add_parser(
        "coordination",
        help=(
            "Analyze Na–IP6 coordination across microstates (counts, RDFs, optional projections)."
        ),
    )
    p_coord.add_argument(
        "-p", "--project", required=True, help="Path to the project root directory."
    )
    p_coord.add_argument(
        "--subdir",
        default="process",
        help="Relative subdirectory under the project root to traverse (default: process).",
    )
    p_coord.add_argument(
        "--determine-global-y",
        action="store_true",
        help="Scan all microstates to determine a common RDF y-limit.",
    )
    p_coord.add_argument(
        "--rdf-y-max",
        type=float,
        default=None,
        help="Override RDF y-limit with a fixed value (skips scanning).",
    )
    p_coord.add_argument(
        "--plot-projection",
        action="store_true",
        help="Also generate 2D/3D network projection snapshots.",
    )
    p_coord.set_defaults(func=_cmd_coordination)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
