"""Unified command-line interface for RGAnnotator.

Usage::

    rgannotator identify    -i <proteome.fa | pep_dir> [-o outdir]      # module 1
    rgannotator rhythm      -i expr.csv -o result.csv [-m jtk|cosinor]  # module 2
    rgannotator pseudotime  -i expr.csv -o pseudotime.png               # module 3
    rgannotator radar       -i expr.csv -o radar.svg [-g CCA1,LHY]
    rgannotator line        -i expr.csv -o line.svg [-g CCA1,LHY]

Every subcommand can also run standalone, e.g. ``python -m rgannotator.rhythm``.
"""

import argparse

try:
    from . import __version__
    from . import identify, rhythm, pseudotime, radar, lineplot
except ImportError:  # allow running as a plain script
    import __init__ as _pkg
    __version__ = _pkg.__version__
    import identify, rhythm, pseudotime, radar, lineplot


# command name -> (parser factory, runner)
COMMANDS = {
    "identify": (identify.build_parser, identify.run_identify),
    "rhythm": (rhythm.build_parser, rhythm.run_rhythm),
    "pseudotime": (pseudotime.build_parser, pseudotime.run_pseudotime),
    "radar": (radar.build_parser, radar.run_radar),
    "line": (lineplot.build_parser, lineplot.run_lineplot),
}

SUMMARY = {
    "identify": "identify circadian clock genes via reciprocal BLASTP (module 1)",
    "rhythm": "detect rhythmic expression with JTK_CYCLE or Cosinor (module 2)",
    "pseudotime": "circadian pseudotime expression plot, CT/ZT as pseudotime (module 3)",
    "radar": "24h circadian radar/polar plot (module 3)",
    "line": "time-course line plot with fitted curves (module 3)",
}


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="rgannotator",
        description=("RGAnnotator: identification, rhythm detection and visualization "
                     "of circadian clock genes."),
    )
    parser.add_argument("-V", "--version", action="version",
                        version="RGAnnotator " + __version__)

    sub = parser.add_subparsers(dest="command", metavar="command")
    for name, (build, _) in COMMANDS.items():
        sub.add_parser(name, parents=[build(add_help=False)], help=SUMMARY[name], description=SUMMARY[name])

    args = parser.parse_args(argv)

    if not args.command:
        parser.print_help()
        raise SystemExit(2)

    COMMANDS[args.command][1](args)


if __name__ == "__main__":
    main()
