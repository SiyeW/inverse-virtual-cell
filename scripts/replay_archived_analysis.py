"""Explicit historical adapter; NOT the adopted analysis entry point.

For exact reproduction, use the run's recorded Git commit in a separate worktree.
This adapter preserves imports when inspecting the archived source at current HEAD.
"""
import argparse
import runpy
import sys
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    archive = root/'archive/20261023_legacy_chip_analysis/scripts'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('script', choices=sorted(p.name for p in archive.glob('*.py')))
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    remaining = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
    sys.path.insert(0, str(archive))
    sys.argv = [str(archive/args.script), *remaining]
    print('HISTORICAL ROUTE: not a validated default or functional-hit pipeline.', file=sys.stderr)
    runpy.run_path(str(archive/args.script), run_name='__main__')


if __name__ == '__main__':
    main()
