#!/usr/bin/env python3
"""
Quick checker for ROOT ntuples: list trees and entry counts for each file.
Default directory is /group/belle/users/amubarak/03-KEKCC.
"""

import argparse
import sys
from pathlib import Path

try:
    import uproot
except ImportError:
    print(
        "Missing dependency: uproot. Activate the Belle II python environment "
        "before running this script.",
        file=sys.stderr,
    )
    sys.exit(1)

DEFAULT_DIR = Path("/group/belle/users/amubarak/03-KEKCC")


def summarize_file(path: Path):
    """Return (tree_name -> entries) for a ROOT file, or an error string."""
    try:
        with uproot.open(path) as f:
            out = {}
            for name, obj in f.items():
                # uproot gives names like 'DstreeCh1;1'; strip cycle numbers
                clean_name = name.split(";")[0]
                if hasattr(obj, "num_entries"):
                    out[clean_name] = obj.num_entries
            return out, None
    except Exception as exc:  # pragma: no cover - diagnostics only
        return {}, str(exc)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dir",
        type=Path,
        default=DEFAULT_DIR,
        help="Directory to scan (default: %(default)s)",
    )
    parser.add_argument(
        "--pattern",
        default="*.root",
        help="Glob pattern to match files (default: %(default)s)",
    )
    args = parser.parse_args()

    if not args.dir.exists():
        print(f"Directory not found: {args.dir}", file=sys.stderr)
        sys.exit(1)

    files = sorted(args.dir.glob(args.pattern))
    if not files:
        print(f"No files found under {args.dir} matching {args.pattern}")
        return

    print(f"Scanning {len(files)} ROOT files under {args.dir}")
    for path in files:
        tree_info, err = summarize_file(path)
        print(path)
        if err:
            print(f"  ERROR: {err}")
            continue
        if not tree_info:
            print("  No TTrees found")
            continue
        for tname, n in sorted(tree_info.items()):
            print(f"  {tname}: {n}")
        print()


if __name__ == "__main__":
    main()
