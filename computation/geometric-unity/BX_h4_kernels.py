"""Run the full exact H4 rank, closure, and character checks."""

from pathlib import Path
import subprocess
import sys


DATA = Path(__file__).resolve().parent / "BX_data"


def main():
    for filename, flags in (
        ("verify_h4_an.py", ("--exhaustive",)),
        ("verify_h4_dual_curvature.py", ("--exhaustive",)),
        ("verify_h4_character.py", ()),
    ):
        subprocess.run([sys.executable, "-X", "utf8", str(DATA / filename), *flags], check=True)


if __name__ == "__main__":
    main()
