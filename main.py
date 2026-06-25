"""
Entry point for the Three Rivers Metagenomics analysis.

The full pipeline lives in analyze.py; this just runs it so that either
`python main.py` or `python analyze.py` does the same thing.
"""
from analyze import main as run_analysis


def main() -> None:
    run_analysis()


if __name__ == "__main__":
    main()
