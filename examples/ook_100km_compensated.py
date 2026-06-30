"""Run the first OOK fiber-compensation example."""

from pathlib import Path

from cnlse.examples import run_ook_100km_compensated


if __name__ == "__main__":
    output = run_ook_100km_compensated(Path(__file__).resolve().parents[1] / "outputs" / "ook_100km_compensated.png")
    print(f"Wrote {output}")
