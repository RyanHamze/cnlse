# Contributing

Thank you for helping improve `cnlse`.

## Development Setup

```bash
python -m pip install -e ".[dev,app]"
python -m pytest
```

## Guidelines

- Keep public APIs small and documented.
- Add tests for numerical helpers and solver behavior.
- Preserve attribution to the original 2009 MATLAB project when porting logic.
- Prefer clear educational names over terse research-code abbreviations unless the abbreviation is standard in fiber optics.
- Keep examples reproducible.

## Areas To Help

- Adaptive SSFM.
- CNLSE and PMD/polarization propagation.
- More modulation formats.
- Documentation and tutorials.
- Validation against known optical communication examples.
