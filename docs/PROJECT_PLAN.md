# Project Plan

## Batch 1: Foundation

- Create the `cnlse` package scaffold.
- Add README, license, citation, notice, and attribution.
- Add initial examples and placeholder modules.

## Batch 2: MATLAB Port Map

- Map original MATLAB files to Python modules.
- Preserve original authorship notes in ported modules.
- Add tests that compare simple MATLAB-equivalent behavior where practical.

## Batch 3: Core Signal Model

- Port pattern generation.
- Port electric source generation.
- Port laser source generation.
- Port Mach-Zehnder modulation.
- Port average power calculation.

## Batch 4: Fiber Solver

- Port attenuation and effective length calculations.
- Port dispersion phase operator.
- Port scalar SSFM propagation.
- Add nonlinear phase rotation and step-control options.

## Batch 5: CNLSE/Polarization

- Add coupled-field data model.
- Add PMD/birefringence structures.
- Add Manakov-oriented propagation option.

## Batch 6: Examples And Tutorials

- Recreate the original 100 km OOK compensated-fiber simulation.
- Add notebooks and small guides.
- Add labeled plots explaining input pulse, distorted output, and compensated output.

## Batch 7: Friendly Interface

- Add a small local GUI or Streamlit app.
- Add labeled controls and hover explanations.
- Add export buttons for plots and CSV data.

## Batch 8: Release

- Add GitHub Actions tests.
- Publish documentation.
- Tag `v0.1.0`.
