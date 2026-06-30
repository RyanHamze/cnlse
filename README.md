# Fiber Optics Modulation Using Coupled Nonlinear Schrodinger Equation

`cnlse` is an open-source Python package for learning, simulating, and visualizing optical fiber modulation systems based on the nonlinear Schrodinger equation (NLSE), the coupled nonlinear Schrodinger equation (CNLSE), and split-step Fourier method (SSFM) propagation.

The first release focuses on approachable examples: generating an OOK optical pulse train, propagating it through fiber with attenuation and chromatic dispersion, adding optional nonlinear phase rotation, and comparing uncompensated and dispersion-compensated outputs.

## Origin And Credit

This project is based on the original 2009 senior project:

**Design and Implementation of Fiber Optics Using MATLAB**  
Arts, Sciences, and Technology University in Lebanon  
Beirut, Lebanon, October 2009

Original project authors:

- Ryan Maher Hamze
- Krayem Ahmad Ajjawi
- Meriam Hawwa

Related ebook:

**Design, Implementation and Simulation of Fiber Optics Modulation Using MATLAB**  
ASIN: `B077Z6FMF3`  
2nd edition, published December 6, 2017

The original MATLAB code, report structure, simulation approach, and educational framing are credited to the three original authors above. The Python package implementation and current open-source maintenance are led by **Ryan Maher Hamze**. Contact/reference: [ryanhamze.com](https://ryanhamze.com).

## Planned Features

- Python implementation of the original MATLAB simulation architecture.
- OOK pulse generation and optical source/modulator helpers.
- Fiber propagation using attenuation, GVD, and optional nonlinear phase effects.
- CNLSE/Manakov-oriented extension points for polarization-aware propagation.
- Example simulations with saved plots.
- Beginner guides and tutorials.
- Optional GUI or notebook widgets with clear labels, buttons, and hover explanations.

## Quick Start

Install the package locally and run the first OOK compensation simulation:

```bash
python -m pip install -e .
python examples/ook_100km_compensated.py
```

Expected output:

- Input OOK pulse plot.
- Output after fiber without compensation.
- Output after dispersion-compensating fiber.

## Tutorials

- [Tutorial 1: OOK Over 100 km With Dispersion Compensation](docs/tutorials/01_ook_100km_compensated.md)
- [Tutorial 2: Interactive Streamlit App](docs/tutorials/02_streamlit_app.md)

## Interactive App

Install the app extras:

```bash
python -m pip install -e ".[app]"
```

Run the local tool:

```bash
streamlit run app.py
```

For Streamlit Community Cloud, deploy:

- Repository: `RyanHamze/cnlse`
- Branch: `main`
- Main file path: `streamlit_app.py`

## License

This project is released under the MIT License. See [LICENSE](LICENSE).

## Citation

If this project helps your work, please cite the original thesis project and this Python package. See [CITATION.cff](CITATION.cff).

## Legacy MATLAB Source

The original MATLAB files are preserved in [legacy/matlab](legacy/matlab) for historical reference and attribution to the 2009 project authors.
