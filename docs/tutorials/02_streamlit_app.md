# Tutorial 2: Interactive Streamlit App

The Streamlit app turns the first simulation into a small local web tool with sliders, labels, hover help, plots, and an export button.

## Install

From the repository root:

```bash
python -m pip install -e ".[app]"
```

For development and tests:

```bash
python -m pip install -e ".[dev,app]"
```

## Run

```bash
streamlit run app.py
```

Or, after installation:

```bash
cnlse-app
```

## Controls

| Control | What it does |
| --- | --- |
| Number of bits | Chooses the length of the simulated bit pattern. |
| Detail per bit | Controls samples per symbol; higher is smoother and slower. |
| Symbol rate | Sets the communication symbol rate in Gbaud. |
| Bit pattern | Chooses De Bruijn or seeded random bits. |
| Laser wavelength | Sets the optical wavelength in nanometers. |
| Nonlinear phase target | Chooses launch power through the nonlinear phase target. |
| Extinction ratio | Controls Mach-Zehnder on/off contrast. |
| Fiber length | Length of the transmit fiber in kilometers. |
| Fiber loss | Attenuation in dB/km. |
| Fiber dispersion | Chromatic dispersion in ps/nm/km. |
| Compensating dispersion | Negative dispersion used to counteract the transmit fiber. |
| Include nonlinear SPM | Toggles self-phase modulation in the scalar solver. |
| Export Plot | Downloads the current graph as a PNG file. |
