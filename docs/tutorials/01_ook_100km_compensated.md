# Tutorial 1: OOK Over 100 km With Dispersion Compensation

This tutorial recreates the first public demo from the original fiber optics senior project in a Python-first style.

## What This Simulation Shows

The model sends a simple on-off-keying signal through an optical fiber link:

1. Generate a binary bit pattern.
2. Turn that bit pattern into an electrical OOK drive signal.
3. Create a laser optical field.
4. Modulate the laser with a Mach-Zehnder modulator.
5. Propagate the optical field through 100 km of fiber.
6. Propagate a second copy through 100 km of fiber plus a dispersion-compensating fiber.
7. Plot the input, uncompensated output, and compensated output.

## Install

From the repository root:

```bash
python -m pip install -e ".[dev]"
```

## Run

```bash
python examples/ook_100km_compensated.py
```

The output plot is written to:

```text
outputs/ook_100km_compensated.png
```

## Main Concepts

### Bit Pattern

The bit pattern is the digital information being sent. In this demo, a compact De Bruijn pattern is used because it gives useful bit transitions in a short sequence.

### Electrical OOK Signal

OOK means "on-off keying." A `1` turns the optical power on, and a `0` turns it off.

### Laser Source

The laser creates a constant optical field at the selected wavelength, usually around 1550 nm for fiber communication examples.

### Mach-Zehnder Modulator

The modulator turns the laser field into a data-carrying optical signal using the electrical drive signal.

### Fiber

The fiber applies:

- attenuation, which reduces signal power;
- chromatic dispersion, which spreads pulses in time;
- nonlinear phase rotation, which models self-phase modulation.

### Compensating Fiber

The compensating fiber has negative dispersion. It is selected to counteract the positive dispersion accumulated in the transmit fiber.

## Future GUI Control Labels

These labels are intended for the future small web/desktop interface:

| Control | Plain-language label | Tooltip |
| --- | --- | --- |
| `symbols` | Number of bits | How many digital symbols to simulate. |
| `samples_per_symbol` | Detail per bit | More samples makes smoother plots but runs slower. |
| `symbol_rate_gbaud` | Symbol rate | Transmission speed in billions of symbols per second. |
| `wavelength_nm` | Laser wavelength | Optical carrier wavelength, usually near 1550 nm. |
| `alpha_db_per_km` | Fiber loss | Power loss per kilometer. |
| `dispersion_ps_nm_km` | Fiber dispersion | How strongly the fiber spreads pulses by wavelength. |
| `dphi_max_rad` | Nonlinear step limit | Smaller values are more accurate but slower. |
| `dz_max_m` | Maximum step length | Caps each split-step propagation segment. |
| `Run Simulation` | Run Simulation | Generate signals, propagate through fiber, and update plots. |
| `Export Plot` | Export Plot | Save the current plot to PNG or SVG. |

## Code Sketch

```python
from cnlse.examples import run_ook_100km_compensated

run_ook_100km_compensated("outputs/ook_100km_compensated.png")
```
