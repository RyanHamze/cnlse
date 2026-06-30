"""Fiber optics modulation simulations using NLSE/CNLSE models."""

from .constants import C, E_CHARGE, H_PLANCK, K_BOLTZMANN
from .electric import ook_signal, pulse_shape
from .fields import create_field
from .fiber import FiberParameters, FiberResult, propagate_scalar
from .laser import channel_wavelengths, laser_source
from .link import phi_to_power
from .mathutils import fastexp, fastshift, nmod
from .modulators import mach_zehnder
from .patterns import de_bruijn_binary, random_bits
from .power import average_power
from .plotting import plot_power_traces
from .state import SimulationState, create_state

__all__ = [
    "C",
    "E_CHARGE",
    "FiberParameters",
    "FiberResult",
    "H_PLANCK",
    "K_BOLTZMANN",
    "SimulationState",
    "average_power",
    "channel_wavelengths",
    "create_field",
    "create_state",
    "de_bruijn_binary",
    "fastexp",
    "fastshift",
    "laser_source",
    "mach_zehnder",
    "nmod",
    "ook_signal",
    "phi_to_power",
    "pulse_shape",
    "plot_power_traces",
    "propagate_scalar",
    "random_bits",
]

__version__ = "0.1.0"
