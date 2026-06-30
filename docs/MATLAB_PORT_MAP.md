# MATLAB Port Map

This document tracks how the original MATLAB thesis simulator maps into the Python `cnlse` package.

Original MATLAB project:

`09_Fiber_Optics_Software_Matlab_Files/`

Original authors:

- Ryan Maher Hamze
- Krayem Ahmad Ajjawi
- Meriam Hawwa

Current Python package implementation:

- Ryan Maher Hamze

## Porting Policy

- Preserve credit to the original MATLAB authors in documentation.
- Prefer readable Python and NumPy idioms over line-for-line translation.
- Add tests for every small utility before moving to larger solvers.
- Keep examples educational and reproducible.

## File Mapping

| MATLAB file | Python target | Status |
| --- | --- | --- |
| `Reset_All.m` | `cnlse.state`, `cnlse.constants` | Started |
| `Pattern.m` | `cnlse.patterns` | Started |
| `Nmod.m` | `cnlse.mathutils.nmod` | Ported |
| `Fastshift.m` | `cnlse.mathutils.fastshift` | Ported |
| `Fastexp.m` | `cnlse.mathutils.fastexp` | Ported |
| `Phi2pow.m` | `cnlse.link.phi_to_power` | Ported |
| `Electricsource.m` | `cnlse.electric` | Started |
| `Lasersource.m` | `cnlse.laser` | Started |
| `Mz_Modulator.m` | `cnlse.modulators` | Started |
| `Create_Field.m` | `cnlse.fields` | Started |
| `Avg_Power.m` | `cnlse.power` | Started |
| `Ampliflat.m` | `cnlse.amplifiers` | Planned |
| `Fiber.m` | `cnlse.fiber` | Started: scalar attenuation/GVD/SPM |
| `Plotfield.m` | `cnlse.plotting` | Planned |
| `Checkfields.m` | Python dataclass validation | Planned |
| `Examples/Ex01.m` | `examples/ook_100km_compensated.py` | Started |

## Solver Port Order

1. Small math helpers.
2. Pattern and waveform generation.
3. Laser source and modulator.
4. Field container/state model.
5. Average power and plotting.
6. Scalar NLSE fiber propagation.
7. CNLSE/polarization propagation.
8. GUI/tutorial controls.
