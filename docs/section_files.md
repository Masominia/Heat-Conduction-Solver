# Cross-section files

One statement per line. `#` starts a comment; the comment lines at the top of a file are its description and appear in the report. Lengths are in micrometres; the window is centred on x = 0 and its bottom edge is the heat sink.

## Statements

| statement | meaning | default |
|---|---|---|
| `title <text>` | name used in reports | empty |
| `grid dx=<um> dy=<um>` | cell size | 0.25, 0.25 |
| `window width=<um> ymin=<um> ymax=<um>` | must be a whole number of cells | 60, -20, 5 |
| `material <name> k=<W/(m K)> c=<J/(m3 K)>` | define or redefine a material | |
| `background <material>` | what fills the window before painting | sio2 |
| `layer <material> <ymin> <ymax>` | a layer across the whole window | |
| `rect <material> <xmin> <xmax> <ymin> <ymax>` | a rectangle; later statements paint over earlier ones | |
| `heater <name> xmin= xmax= ymin= ymax= power=<mW/um>` | a region that dissipates heat uniformly | power 0.1 |
| `probe <name> xmin= xmax= ymin= ymax=` | mean temperature over a rectangle | |
| `probe <name> x= y=` | temperature of one cell | |
| `top adiabatic` or `top convective h=<W/(m2 K)>` | top surface | adiabatic |
| `solver direct` or `solver cg` | banded Cholesky or conjugate gradients | direct |
| `optics wavelength=<um> dndt=<1/K>` | used for the power for pi | 1.55, 1.86e-4 |
| `transient heater=<name> dt=<us> steps=<n>` | step response of one heater | none |
| `expect <kind> <args>` | a check, see below | |

Built-in materials: `si` (148, 1.63e6), `sio2` (1.38, 1.61e6), `sin` (20, 2.2e6), `tin` (28, 3.2e6), `al` (237, 2.42e6), `air` (0.026, 1.2e3).

A heater statement only adds heat; paint the heater metal with a `rect` of the same extent if its conductivity matters. A cell belongs to a rectangle when its centre does, so choose cell sizes that divide the dimensions of the structure. Heater and probe names share one namespace.

## Expectations

| kind | arguments | passes when |
|---|---|---|
| `temperature` | probe, K, relative tolerance | with every heater at its stated power |
| `response` | probe, heater, K per mW/um, relative tolerance | |
| `pi_power` | probe, heater, mW, relative tolerance | |
| `crosstalk` | probe, reference probe, heater, ratio, relative tolerance | ratio of the two responses |
| `hotter` | probe, probe | the first is warmer than the second |
| `peak` | K, relative tolerance | highest temperature anywhere |
| `balance` | tolerance | heat out equals heat in |
| `surface_share` | fraction, tolerance | share of the heat leaving through the top |
| `tau` | probe, us, relative tolerance | time to 63% of the step response |
| `rise` | probe, us, relative tolerance | 10 to 90% rise time |
| `cells` | count | the mesh has that many cells |

A malformed expectation, or one naming a probe or heater that does not exist, is reported as a failed check with its line number; the other checks still run.

## Output

`heatphase check <file>` prints `key value` lines, one `probe` line per probe, one `response` line per probe and heater (with the power for pi), one `step` line per probe if there is a transient, and one `check PASS|FAIL <name> :: <detail>` line per expectation. Exit status: 0 if all checks pass, 1 if one fails, 2 if the file cannot be read or is invalid. `heatphase run <file>` adds the conductivity and temperature of every cell row by row from the bottom, and the step-response traces.

## Example

    # Two phase shifters 10 um apart with a trench between them.
    title Trench to the substrate
    grid dx=0.25 dy=0.25
    window width=80 ymin=-10 ymax=5
    layer si -10 0
    rect si -0.25 0.25 2 2.25
    rect si 9.75 10.25 2 2.25
    rect tin -1 1 4 4.25
    rect air 4 6 0 5
    heater h1 xmin=-1 xmax=1 ymin=4 ymax=4.25 power=0.1
    probe wg1 xmin=-0.25 xmax=0.25 ymin=2 ymax=2.25
    probe wg2 xmin=9.75 xmax=10.25 ymin=2 ymax=2.25
    expect hotter wg1 wg2
    expect balance 1e-9
