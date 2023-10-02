# Methodology

HeatPhase solves heat conduction in the cross-section of a photonic chip: a heater above a waveguide, the oxide around them, the silicon substrate below. It answers what a thermo-optic phase shifter costs in power, how fast it is, and how much it disturbs its neighbours.

## 1. Model

The cross-section is the plane (x, y) perpendicular to the waveguides; everything is uniform along them. The steady temperature rise T above the heat sink obeys

    div(k grad T) + q = 0

with k the local conductivity and q the heat generated in the heaters. The bottom of the window is the heat sink (T = 0), the sides are adiabatic, the top is adiabatic or loses heat through a film coefficient h. Because the problem is linear, the temperature is the sum of the fields of the individual heaters, each proportional to its power.

Units: micrometres, kelvin, heater power in mW per micrometre of heater length (0.1 mW/um is 20 mW in a 200 um heater), microseconds.

## 2. Discretisation

Finite volumes on a uniform mesh of rectangular cells, one material per cell (a cell takes the material of the last rectangle that contains its centre). The conductance between neighbouring cells uses the harmonic mean of their conductivities, which is exact for heat flowing across a material interface that lies on a cell face. The bottom row is tied to the sink through half a cell; with a film coefficient the top row is tied to the surroundings through half a cell in series with the film.

A heater is a rectangle that dissipates its power uniformly; each cell receives the share of the area it covers, so the total is exact whatever the mesh. A probe is the area-weighted mean temperature over a rectangle (for a waveguide core, its cross-section) or the value of one cell.

For a layered stack heated over its whole width the scheme reproduces the one-dimensional result exactly, to rounding.

## 3. Solvers

The matrix is symmetric positive definite. The default solver stores its lower band (cells are numbered along the shorter side of the mesh), factorises it once by Cholesky and solves for each heater in turn. A preconditioned conjugate-gradient solver works on the same matrix without factorising it; `solver cg` selects it, and the library uses it on three cross-sections as an internal cross-check.

## 4. Step response

Switching a heater on at t = 0 from a cold chip:

    C dT/dt = -A T + q

with C the heat capacity of each cell. Implicit Euler: (C / dt + A) T_new = C / dt T_old + q, factorised once. The solver reports, for every probe, the time to reach 1 - 1/e of its steady value and the 10 to 90% rise time. Implicit Euler is first order: with a step of a twentieth of the time constant the reported time is a few percent long.

## 5. What is reported

- response R: temperature rise of a probe in K per (mW/um) dissipated in a heater;
- power for a phase shift of pi: P_pi = lambda / (2 dn/dT R), independent of the heater length, with dn/dT = 1.86e-4 per K for silicon at 1.55 um by default;
- crosstalk: the response of a neighbouring waveguide over that of the heated one;
- the heat balance: power in, power to the sink, power through the top;
- time constants of the step response.

## 6. Validation

1. Rust unit tests: Cholesky against known solutions, conjugate gradients against Cholesky, conductances and boundary terms of the matrix by hand, exactness on a layered stack, heat balance, the slowest decay of a slab against 4 d^2 / (pi^2 alpha) within 1%.
2. Closed forms in the library: layered stacks (exact to 1e-9); a line heater on a slab over a sink, T = P / (pi k) ln coth(pi x / 4d), reproduced within 0.2% from half a thickness to three thicknesses away; the step response of a slab from its series, within 3%.
3. Independent solution: every cross-section is read, meshed, assembled and solved again by `reference.py` with SciPy's sparse direct solver, which shares no code with the Rust solver. Responses, temperatures, heat flows, time constants and the whole field are compared. The reference field is one solve with all heaters on, so the comparison also tests superposition.
4. Consistency across cross-sections: crosstalk falling with pitch at the rate set by the oxide thickness, the product of power and time constant staying put when the buried oxide changes, insensitivity to window and substrate depth, mesh refinement.

What this does not validate is the model itself: see the README.
