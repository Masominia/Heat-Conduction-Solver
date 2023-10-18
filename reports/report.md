# HeatPhase report

80 cross-sections, 1266976 cells, 1294 checks: 1294 PASS, 0 FAIL.

Every cross-section is a text file in `cases/`. Its expectations come from an independent solution of the same mesh or from a closed form, stated in the comments at the top of the file.
Per-section results are in `reports/cases/`, temperature maps in `reports/figures/cases/`.

## Families

| family | sections | checks | passed |
|---|---|---|---|
| arrays | 6 | 651 | 651 |
| line heater | 5 | 60 | 60 |
| mesh and solver | 12 | 131 | 131 |
| metal | 4 | 24 | 24 |
| crosstalk | 10 | 70 | 70 |
| phase shifter | 18 | 94 | 94 |
| layered stacks | 6 | 42 | 42 |
| top surface | 3 | 15 | 15 |
| step response | 8 | 59 | 59 |
| trenches and undercut | 8 | 60 | 60 |

## Agreement with the independent solution

Each cross-section is read, meshed, assembled and solved again by `reference.py` with SciPy's sparse direct solver. The field is compared with one solve with all heaters on, relative to its peak; the printed field has seven digits.

| quantity | largest relative difference over the library |
|---|---|
| response | 4.82e-12 |
| temperature | 1.13e-12 |
| time | 1.14e-08 |
| peak | 5.12e-12 |
| balance | 1.00e-12 |
| outflow | 1.06e-12 |
| field | 2.95e-07 |

## Studies

### Crosstalk against pitch, 2 um buried oxide

Fitted decay length between 6 and 20 um: 3.54 um (2 d / pi for 5 um of oxide over a perfect sink: 3.18 um).

| pitch, um | own waveguide, K/(mW/um) | neighbour, K/(mW/um) | crosstalk |
|---|---|---|---|
| 4 | 171.5346 | 79.3490 | 0.4626 |
| 6 | 171.7566 | 44.1460 | 0.2570 |
| 8 | 171.8199 | 24.2980 | 0.1414 |
| 10 | 171.8381 | 13.4208 | 0.0781 |
| 15 | 171.8453 | 3.1929 | 0.0186 |
| 20 | 171.8456 | 0.8477 | 4.93e-03 |
| 30 | 171.8457 | 0.0956 | 5.56e-04 |

### Crosstalk against pitch, 3 um buried oxide

Fitted decay length between 6 and 20 um: 4.03 um (2 d / pi for 6 um of oxide over a perfect sink: 3.82 um).

| pitch, um | own waveguide, K/(mW/um) | neighbour, K/(mW/um) | crosstalk |
|---|---|---|---|
| 6 | 218.9626 | 70.2240 | 0.3207 |
| 10 | 219.0905 | 25.5973 | 0.1168 |
| 20 | 219.1090 | 2.1714 | 9.91e-03 |

### Efficiency against speed

| buried oxide, um | power for pi, mW | time constant, us | product, nJ |
|---|---|---|---|
| 1 | 37.936 | 6.339 | 240.471 |
| 2 | 24.246 | 9.310 | 225.737 |
| 3 | 19.016 | 12.593 | 239.476 |

### Trenches and undercut (neighbour at 10 um)

| variant | power for pi, mW | crosstalk at 10 um | time constant, us |
|---|---|---|---|
| none | 24.2476 | 0.0781 | 9.6047 |
| clad | 22.4500 | 0.0212 | n/a |
| box | 21.5627 | 0.0108 | n/a |
| deep | 21.4613 | 7.43e-03 | n/a |
| both | 19.0546 | 0.0105 | 12.6040 |
| narrow | 13.6046 | 8.19e-03 | n/a |
| undercut | 8.4180 | 0.0906 | 41.2962 |
| undercut_trench | 0.8073 | 6.77e-03 | n/a |

### Mesh (relative to 125 nm cells)

| variant | response, K/(mW/um) | change |
|---|---|---|
| 500 nm cells | 179.7810 | 0.0496 |
| 250 nm cells | 171.8771 | 3.41e-03 |

### Window (relative to 60 um wide, 10 um of substrate)

| variant | response, K/(mW/um) | change |
|---|---|---|
| 30 um wide | 172.0322 | 1.08e-03 |
| 40 um wide | 171.8771 | 1.76e-04 |
| 100 um wide | 171.8456 | -7.57e-06 |
| 2 um of substrate | 169.5707 | -0.0132 |
| 5 um of substrate | 170.6571 | -6.92e-03 |
| 20 um of substrate | 173.3141 | 8.54e-03 |

### Line heater on a slab: solver over closed form

| section | distance over thickness | solver over closed form |
|---|---|---|
| line_d2 | 0.5312 | 0.9991 |
| line_d2 | 1.0312 | 0.9996 |
| line_d2 | 2.0312 | 1.0006 |
| line_d2 | 3.0312 | 1.0018 |
| line_d3_k02 | 0.5208 | 0.9996 |
| line_d3_k02 | 1.0208 | 0.9998 |
| line_d3_k02 | 2.0208 | 1.0003 |
| line_d3_k02 | 3.0208 | 1.0008 |
| line_d4 | 0.5156 | 0.9998 |
| line_d4 | 1.0156 | 0.9999 |
| line_d4 | 2.0156 | 1.0001 |
| line_d4 | 3.0156 | 1.0004 |
| line_d4_k5 | 0.5156 | 0.9998 |
| line_d4_k5 | 1.0156 | 0.9999 |
| line_d4_k5 | 2.0156 | 1.0001 |
| line_d4_k5 | 3.0156 | 1.0004 |
| line_d6 | 0.5104 | 0.9999 |
| line_d6 | 1.0104 | 0.9999 |
| line_d6 | 2.0104 | 1.0001 |
| line_d6 | 3.0104 | 1.0002 |

### Heat leaving through the top surface

| film coefficient, W/(m2 K) | share of heat leaving through the top | power for pi, mW |
|---|---|---|
| 10.00000 | 3.01e-05 | 24.24686 |
| 1,000.00000 | 3.00e-03 | 24.29412 |
| 100,000.00000 | 0.22075 | 28.60870 |

## All cross-sections

| section | cells | heaters | peak, K | power for pi, mW | checks | passed |
|---|---|---|---|---|---|---|
| array_all_on | 19,200 | 8 | 65.313 | n/a | 76 | 76 |
| array_alternate | 19,200 | 8 | 45.518 | n/a | 75 | 75 |
| array_one_on | 19,200 | 8 | 40.751 | n/a | 75 | 75 |
| array_ramp | 19,200 | 8 | 73.687 | n/a | 75 | 75 |
| array_ramp_p10 | 28,800 | 8 | 60.632 | n/a | 75 | 75 |
| array_sixteen | 24,960 | 16 | 48.083 | n/a | 275 | 275 |
| line_d2 | 4,096 | 1 | 77.168 | 19.444 | 12 | 12 |
| line_d3_k02 | 9,216 | 1 | 597.006 | 2.764 | 12 | 12 |
| line_d4 | 16,384 | 1 | 93.159 | 18.888 | 12 | 12 |
| line_d4_k5 | 16,384 | 1 | 25.712 | 68.434 | 12 | 12 |
| line_d6 | 36,864 | 1 | 102.512 | 18.707 | 12 | 12 |
| mesh_d125 | 38,400 | 1 | 41.084 | 24.325 | 5 | 5 |
| mesh_d250 | 9,600 | 1 | 41.451 | 24.242 | 5 | 5 |
| mesh_d500 | 2,400 | 1 | 40.581 | 23.176 | 5 | 5 |
| metal_doped | 14,400 | 2 | 19.981 | n/a | 6 | 6 |
| metal_line | 19,200 | 1 | 40.035 | 25.412 | 6 | 6 |
| metal_plane | 21,760 | 1 | 21.889 | 49.471 | 6 | 6 |
| metal_via | 19,200 | 1 | 39.044 | 26.657 | 6 | 6 |
| pitch_box2_p04 | 19,200 | 1 | 41.404 | 24.291 | 7 | 7 |
| pitch_box2_p06 | 19,200 | 1 | 41.434 | 24.259 | 7 | 7 |
| pitch_box2_p08 | 19,200 | 1 | 41.444 | 24.250 | 7 | 7 |
| pitch_box2_p10 | 19,200 | 1 | 41.447 | 24.248 | 7 | 7 |
| pitch_box2_p15 | 19,200 | 1 | 41.448 | 24.247 | 7 | 7 |
| pitch_box2_p20 | 19,200 | 1 | 41.448 | 24.247 | 7 | 7 |
| pitch_box2_p30 | 19,200 | 1 | 41.448 | 24.247 | 7 | 7 |
| pitch_box3_p06 | 20,480 | 1 | 45.659 | 19.029 | 7 | 7 |
| pitch_box3_p10 | 20,480 | 1 | 45.676 | 19.018 | 7 | 7 |
| pitch_box3_p20 | 20,480 | 1 | 45.678 | 19.016 | 7 | 7 |
| shifter_aluminium | 14,400 | 1 | 40.786 | 24.299 | 5 | 5 |
| shifter_box1 | 13,440 | 1 | 36.231 | 37.936 | 5 | 5 |
| shifter_box3 | 15,360 | 1 | 45.678 | 19.016 | 5 | 5 |
| shifter_box4 | 16,320 | 1 | 49.242 | 16.173 | 5 | 5 |
| shifter_clad045 | 13,920 | 1 | 45.604 | 22.459 | 5 | 5 |
| shifter_clad080 | 17,280 | 1 | 34.251 | 29.204 | 5 | 5 |
| shifter_clad120 | 21,120 | 1 | 32.516 | 31.043 | 5 | 5 |
| shifter_gap075 | 14,400 | 1 | 32.558 | 20.397 | 5 | 5 |
| shifter_gap125 | 14,400 | 1 | 36.589 | 22.627 | 5 | 5 |
| shifter_gap250 | 14,400 | 1 | 52.661 | 25.439 | 5 | 5 |
| shifter_no_core | 14,400 | 1 | 41.651 | 24.131 | 5 | 5 |
| shifter_offset | 14,400 | 1 | 41.470 | 25.905 | 5 | 5 |
| shifter_slab | 57,600 | 1 | 36.076 | 43.925 | 6 | 6 |
| shifter_standard | 14,400 | 1 | 41.448 | 24.246 | 8 | 8 |
| shifter_width01 | 14,400 | 1 | 48.990 | 23.520 | 5 | 5 |
| shifter_width04 | 14,400 | 1 | 32.885 | 26.832 | 5 | 5 |
| shifter_width06 | 14,400 | 1 | 27.580 | 30.244 | 5 | 5 |
| shifter_width10 | 14,400 | 1 | 21.006 | 38.040 | 5 | 5 |
| solver_cg_array | 19,200 | 8 | 45.518 | n/a | 75 | 75 |
| solver_cg_standard | 14,400 | 1 | 41.448 | 24.246 | 5 | 5 |
| solver_cg_trench | 19,200 | 1 | 48.782 | 19.055 | 6 | 6 |
| stack_custom | 1,728 | 1 | 45.242 | 17.071 | 7 | 7 |
| stack_metal | 1,728 | 1 | 8.839 | 71.597 | 7 | 7 |
| stack_nitride | 1,728 | 1 | 10.435 | 56.186 | 7 | 7 |
| stack_oxide2 | 2,304 | 1 | 5.943 | 142.526 | 7 | 7 |
| stack_oxide5 | 2,880 | 1 | 15.001 | 34.776 | 7 | 7 |
| stack_oxide_only | 576 | 1 | 8.681 | 73.600 | 7 | 7 |
| top_forced_air | 14,400 | 1 | 41.379 | 24.294 | 5 | 5 |
| top_liquid | 14,400 | 1 | 36.030 | 28.609 | 5 | 5 |
| top_still_air | 14,400 | 1 | 41.447 | 24.247 | 5 | 5 |
| transient_box1 | 13,440 | 1 | 36.231 | 37.936 | 7 | 7 |
| transient_box3 | 15,360 | 1 | 45.678 | 19.016 | 7 | 7 |
| transient_close | 14,400 | 1 | 32.558 | 20.397 | 7 | 7 |
| transient_neighbour | 19,200 | 1 | 41.448 | 24.247 | 10 | 10 |
| transient_slab_d1 | 256 | 1 | 17.833 | 23.365 | 7 | 7 |
| transient_slab_d2 | 256 | 1 | 35.666 | 11.683 | 7 | 7 |
| transient_slab_d3 | 256 | 1 | 53.499 | 7.788 | 7 | 7 |
| transient_standard | 14,400 | 1 | 41.448 | 24.246 | 7 | 7 |
| trench_both | 19,200 | 1 | 48.782 | 19.055 | 10 | 10 |
| trench_box | 19,200 | 1 | 44.838 | 21.563 | 6 | 6 |
| trench_clad | 19,200 | 1 | 43.762 | 22.450 | 6 | 6 |
| trench_deep | 19,200 | 1 | 44.929 | 21.461 | 6 | 6 |
| trench_narrow | 19,200 | 1 | 62.932 | 13.605 | 6 | 6 |
| trench_none | 19,200 | 1 | 41.447 | 24.248 | 10 | 10 |
| trench_undercut | 19,200 | 1 | 69.590 | 8.418 | 10 | 10 |
| trench_undercut_trench | 19,200 | 1 | 538.232 | 0.807 | 6 | 6 |
| window_sub02 | 6,720 | 1 | 41.232 | 24.572 | 5 | 5 |
| window_sub05 | 9,600 | 1 | 41.333 | 24.415 | 5 | 5 |
| window_sub20 | 24,000 | 1 | 41.593 | 24.041 | 5 | 5 |
| window_w030 | 7,200 | 1 | 41.470 | 24.220 | 5 | 5 |
| window_w040 | 9,600 | 1 | 41.451 | 24.242 | 5 | 5 |
| window_w100 | 24,000 | 1 | 41.448 | 24.247 | 5 | 5 |
