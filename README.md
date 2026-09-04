# Relaxation of Sliding Friction: Data and Source Code

## Overview

This repository contains the numerical datasets, C source code, Gnuplot
scripts, and base PDF figures used for Figures 2-6 of the manuscript:

**"Relaxation of sliding friction from a statistical model of aging contacts"**

**Data and code prepared by:**  
Ryudo Suzuki  
Email: suzuki.ryudo.55e@st.kyoto-u.ac.jp

## Repository structure

```text
.
|-- figures/
|   |-- fig2/
|   |   |-- data/
|   |   |-- fig/
|   |   `-- load.plt
|   |-- fig3/
|   |   |-- data/
|   |   |-- fig/
|   |   `-- load.plt
|   |-- fig4/
|   |   |-- data/
|   |   |-- fig/
|   |   `-- load.plt
|   |-- fig5/
|   |   |-- data/
|   |   |-- fig/
|   |   `-- load.plt
|   `-- fig6/
|       |-- data/
|       |-- fig/
|       `-- load.plt
`-- src/
    |-- fig2/
    |-- fig3/
    |-- fig4/
    |-- fig5/
    `-- fig6/
```

Each `data/` directory contains plain-text numerical data. Comment lines at
the beginning of each file record the parameters and column definitions. Each
`fig/` directory contains the corresponding PDF base plot. The final figures
in the manuscript may include additional layout adjustments.

## Contents

| Figure | Numerical content | Source code |
| --- | --- | --- |
| Fig. 2 | Finite-N and analytical friction relaxation after a velocity step | `finite_n_relaxation.c`, `analytic_relaxation.c` |
| Fig. 3 | Contact-time distributions at selected times | `contact_time_distribution.c` |
| Fig. 4 | Relaxation functions for several exponents and upper cutoffs | `relaxation_curve.c` |
| Fig. 5 | Steady-state friction as a function of velocity | `steady_state_friction.c` |
| Fig. 6 | Alpha dependence of relaxation and characteristic slip | `relaxation_alpha_scan.c`, `characteristic_slip.c` |

## Requirements

- A C11 compiler, such as GCC or Clang
- Gnuplot with the `pdfcairo` terminal

The plotting scripts request Times New Roman. Gnuplot may substitute another
font if it is unavailable.

## Reproducing the plots

The numerical data are included. Run each plotting script from its figure
directory, for example:

```sh
cd figures/fig2
gnuplot load.plt
```

The resulting PDF is written to the local `fig/` directory. Repeat the same
command in `fig3`, `fig4`, `fig5`, and `fig6` to reproduce the other plots.

## Regenerating the numerical data

Each C file begins with its build command, run command, output location, and a
short description of the calculation. Build programs from the repository root
and run them with the figure directory as the working directory.

For Figure 4, `relaxation_curve.c` generates one parameter combination per
run. Set `PARAM_ALPHA` and `PARAM_LMAX` to the combinations referenced by
`figures/fig4/load.plt` to regenerate the complete dataset.

The finite-N calculation for Figure 2 may require substantial computation
time. It uses a fixed random seed for reproducibility.
