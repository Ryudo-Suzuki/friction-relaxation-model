# Relaxation of Sliding Friction: Data and Source Code

## Overview

This repository contains the numerical datasets and source code used in the paper:

**“Relaxation of sliding friction from a statistical model of aging contacts”**

**Data and code prepared by:**  
Ryudo Suzuki  
Email: suzuki.ryudo.55e@st.kyoto-u.ac.jp

---

## Repository Structure

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

---

## Description of Contents

### figures/

Each directory (`fig2`--`fig6`) contains:

- **data/**  
  Numerical data used to generate the corresponding figure.

- **fig/**  
  Base PDF plots generated from the numerical data.

- **load.plt**  
  Gnuplot script used to reproduce the base plots.

The final figures in the manuscript may include additional layout adjustments
to these base plots.

### src/

Contains the C source code used to generate the numerical data for each figure.

- **fig2/**  
  Finite-N simulations and analytical calculation of friction relaxation
  after a velocity step.

- **fig3/**  
  Calculation of contact-time distributions at selected times.

- **fig4/**  
  Calculation of relaxation functions for different power-law exponents and
  upper cutoffs.

- **fig5/**  
  Calculation of the steady-state friction as a function of sliding velocity.

- **fig6/**  
  Calculation of the dependence of relaxation and characteristic slip distance
  on the power-law exponent.

---

## Reproducing the Plots

The numerical data required for the figures are included in the repository.

For example, Fig. 2 can be reproduced by running

```sh
cd figures/fig2
gnuplot load.plt
```

The resulting PDF is written to the corresponding `fig/` directory.
The same procedure can be used for Figs. 3--6.

---

## Regenerating the Numerical Data

The C source files in `src/` generate the numerical datasets contained in the
corresponding `figures/fig*/data/` directories.

Each source file includes information on how to compile and run the program.

A C11 compiler such as GCC or Clang is required.
Gnuplot with the `pdfcairo` terminal is used to generate the plots.

The finite-N simulation for Fig. 2 may require relatively long computation
time. A fixed random seed is used for reproducibility.

