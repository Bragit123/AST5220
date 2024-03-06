# AST5220
This repo contains the code used to compute the CMB spectrum in my project in
AST5220.

## How to run the code
To run the computations you must first compile the C++ code. Do this by running ´make´
in the terminal. After doing this you can run the code using

> ./cmb

The plots are produced by running the different python scripts. The figures are
stored in their respective milestone directory in the Figures/ directory.

### Milestone 1
To produce the plots for milestone 1, run

> python3 plot_BC.py

and

> python3 plot_supernova.py

Note: In Utils.h, there are two options to choose from for x_start and x_end.
These values are the range of x-values that will be considered when doing the
calculations. I have used different x-values when
considering supernova data that when not considering it, in order to fit with
the values given in the supernova dataset.
When not considering the supernova dataset, comment line 55-56, and uncomment
line 52-53.
When considering the supernova dataset, comment line 52-53, and uncomment line 55-56.
Similarly, you must uncomment line 46 in Main.cpp to include the supernova data,
and should uncomment line 42 as well to get the output (and comment line 41 to
avoid outputting twice).

### Milestone 2

### Milestone 3

### Milestone 4