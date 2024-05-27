# AST5220
This repo contains the code used to compute the CMB spectrum in my project in
AST5220.

## How to run the code
To run the computations you must first compile the C++ code. Do this by running `make`
in the terminal. After doing this you can run the code using `./cmb`
The plots are produced by running the different python scripts. The figures are
stored in their respective milestone directory in the Figures/ directory.

To produce the plots for milestone 1, run `python3 plot_BC.py` and
`python3 plot_supernova.py`. To produce the plots from milestone 2, run
`python3 plot_recombination.py`. To produce the plots from milestone 3, run
`python3 plot_perturbations.py`. To produce the plots from milestone 4, run `python3 plot_power.py`.

## Note on the supernova part
Note: In BackgroundCosmology.h, there are two options to choose from for n_x, x_start and x_end.
These values are the range of x-values that will be considered when doing the
calculations. I have used different x-values when
considering supernova data that when not considering it, in order to fit with
the values given in the supernova dataset.
When not considering the supernova dataset, comment line 32-34, and uncomment
line 28-30.
When considering the supernova dataset, comment line 28-30, and uncomment line 32-34.
Similarly, you must uncomment line 46 in Main.cpp to include the supernova data,
and you should uncomment line 42 as well to get the output (and comment line 41 to
avoid outputting twice).
When running the supernova code, the other milestones will not work properly, as
they depend on the x-values I have chosen when not considering the supernova.
Therefore, run the code for the supernova first, and quit the program after it
is finished. Then run the code for all milestones without considering the
supernova to get everything to work.

### How to install GSL
In order to run the code presented in this repository, you will need to install
GSL. The following recipe on how to do this is directly copied from the [github
repo](https://github.com/HAWinther/AST5220-Cosmology) by Hans A. Winther, where the template code this project is built on was provided.

See [this](https://solarianprogrammer.com/) for how to install it on a Windows machine. On Linux or a Mac you can either use a package manager or install it directly as follows:

- Go the the home directory:

cd $HOME

- Make a local folder to hold libraries:

mkdir local

- Enter this directory:

cd local

- Download the code (if you don't have wget you need to get the file to this dir by other means):

wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz

- Untar the code:

tar -xvf gsl-2.6.tar.gz

- You should now have the gsl-2.6 folder. Enter it:

cd gsl-2.6

- Run the configure script:

./configure --prefix=$HOME/local

- Compile and install it:

make ; make install

- In the CMB code Makefile change the include and lib paths to point to the library:

INC  = -I$(HOME)/local/include
LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas

- If this fails with "libgsl.so not found" then run the command:

export LD\_LIBRARY\_PATH="$LD\_LIBRARY\_PATH:$HOME/local/lib"

and try to run ./cmb again and it should work. To avoid having
to run this command every time you open a new terminal open
the $HOME/.bashrc file and add this line to the end of the file
and it will load everytime you open a new window.