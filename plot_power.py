import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## Load each datafile
n_k = 3

data = pd.read_csv("cells.txt", sep=" ", header=None)
data.columns = [
    "l",
    "cell",
    "cellN",
    "skrrt"
]

## Recombination and tight coupling scalefactor
x_recomb = -6.98549 # From output of ./cmb
a_recomb = np.exp(x_recomb)

x_tc = -6.99409
a_tc = np.exp(x_tc)

## Density perturbations
plt.figure()
plt.title("Evolution of density perturbations")
plt.plot(data["l"], data["cell"], label="Theoretic prediction")
plt.xscale("log")
# plt.yscale("log")
plt.xlabel("Multipole $l$")
plt.ylabel("$C_{\\ell}$")
plt.legend()
plt.savefig("Figures/Milestone_4/cmb.pdf")
