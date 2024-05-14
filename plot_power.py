import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

cells = pd.read_csv("cells.txt", sep=" ", header=None)
cells.columns = [
    "l",
    "cell",
    "cellN",
    "skrrt"
]

thetas = pd.read_csv("theta.txt", sep=" ", header=None)
thetas.columns = [
    "k",
    "theta6",
    "theta100",
    "theta200",
    "theta500",
    "theta1000",
    "skrrt",
]


## Recombination and tight coupling scalefactor
x_recomb = -6.98549 # From output of ./cmb
a_recomb = np.exp(x_recomb)

x_tc = -6.99409
a_tc = np.exp(x_tc)

## Density perturbations
plt.figure()
plt.xlim((1, 3000))
plt.ylim((-1000, 8000))
plt.title("Evolution of density perturbations")
plt.plot(cells["l"], cells["cell"], label="Theoretic prediction")
plt.xscale("log")
# plt.yscale("log")
plt.xlabel("Multipole $\ell$")
plt.ylabel("$C_{\\ell}$")
plt.legend()
plt.savefig("Figures/Milestone_4/cmb.pdf")

n = int(len(thetas["k"]) / 3)
## Thetas
plt.figure()
plt.title("Theta")
plt.plot(thetas["k"][:n], thetas["theta6"][:n], label="Theoretic prediction")
plt.plot(thetas["k"][:n], thetas["theta100"][:n], label="Theoretic prediction")
plt.plot(thetas["k"][:n], thetas["theta200"][:n], label="Theoretic prediction")
plt.plot(thetas["k"][:n], thetas["theta500"][:n], label="Theoretic prediction")
plt.plot(thetas["k"][:n], thetas["theta1000"][:n], label="Theoretic prediction")
# plt.xscale("log")
# plt.yscale("log")
plt.xlabel("$ck/H_0$")
plt.ylabel("$\\Theta_l^2 H_0 / (ck)$")
plt.legend()
plt.savefig("Figures/Milestone_4/thetas.pdf")