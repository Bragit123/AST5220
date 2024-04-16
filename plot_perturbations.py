import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

pert0_01 = pd.read_csv("perturbations_k0.01.txt", sep=" ", header=None)
pert0_01.columns = [
    "x",
    "Phi",
    "siu"
]

a = np.exp(pert0_01["x"])

x_recomb = -6.98549 # From output of ./cmb
a_recomb = np.exp(x_recomb)

x_tc = -6.99409
a_tc = np.exp(x_tc)

## Xe
plt.figure()
plt.title("Evolution of $\\Phi$")
plt.plot(a, pert0_01["Phi"], label="$\\Phi$")
plt.axvline(x=a_tc, linestyle="dashed", color="black")
plt.xscale("log")
# plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\Phi$")
plt.legend()
plt.savefig("Figures/Milestone_3/Phi.pdf")
