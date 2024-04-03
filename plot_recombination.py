import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

recomb_df = pd.read_csv("recombination.txt", sep=" ", header=None)
recomb_df.columns = [
    "x",
    "Xe",
    "ne",
    "Xe_saha",
    "ne_saha",
    "tau",
    "dtaudx",
    "ddtauddx",
    "g_tilde",
    "dgdx_tilde",
    "ddgddx_tilde"
]

a = np.exp(recomb_df["x"])

## Xe
plt.figure()
plt.title("Evolution of free electron fraction")
plt.plot(a, recomb_df["Xe"], label="$X_e$")
plt.plot(a, recomb_df["Xe_saha"], linestyle="dashed", label="$X_e$ (Saha-approximation)")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$X_e$")
plt.legend()
plt.savefig("Figures/Milestone_2/recomb_Xe.pdf")

## tau
plt.figure()
plt.title("Evolution of free electron fraction")
plt.plot(a, recomb_df["tau"], label="$\\tau$")
plt.plot(a, recomb_df["dtaudx"], label="$\\frac{d\\tau}{dx}$")
plt.plot(a, recomb_df["ddtauddx"], label="$\\frac{d^2\\tau}{dx^2}$")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Optical depth $\\tau$")
plt.legend()
plt.savefig("Figures/Milestone_2/recomb_tau.pdf")