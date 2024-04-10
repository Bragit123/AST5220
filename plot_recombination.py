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

x_recomb = -6.98549 # From output of ./cmb
a_recomb = np.exp(x_recomb)
print(f"a_recomb = {a_recomb}")

## Xe
plt.figure()
plt.title("Evolution of free electron fraction")
plt.plot(a, recomb_df["Xe"], label="$X_e$")
plt.plot(a, recomb_df["Xe_saha"], linestyle="dashed", label="$X_e$ (Saha-approximation)")
plt.axvline(x=a_recomb, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$X_e$")
plt.legend()
plt.savefig("Figures/Milestone_2/recomb_Xe.pdf")

## Xe (y_log)
plt.figure()
plt.title("Evolution of free electron fraction")
plt.plot(a, recomb_df["Xe"], label="$X_e$")
plt.plot(a, recomb_df["Xe_saha"], linestyle="dashed", label="$X_e$ (Saha-approximation)")
plt.axvline(x=a_recomb, linestyle="dashed", color="black")
plt.xscale("log")
plt.yscale("log")
plt.ylim(bottom=1e-10, top=10)
plt.xlabel("Scalefactor $a$")
plt.ylabel("$X_e$")
plt.legend()
plt.savefig("Figures/Milestone_2/recomb_Xe_log.pdf")

## tau
plt.figure()
plt.title("Evolution of optical depth")
plt.plot(a, recomb_df["tau"], label="$\\tau$")
plt.plot(a, -recomb_df["dtaudx"], label="$-\\frac{d\\tau}{dx}$")
plt.plot(a, recomb_df["ddtauddx"], label="$\\frac{d^2\\tau}{dx^2}$")
plt.axvline(x=a_recomb, linestyle="dashed", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Optical depth $\\tau$")
plt.legend()
plt.savefig("Figures/Milestone_2/recomb_tau.pdf")

## g_tilde
plt.figure()
plt.title("Evolution of the visibility function")
plt.plot(a[1500:2500], recomb_df["g_tilde"][1500:2500], label="$\\tilde{g}$")
plt.plot(a[1500:2500], recomb_df["dgdx_tilde"][1500:2500]/10, label="$\\frac{1}{10}\\frac{d\\tilde{g}}{dx}$")
plt.plot(a[1500:2500], recomb_df["ddgddx_tilde"][1500:2500]/100, label="$\\frac{1}{100}\\frac{d^2\\tilde{g}}{dx^2}$")
plt.axvline(x=a_recomb, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Visibility function $\\tilde{g}$")
plt.legend()
plt.savefig("Figures/Milestone_2/recomb_g_tilde.pdf")