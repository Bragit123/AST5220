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
#Splines are acting up at beginning and end, which are not too interesting
#anyway. We restrict to get a better view of the interesting part.
ind_start = 10
ind_end = -10
plt.figure()
plt.title("Evolution of optical depth")
plt.plot(a[ind_start:ind_end], recomb_df["tau"][ind_start:ind_end], label="$\\tau$")
plt.plot(a[ind_start:ind_end], -recomb_df["dtaudx"][ind_start:ind_end], label="$-\\frac{d\\tau}{dx}$")
plt.plot(a[ind_start:ind_end], recomb_df["ddtauddx"][ind_start:ind_end], label="$\\frac{d^2\\tau}{dx^2}$")
plt.axvline(x=a_recomb, linestyle="dashed", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Optical depth $\\tau$")
plt.legend()
plt.savefig("Figures/Milestone_2/recomb_tau.pdf")

## g_tilde
# Zero for most a-values, so we restrict the view to easier see the interesting part.
n = len(recomb_df["g_tilde"])
ind_start = int(n/2)
ind_end = int(7*n/11)
plt.figure()
plt.title("Evolution of the visibility function")
plt.plot(a[ind_start:ind_end], recomb_df["g_tilde"][ind_start:ind_end], label="$\\tilde{g}$")
plt.plot(a[ind_start:ind_end], recomb_df["dgdx_tilde"][ind_start:ind_end]/10, label="$\\frac{1}{10}\\frac{d\\tilde{g}}{dx}$")
plt.plot(a[ind_start:ind_end], recomb_df["ddgddx_tilde"][ind_start:ind_end]/100, label="$\\frac{1}{100}\\frac{d^2\\tilde{g}}{dx^2}$")
plt.axvline(x=a_recomb, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Visibility function $\\tilde{g}$")
plt.legend()
plt.savefig("Figures/Milestone_2/recomb_g_tilde.pdf")