import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

cosmo_df = pd.read_csv("cosmology.txt", sep=" ", header=None)
cosmo_df.columns = [
    "x",
    "z",
    "eta",
    "t",
    "H_over_H0",
    "Hp",
    "dHpdx",
    "ddHpddx",
    "dL",
    "chi",
    "dA",
    "OmegaB",
    "OmegaCDM",
    "OmegaLambda",
    "OmegaR",
    "OmegaNu",
    "OmegaK"
]

# Units (copied from utils.h)
Mpc = 3.08567758e22
H0_over_h = 100 * 1e3 / Mpc
c = 3e8

# Things to plot
a = np.exp(cosmo_df["x"])
cosmo_df["etaHp_c"] = cosmo_df["eta"] * cosmo_df["Hp"] / c
cosmo_df["OmegaRel"] = cosmo_df["OmegaR"] + cosmo_df["OmegaNu"]
cosmo_df["OmegaMat"] = cosmo_df["OmegaB"] + cosmo_df["OmegaCDM"]
cosmo_df["dHpdx_Hp"] = cosmo_df["dHpdx"] / cosmo_df["Hp"]
cosmo_df["ddHpddx_Hp"] = cosmo_df["ddHpddx"] / cosmo_df["Hp"]

# Scale data to proper units
cosmo_df["eta"] = cosmo_df["eta"] / Mpc
cosmo_df["Hp"] = cosmo_df["Hp"] / H0_over_h
cosmo_df["dL"] = cosmo_df["dL"] / Mpc
cosmo_df["chi"] = cosmo_df["chi"] / Mpc
cosmo_df["dA"] = cosmo_df["dA"] / Mpc

## eta
plt.figure()
plt.title("Evolution of conformal time")
plt.plot(a, cosmo_df["eta"])
plt.axvline(x=1, linestyle="dashed", color="black", label="Today")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\eta$ (Mpc)")
plt.savefig("Figures/BC_eta.pdf")

## Hp
plt.figure()
plt.title("$\\mathcal{H}(x) \\; (\\frac{100 km/s}{Mpc})$")
sns.lineplot(data=cosmo_df, x="x", y="Hp")
plt.axvline(x=0, linestyle="dashed", color="black")
plt.yscale("log")
plt.savefig("Figures/BC_Hp.pdf")

## dHpdx / Hp
plt.figure()
plt.title("$\\frac{1}{\\mathcal{H}} \\frac{d\\mathcal{H}}{dx} \\; , \\; \\frac{1}{\\mathcal{H}} \\frac{d^2\\mathcal{H}}{dx^2}$")
sns.lineplot(data=cosmo_df, x="x", y="dHpdx_Hp")
sns.lineplot(data=cosmo_df, x="x", y="ddHpddx_Hp")
plt.axhline(y=0, linestyle="dashed", color="black")
plt.axvline(x=0, linestyle="dashed", color="black")
plt.savefig("Figures/BC_Hp_derivative.pdf")

## H/H0
plt.figure()
plt.title("$H/H_0$")
sns.lineplot(data=cosmo_df, x="x", y="H_over_H0")
plt.axvline(x=0, linestyle="dashed", color="black")
plt.savefig("Figures/BC_H_H0.pdf")

## eta*Hp/c
plt.figure()
plt.title("$\\frac{\\eta(x) \\mathcal{H}(x)}{c}$")
sns.lineplot(data=cosmo_df, x="x", y="etaHp_c")
plt.axvline(x=0, linestyle="dashed", color="black")
plt.yscale("log")
plt.savefig("Figures/BC_etaHp_c.pdf")

## Distance measures
plt.figure()
plt.title("$d_L \; , \; \\chi \; , \; d_A$")
sns.lineplot(data=cosmo_df, x="z", y="dL", label="$d_L$")
sns.lineplot(data=cosmo_df, x="z", y="chi", label="$\\chi$")
sns.lineplot(data=cosmo_df, x="z", y="dA", label="$d_A$")
plt.ylabel("Distance")
plt.xscale("log")
plt.yscale("log")
plt.savefig("Figures/BC_distance.pdf")

## Omega
plt.figure()
plt.title("$\\Omega_i$")
sns.lineplot(data=cosmo_df, x="x", y="OmegaRel", label="$\\Omega_{Relativistic}=\\Omega_\\gamma+\\Omega_\\nu$")
sns.lineplot(data=cosmo_df, x="x", y="OmegaMat", label="$\\Omega_{Matter}=\\Omega_b+\\Omega_{CDM}$")
sns.lineplot(data=cosmo_df, x="x", y="OmegaLambda", label="$\\Omega_\\Lambda$")
sns.lineplot(data=cosmo_df, x="x", y="OmegaK", label="$\\Omega_K$", linestyle="dashed")
plt.axvline(x=0, linestyle="dashed", color="black")
plt.legend()
plt.savefig("Figures/BC_Omega.pdf")