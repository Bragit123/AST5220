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
h = 0.67
H0 = H0_over_h * h

# Things to plot
a = np.exp(cosmo_df["x"])
etaHp_c = cosmo_df["eta"] * cosmo_df["Hp"] / c
OmegaRad = cosmo_df["OmegaR"] + cosmo_df["OmegaNu"]
OmegaMat = cosmo_df["OmegaB"] + cosmo_df["OmegaCDM"]
dHpdx_Hp = cosmo_df["dHpdx"] / cosmo_df["Hp"]
ddHpddx_Hp = cosmo_df["ddHpddx"] / cosmo_df["Hp"]
d_H0 = H0 * cosmo_df["z"]

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
plt.axvline(x=1, linestyle="dashed", color="black", label="Today ($a=1$)")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\eta$ (Mpc)")
plt.legend()
plt.savefig("Figures/BC_eta.pdf")

## Hp
plt.figure()
plt.title("Evolution of conformal Hubble factor")
plt.plot(a, cosmo_df["Hp"])
plt.axvline(x=1, linestyle="dashed", color="black", label="Today ($a=1$)")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\mathcal{H}(x) \\; (\\frac{100 km/s}{Mpc})$")
plt.legend()
plt.savefig("Figures/BC_Hp.pdf")

## dHpdx / Hp
plt.figure()
plt.title("Evolution of Hubble factors (derivatives of $\\mathcal{H}=aH$)")
plt.plot(a, dHpdx_Hp, label="$\\frac{1}{\\mathcal{H}} \\frac{d\\mathcal{H}}{dx}$")
plt.plot(a, ddHpddx_Hp, label="$\\frac{1}{\\mathcal{H}} \\frac{d^2\\mathcal{H}}{dx^2}$")
plt.axhline(y=0, linestyle="dashed", color="black")
plt.axvline(x=1, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\frac{1}{\\mathcal{H}} \\frac{d\\mathcal{H}}{dx}$ , $\\frac{1}{\\mathcal{H}} \\frac{d^2\\mathcal{H}}{dx^2}$")
plt.legend()
plt.savefig("Figures/BC_Hp_derivative.pdf")

## H/H0
plt.figure()
plt.title("Evolution of Hubble factor")
plt.plot(a, cosmo_df["H_over_H0"])
plt.axvline(x=1, linestyle="dashed", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$H/H_0$")
plt.savefig("Figures/BC_H_H0.pdf")

## eta*Hp/c
# Ignore first and last part:
condition = (1e-6 <= a) & (a <= 1)
a_condition = a[condition]
etaHp_c_condition = etaHp_c[condition]
plt.figure()
plt.title("Evolution of the product of conformal time and conformal Hubble factor.")
plt.plot(a_condition, etaHp_c_condition)
plt.axvline(x=1, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\frac{\\eta(x) \\mathcal{H}(x)}{c}$")
plt.savefig("Figures/BC_etaHp_c.pdf")

## Distance measures
plt.figure()
plt.title("Evolution of cosmological distance measures")
plt.plot(cosmo_df["z"], d_H0, color="black", linestyle="dashed", label="Naive Hubble distance ($d = H_0 z$)")
plt.plot(cosmo_df["z"], cosmo_df["dL"], label="Luminosity distance $d_L$")
plt.plot(cosmo_df["z"], cosmo_df["chi"], label="Comoving distance $\\chi$")
plt.plot(cosmo_df["z"], cosmo_df["dA"], label="Angular diameter distance $d_A$")
plt.xlabel("Redshift $z$")
plt.ylabel("Distance (Mpc)")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.savefig("Figures/BC_distance.pdf")

## Omega
plt.figure()
plt.title("Evolution of density parameters")
plt.plot(a, OmegaMat, label="$\\Omega_{Matter}=\\Omega_b+\\Omega_{CDM}$")
plt.plot(a, OmegaRad, label="$\\Omega_{Radiation}=\\Omega_\\gamma+\\Omega_\\nu$")
plt.plot(a, cosmo_df["OmegaLambda"], label="$\\Omega_{DarkEnergy}$")
plt.plot(a, cosmo_df["OmegaK"], label="$\\Omega_{Curvature}$")
plt.axvline(x=1, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\Omega$")
plt.legend()
plt.savefig("Figures/BC_Omega.pdf")