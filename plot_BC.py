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
Gyr = 1e9*365*24*60*60

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
cosmo_df["t"] = cosmo_df["t"] / Gyr

### Time-values at different epochs
## Radiation-Matter equality
arg_RM_eq = np.argmin(np.abs(OmegaRad[a<1] - OmegaMat[a<1])) # Radiation-Matter equality
a_RM_eq = a[arg_RM_eq]
z_RM_eq = cosmo_df["z"][arg_RM_eq]
t_RM_eq = cosmo_df["t"][arg_RM_eq]
eta_RM_eq = cosmo_df["eta"][arg_RM_eq]
H_over_H0_RM_eq = cosmo_df["H_over_H0"][arg_RM_eq]

## Matter-Dark energy equality
ind_start = 500 # Matter and dark energy have equal densities 0 in the beginning, so must ignore first indices.
ML_diff = np.abs(OmegaMat[ind_start:]-cosmo_df["OmegaLambda"][ind_start:])
arg_ML_eq = np.argmin(ML_diff) + ind_start # Matter-Dark energy equality
a_ML_eq = a[arg_ML_eq]
z_ML_eq = cosmo_df["z"][arg_ML_eq]
t_ML_eq = cosmo_df["t"][arg_ML_eq]
eta_ML_eq = cosmo_df["eta"][arg_ML_eq]


## Universe starts accelerating
arg_acc = np.argmin(np.abs(cosmo_df["dHpdx"])) # Matter-Dark energy equality
a_acc = a[arg_acc]
z_acc = cosmo_df["z"][arg_acc]
t_acc = cosmo_df["t"][arg_acc]
eta_acc = cosmo_df["eta"][arg_acc]

## Today
arg_today = np.argmin(np.abs(a-1))
a_today = a[arg_today]
z_today = cosmo_df["z"][arg_today]
t_today = cosmo_df["t"][arg_today]
eta_today = cosmo_df["eta"][arg_today]


## eta
plt.figure()
plt.title("Evolution of conformal time")
plt.plot(a, cosmo_df["eta"]/c)
plt.axvline(x=1, linestyle="dashed", color="black")
plt.axvline(x=a_acc, linestyle="dotted", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\eta$ (Mpc)")
plt.savefig("Figures/Milestone_1/BC_eta.pdf")

## Hp
plt.figure()
plt.title("Evolution of conformal Hubble factor")
plt.plot(a, cosmo_df["Hp"])
plt.axvline(x=1, linestyle="dashed", color="black", label="Today ($a=1$)")
plt.axvline(x=a_acc, linestyle="dotted", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\mathcal{H}(x) \\; (\\frac{100 km/s}{Mpc})$")
plt.legend()
plt.savefig("Figures/Milestone_1/BC_Hp.pdf")

## dHpdx / Hp
plt.figure()
plt.title("Evolution of Hubble factors (derivatives of $\\mathcal{H}=aH$)")
plt.plot(a, dHpdx_Hp, label="$\\frac{1}{\\mathcal{H}} \\frac{d\\mathcal{H}}{dx}$")
plt.plot(a, ddHpddx_Hp, label="$\\frac{1}{\\mathcal{H}} \\frac{d^2\\mathcal{H}}{dx^2}$")
plt.axhline(y=0, linestyle="dashed", color="black")
plt.axvline(x=1, linestyle="dashed", color="black")
# plt.axvline(x=a_acc, linestyle="dotted", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\frac{1}{\\mathcal{H}} \\frac{d\\mathcal{H}}{dx}$ , $\\frac{1}{\\mathcal{H}} \\frac{d^2\\mathcal{H}}{dx^2}$")
plt.legend()
plt.savefig("Figures/Milestone_1/BC_Hp_derivative.pdf")

## H/H0
plt.figure()
plt.title("Evolution of Hubble factor")
plt.plot(a, cosmo_df["H_over_H0"])
plt.axvline(x=1, linestyle="dashed", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$H/H_0$")
plt.savefig("Figures/Milestone_1/BC_H_H0.pdf")

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
plt.savefig("Figures/Milestone_1/BC_etaHp_c.pdf")

## Distance measures
plt.figure()
plt.title("Evolution of cosmological distance measures")
plt.plot(cosmo_df["z"], cosmo_df["dL"], label="Luminosity distance $d_L$")
plt.plot(cosmo_df["z"], cosmo_df["chi"], label="Comoving distance $\\chi$")
plt.plot(cosmo_df["z"], cosmo_df["dA"], label="Angular diameter distance $d_A$")
plt.plot(cosmo_df["z"], d_H0, color="black", linestyle="dashed", label="Naive Hubble distance ($d = H_0 z$)")
# plt.axline(xy1=(0,0), slope=H0)
plt.xlabel("Redshift $z$")
plt.ylabel("Distance (Mpc)")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.savefig("Figures/Milestone_1/BC_distance.pdf")

## Cosmic time
plt.figure()
plt.title("Evolution of cosmic time")
plt.plot(a, cosmo_df["t"])
plt.axvline(x=1, linestyle="dashed", color="black", label="Today ($a=1$)")
plt.axvline(x=a_acc, linestyle="dotted", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Cosmic time $t$ (Gyr)")
plt.legend()
plt.savefig("Figures/Milestone_1/BC_t.pdf")

## Omega
plt.figure()
plt.title("Evolution of density parameters")
plt.plot(a, OmegaMat, label="$\\Omega_M=\\Omega_b+\\Omega_{CDM}$")
plt.plot(a, OmegaRad, label="$\\Omega_R=\\Omega_\\gamma+\\Omega_\\nu$")
plt.plot(a, cosmo_df["OmegaLambda"], label="$\\Omega_\\Lambda$")
# plt.plot(a, cosmo_df["OmegaK"], label="$\\Omega_K$")
plt.axvline(x=1, linestyle="dashed", color="black", label="Today ($a=1$)")
# plt.axvline(x=a_RM_eq, linestyle="dashed", color="gray", label="$\\Omega_M = \\Omega_R$")
# plt.axvline(x=a_ML_eq, linestyle="dotted", color="gray", label="$\\Omega_M = \\Omega_\\Lambda$")
# plt.axvline(x=a_acc, linestyle="dotted", color="black", label="Acceleration starts")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\Omega$")
plt.legend()
plt.savefig("Figures/Milestone_1/BC_Omega.pdf")



print("Radiation-matter equality")
print(f"  a = {a_RM_eq:.5f}")
print(f"  z = {z_RM_eq:4.5f}")
print(f"  t = {t_RM_eq:.5f}")
print(f"  eta = {eta_RM_eq:.5f}")
print(f"  H/H0 = {H_over_H0_RM_eq:.5f}")

print("Matter-dark energy equality")
print(f"  a = {a_ML_eq:.5f}")
print(f"  z = {z_ML_eq:4.5f}")
print(f"  t = {t_ML_eq:.5f}")
print(f"  eta = {eta_ML_eq:.5f}")

print("Universe begins to accelerate")
print(f"  a = {a_acc:.5f}")
print(f"  z = {z_acc:4.5f}")
print(f"  t = {t_acc:.5f}")
print(f"  eta = {eta_acc:.5f}")

print("Today")
print(f"  a = {a_today:.5f}")
print(f"  z = {z_today:4.5f}")
print(f"  t = {t_today:.5f}")
print(f"  eta = {eta_today:.5f}")

print("Dark energy non-zero")
print(f"  a = {a_DE:.5f}")