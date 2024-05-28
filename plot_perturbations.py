import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## Need conformal Hubble parameter to compare k to horizon
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
x_cosmo = cosmo_df["x"]
a_cosmo = np.exp(x_cosmo)
Hp = cosmo_df["Hp"]
c = 3e8


## Load each datafile
n_k = 3

k0_1 = pd.read_csv("perturbations_k0.1.txt", sep=" ", header=None)
k0_1.columns = [
    "x",
    "delta_cdm",
    "delta_b",
    "v_cdm",
    "v_b",
    "Theta0",
    "Theta1",
    "Theta2",
    "Phi",
    "Psi",
    "S_T",
    "S_T_j5",
    "S_T_j50",
    "S_T_j500"
]

k0_01 = pd.read_csv("perturbations_k0.01.txt", sep=" ", header=None)
k0_01.columns = [
    "x",
    "delta_cdm",
    "delta_b",
    "v_cdm",
    "v_b",
    "Theta0",
    "Theta1",
    "Theta2",
    "Phi",
    "Psi",
    "S_T",
    "S_T_j5",
    "S_T_j50",
    "S_T_j500"
]

k0_001 = pd.read_csv("perturbations_k0.001.txt", sep=" ", header=None)
k0_001.columns = [
    "x",
    "delta_cdm",
    "delta_b",
    "v_cdm",
    "v_b",
    "Theta0",
    "Theta1",
    "Theta2",
    "Phi",
    "Psi",
    "S_T",
    "S_T_j5",
    "S_T_j50",
    "S_T_j500"
]

k_values = [0.001, 0.01, 0.1]
df_list = [k0_001, k0_01, k0_1]
a_list = []
delta_gamma_list = []
v_gamma_list = []

for df in df_list:
    a_list.append(np.exp(df["x"]))
    delta_gamma_list.append(4.0 * df["Theta0"])
    v_gamma_list.append(-3.0 * df["Theta1"])

## Decoupling scalefactor
a_decoup = 0.00092396 # From ./cmb

x_tc = -6.99409
a_tc = np.exp(x_tc)

## Horizon
a_horizon = []
Mpc = 3.08567758e22
for k in k_values:
    horizon_ind = np.argmin(np.abs(Hp - c*k / Mpc))
    a_horizon.append(a_cosmo[horizon_ind])

plot_color_k = ["red", "blue", "green", "black"]

## Density perturbations
plt.figure()
plt.title("Evolution of density perturbations")
for i in range(n_k):
    plt.plot(a_list[i], df_list[i]["delta_cdm"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
    plt.plot(a_list[i], df_list[i]["delta_b"], color=plot_color_k[i], linestyle="dashed")
    plt.plot(a_list[i], delta_gamma_list[i], color=plot_color_k[i], linestyle="dotted")
    plt.axvline(a_horizon[i], linestyle="dashed", color="black")
    plt.text(a_horizon[i], 1000, f"$k={k_values[i]}$", rotation=90, horizontalalignment="right")
plt.axvline(a_decoup, linestyle="dashed", color="black")
plt.text(a_decoup, 1000, "Decoupling", rotation=90, horizontalalignment="right")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\delta_{\\text{CDM}}$ (solid), $\\delta_b$ (dashed), $\\delta_\\gamma$ (dotted)")
plt.legend()
plt.savefig("Figures/Milestone_3/densities.pdf")

## Velocity perturbations
plt.figure()
plt.title("Evolution of velocity perturbations")
plt.ylim(bottom=1e-7, top=1e3)
for i in range(n_k):
    plt.plot(a_list[i], df_list[i]["v_cdm"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
    plt.plot(a_list[i], df_list[i]["v_b"], color=plot_color_k[i], linestyle="dashed")
    plt.plot(a_list[i], v_gamma_list[i], color=plot_color_k[i], linestyle="dotted")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$v_{\\text{CDM}}$ (solid), $v_b$ (dashed), $v_\\gamma$ (dotted)")
plt.legend()
plt.savefig("Figures/Milestone_3/velocities.pdf")

## Photon quadrupole
plt.figure()
plt.title("Evolution of photon quadrupole")
for i in range(n_k):
    plt.plot(a_list[i], df_list[i]["Theta2"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Photon quadrupole $\\Theta_2$")
plt.legend()
plt.savefig("Figures/Milestone_3/quadrupole.pdf")

## Phi
plt.figure()
plt.title("Evolution of spatial gravitational potential")
for i in range(n_k):
    plt.plot(a_list[i], df_list[i]["Phi"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Spatial gravitational potential $\\Phi$")
plt.legend()
plt.savefig("Figures/Milestone_3/potential.pdf")

## Phi + Psi
plt.figure()
plt.title("Evolution of the anisotropic stress")
for i in range(n_k):
    plt.plot(a_list[i], df_list[i]["Phi"] + df_list[i]["Psi"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Anisotropic stress $\\Phi + \\Psi$")
plt.legend()
plt.savefig("Figures/Milestone_3/stress.pdf")

hist_inds = np.argwhere((-8 <= df_list[0]["x"]) * (df_list[0]["x"] <= 0))[:,0] # Indices of "history" (times before today)
## Source function
plt.figure()
plt.title("Source function")
plt.plot(df_list[2]["x"][hist_inds], df_list[2]["S_T_j500"][hist_inds] / 1e-3, color=plot_color_k[2], linestyle="solid", label=f"k={k_values[2]} Mpc")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\tilde{S}(k,x) j_l [ k (\\eta_0 - \\eta(x)) ] / 10^{-3}$")
plt.legend()
plt.savefig("Figures/Milestone_3/source.pdf")
