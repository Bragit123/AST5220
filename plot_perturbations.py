import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
    "siu"
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
    "siu"
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
    "siu"
]

# k0_0001 = pd.read_csv("perturbations_k0.0001.txt", sep=" ", header=None)
# k0_0001.columns = [
#     "x",
#     "delta_cdm",
#     "delta_b",
#     "v_cdm",
#     "v_b",
#     "Theta0",
#     "Theta1",
#     "Theta2",
#     "Phi",
#     "Psi",
#     "siu"
# ]

# k_values = [0.0001, 0.001, 0.01, 0.1]
# df_list = [k0_0001, k0_001, k0_01, k0_1]
k_values = [0.001, 0.01, 0.1]
df_list = [k0_001, k0_01, k0_1]
a_list = []
delta_gamma_list = []
v_gamma_list = []

# ## scalefactor from x
# a0_1 = np.exp(k0_1["x"])
# a0_01 = np.exp(k0_01["x"])
# a0_001 = np.exp(k0_001["x"])
# a0_0001 = np.exp(k0_0001["x"])

# ## delta_gamma
# delta_gamma0_1 = 4.0 * k0_1["Theta0"]
# delta_gamma0_01 = 4.0 * k0_01["Theta0"]
# delta_gamma0_001 = 4.0 * k0_001["Theta0"]
# delta_gamma0_0001 = 4.0 * k0_0001["Theta0"]

# ## v_gamma
# v_gamma0_1 = -3.0 * k0_1["Theta1"]
# v_gamma0_01 = -3.0 * k0_01["Theta1"]
# v_gamma0_001 = -3.0 * k0_001["Theta1"]
# v_gamma0_0001 = -3.0 * k0_0001["Theta1"]

for df in df_list:
    a_list.append(np.exp(df["x"]))
    delta_gamma_list.append(4.0 * df["Theta0"])
    v_gamma_list.append(-3.0 * df["Theta1"])

## Recombination and tight coupling scalefactor
x_recomb = -6.98549 # From output of ./cmb
a_recomb = np.exp(x_recomb)

x_tc = -6.99409
a_tc = np.exp(x_tc)

plot_color_k = ["red", "blue", "green", "black"]

## Density perturbations
plt.figure()
plt.title("Evolution of density perturbations")
for i in range(n_k):
    plt.plot(a_list[i], df_list[i]["delta_cdm"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
    plt.plot(a_list[i], df_list[i]["delta_b"], color=plot_color_k[i], linestyle="dashed")
    plt.plot(a_list[i], delta_gamma_list[i], color=plot_color_k[i], linestyle="dotted")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\delta_{\\text{CDM}}$ (solid), $\\delta_b$ (dashed), $\\delta_\\gamma$ (dotted)")
plt.legend()
plt.savefig("Figures/Milestone_3/densities.pdf")

## Velocity perturbations
plt.figure()
plt.title("Evolution of velocity perturbations")
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
    # plt.plot(a_list[i], delta_gamma_list[i], color=plot_color_k[i])
    # plt.plot(a_list[i], v_gamma_list[i], color=plot_color_k[i], linestyle="dashed")
    plt.plot(a_list[i], df_list[i]["Theta2"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
    # plt.axvline(x=a_recomb, color="black", linestyle="dashed")
plt.xscale("log")
# plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Photon quadrupole $\\Theta_2$")
plt.legend()
plt.savefig("Figures/Milestone_3/quadrupole.pdf")

## Phi
plt.figure()
plt.title("Evolution of spatial gravitational potential")
for i in range(n_k):
    plt.plot(a_list[i], df_list[i]["Phi"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
    # plt.axvline(x=a_recomb, color="black", linestyle="dashed")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Spatial gravitational potential $\\Phi$")
plt.legend()
plt.savefig("Figures/Milestone_3/potential.pdf")

## Phi + Psi
plt.figure()
plt.title("Evolution of the anisotropic stress")
for i in range(n_k):
    # plt.plot(a_list[i], df_list[i]["Phi"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
    plt.plot(a_list[i], df_list[i]["Phi"] + df_list[i]["Psi"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
    # plt.axvline(x=a_recomb, color="black", linestyle="dashed")
plt.xscale("log")
# plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("Anisotropic stress $\\Phi + \\Psi$")
plt.legend()
plt.savefig("Figures/Milestone_3/stress.pdf")

# ## Theta0
# plt.figure()
# plt.title("$\\Theta_0$")
# for i in range(n_k):
#     # plt.plot(a_list[i], df_list[i]["Phi"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
#     plt.plot(a_list[i], df_list[i]["Theta0"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
# plt.xscale("log")
# # plt.yscale("log")
# plt.xlabel("Scalefactor $a$")
# plt.ylabel("$\\Theta_0$")
# plt.legend()
# plt.savefig("Figures/Milestone_3/Theta0.pdf")

# ## Theta1
# plt.figure()
# plt.title("$\\Theta_1$")
# for i in range(n_k):
#     # plt.plot(a_list[i], df_list[i]["Phi"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
#     plt.plot(a_list[i], df_list[i]["Theta1"], color=plot_color_k[i], linestyle="solid", label=f"k={k_values[i]} Mpc")
# plt.xscale("log")
# # plt.yscale("log")
# plt.xlabel("Scalefactor $a$")
# plt.ylabel("$\\Theta_1$")
# plt.legend()
# plt.savefig("Figures/Milestone_3/Theta1.pdf")