import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

pert0_01 = pd.read_csv("perturbations_k0.01.txt", sep=" ", header=None)
pert0_01.columns = [
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

a = np.exp(pert0_01["x"])

x_recomb = -6.98549 # From output of ./cmb
a_recomb = np.exp(x_recomb)

x_tc = -6.99409
a_tc = np.exp(x_tc)

## delta_cdm
plt.figure()
plt.title("Evolution of $\\delta_{\\text{CDM}}$")
plt.plot(a, pert0_01["delta_cdm"], label="$\\delta_{\\text{CDM}}$")
plt.axvline(x=a_tc, linestyle="dashed", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\delta_{\\text{CDM}}$")
plt.legend()
plt.savefig("Figures/Milestone_3/delta_cdm.pdf")

## delta_b
plt.figure()
plt.title("Evolution of $\\delta_{\\text{B}}$")
plt.plot(a, pert0_01["delta_b"], label="$\\delta_{\\text{B}}$")
plt.axvline(x=a_tc, linestyle="dashed", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\delta_{\\text{B}}$")
plt.legend()
plt.savefig("Figures/Milestone_3/delta_b.pdf")

## v_cdm
plt.figure()
plt.title("Evolution of $v_{\\text{CDM}}$")
plt.plot(a, pert0_01["v_cdm"], label="$v_{\\text{CDM}}$")
# plt.axvline(x=a_tc, linestyle="dashed", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$v_{\\text{CDM}}$")
plt.legend()
plt.savefig("Figures/Milestone_3/v_cdm.pdf")

## v_b
plt.figure()
plt.title("Evolution of $v_{\\text{B}}$")
plt.plot(a, pert0_01["v_b"], label="$v_{\\text{B}}$")
plt.axvline(x=a_tc, linestyle="dashed", color="black")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$v_{\\text{B}}$")
plt.legend()
plt.savefig("Figures/Milestone_3/v_b.pdf")

## Theta0
plt.figure()
plt.title("Evolution of $\\Theta_0$")
plt.plot(a, pert0_01["Theta0"], label="$\\Theta_0$")
# plt.axvline(x=a_tc, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\Theta0$")
plt.legend()
plt.savefig("Figures/Milestone_3/Theta0.pdf")

## Theta1
plt.figure()
plt.title("Evolution of $\\Theta_1$")
plt.plot(a, pert0_01["Theta1"], label="$\\Theta_0$")
plt.axvline(x=a_tc, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\Theta1$")
plt.legend()
plt.savefig("Figures/Milestone_3/Theta1.pdf")

## Theta2
plt.figure()
plt.title("Evolution of $\\Theta_2$")
plt.plot(a, pert0_01["Theta2"], label="$\\Theta_0$")
plt.axvline(x=a_tc, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\Theta2$")
plt.legend()
plt.savefig("Figures/Milestone_3/Theta2.pdf")

## Phi
plt.figure()
plt.title("Evolution of $\\Phi$")
plt.plot(a, pert0_01["Phi"], label="$\\Phi$")
plt.axvline(x=a_tc, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\Phi$")
plt.legend()
plt.savefig("Figures/Milestone_3/Phi.pdf")

## Psi
plt.figure()
plt.title("Evolution of $\\Psi$")
plt.plot(a, pert0_01["Psi"], label="$\\Psi$")
plt.axvline(x=a_tc, linestyle="dashed", color="black")
plt.xscale("log")
plt.xlabel("Scalefactor $a$")
plt.ylabel("$\\Psi$")
plt.legend()
plt.savefig("Figures/Milestone_3/Psi.pdf")
