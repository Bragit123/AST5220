import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import healpy as hp

cells = pd.read_csv("cells.txt", sep=" ", header=None)
cells.columns = [
    "l",
    "cell",
    "cellN",
    "skrrt"
]

planck_low = pd.read_csv("data/planck_cell_low.txt", sep=r"\s+", comment="#", header=None)
planck_low.columns = [
    "l",
    "Cell",
    "err_up",
    "err_down"
]

planck_high = pd.read_csv("data/planck_cell_high.txt", sep=r"\s+", comment="#", header=None)
planck_high.columns = [
    "l",
    "Cell",
    "err_down",
    "err_up",
    "best_fit"
]

thetas = pd.read_csv("theta.txt", sep=" ", header=None)
thetas.columns = [
    "k",
    "theta6",
    "theta100",
    "theta200",
    "theta500",
    "theta1000",
    "skrrt",
]

matter = pd.read_csv("matter.txt", sep=" ", header=None)
matter.columns = [
    "k",
    "P",
    "skrrt"
]

reid = pd.read_csv("data/reid_DR7.txt", sep=" ", comment="#", header=None)
reid.columns = [
    "k",
    "P",
    "error"
]

wmap = pd.read_csv("data/wmap_act.txt", sep=r"\s+", comment="#", header=None)
wmap.columns = [
    "k",
    "P",
    "P_up"
]

## Recombination and tight coupling scalefactor
x_recomb = -6.98549 # From output of ./cmb
a_recomb = np.exp(x_recomb)

x_tc = -6.99409
a_tc = np.exp(x_tc)

## Density perturbations
plt.figure()
plt.title("Evolution of density perturbations")
plt.plot(cells["l"], cells["cell"], label="Theory prediction")
plt.errorbar(planck_low["l"], planck_low["Cell"], (planck_low["err_down"], planck_low["err_up"]), fmt=".k", label="Planck best fit")
plt.errorbar(planck_high["l"], planck_high["Cell"], (planck_high["err_down"], planck_high["err_up"]), fmt=".k")
plt.xscale("log")
# plt.yscale("log")
plt.xlabel("Multipole $\ell$")
plt.ylabel("$\\frac{\\ell (\\ell + 1)}{2 \\pi}C_{\\ell} (\\mu K)^2$")
plt.legend()
plt.savefig("Figures/Milestone_4/cmb.pdf")

n = int(len(thetas["k"]) / 3)
k = thetas["k"][:n]
ells = [6, 100, 200, 500, 1000]
theta_arrs = [
    thetas["theta6"][:n],
    thetas["theta100"][:n],
    thetas["theta200"][:n],
    thetas["theta500"][:n],
    thetas["theta1000"][:n]
]
theta2_arrs = []
for theta_arr in theta_arrs:
    theta2_arr = theta_arr * theta_arr / k
    theta2_arrs.append(theta2_arr)
## Thetas
plt.figure()
plt.title("Photon multipoles")
for i in range(len(ells)):
    plt.plot(k, theta_arrs[i], label=f"$\\ell = {ells[i]}$")
plt.xlabel("$ck/H_0$")
plt.ylabel("$\\Theta_\\ell$")
plt.legend()
plt.savefig("Figures/Milestone_4/thetas.pdf", bbox_inches="tight")

# n2 = int(n/3)
## Theta^2/k
ymax = [1e-4, 1e-6, 1e-8]
plt.figure()
fig, axs = plt.subplots(3, 1, sharex=True)
for ax_i in range(3):
    ax = axs[ax_i]
    ax.set_ylim(top=ymax[ax_i])
    for i in range(len(ells)):
        # plt.plot(k[:n2], theta2_arrs[i][:n2], label=f"$\\ell = {ells[i]}$")
        if ax_i == 0:
            ax.plot(k, theta2_arrs[i], label=f"$\\ell = {ells[i]}$")
        else:
            # Avoid duplicate legends
            ax.plot(k, theta2_arrs[i])
axs[0].set_title("Square photon multipoles appearing in $C_\\ell$ integral")
axs[2].set_xlabel("$ck/H_0$")
axs[1].set_ylabel("$\\Theta_{\\ell}^2 H_0 / (ck)$")
fig.legend(loc="center right")
fig.savefig("Figures/Milestone_4/theta2s.pdf", bbox_inches="tight")

# # n2 = int(n/3)
# ## Theta^2/k
# plt.figure()
# fig, axs = plt.subplots(3, 1)
# for ax_i in range(3):
#     plt.ylim(0, ymax[ax_i])
#     plt.title("Theta")
#     for i in range(len(ells)):
#         # plt.plot(k[:n2], theta2_arrs[i][:n2], label=f"$\\ell = {ells[i]}$")
#         plt.plot(k, theta2_arrs[i], label=f"$\\ell = {ells[i]}$")
#     plt.xlabel("$ck/H_0$")
#     plt.ylabel("$\\Theta_l^2 H_0 / (ck)$")
#     # plt.xscale("log")
#     # plt.yscale("log")
#     plt.legend()
#     plt.savefig("Figures/Milestone_4/theta2s.pdf", bbox_inches="tight")


Mpc = 3.08567758e22 # From Utils.h
c = 2.99792458e8 # From Utils.h
a_eq = 0.00016 # From plot_BC output
H_over_H0_eq = 391813 # From plot_BC output
H0 = 67e3 / Mpc # From Utils.h
H_eq = H_over_H0_eq * H0
k_eq = a_eq * H_eq / c
k_eq = k_eq * Mpc / 0.67
plt.figure()
plt.title("Matter power spectrum")
ind_start = 10
# plt.plot(matter["k"][ind_start:], matter["P"][ind_start:])
plt.plot(matter["k"][5:], matter["P"][5:], label="Theory prediction")
plt.axvline(x=k_eq, color="black", linestyle="dashed")
# plt.plot(matter["k"], matter["P"] / k, linestyle="dashed", color="blue")
plt.errorbar(reid["k"], reid["P"], yerr=reid["error"], fmt=".", label="Reid")
plt.errorbar(wmap["k"], wmap["P"], yerr=wmap["P_up"]-wmap["P"], fmt=".", label="WMAP")
# plt.errorbar(wmap["k"], wmap_mean, yerr=wmap_err)
# plt.plot(wmap["k"], wmap["P_up"], ".")
plt.xlabel("$k$ $(h/Mpc)$")
plt.ylabel("$P(k, x=0)$ $(Mpc/h)^3$")
plt.xscale("log")
plt.yscale("log")
plt.savefig("Figures/Milestone_4/matter.pdf", bbox_inches="tight")


# ## CMB map
# np.random.seed(100)
# cmb_map = hp.synfast(cells["cell"], 1024)
# plt.figure()
# hp.mollview(cmb_map, min=1e-3, max=1e-1, title="CMB only temperature map", unit="K")
# plt.savefig("Figures/Milestone_4/cmb_map.pdf")