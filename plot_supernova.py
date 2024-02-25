import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import seaborn as sns

supernova_data = pd.read_csv("data/supernovadata.txt", sep="\s+", comment="#", header=None)
supernova_data.columns = ["z", "dL", "error"]

cosmo_df = pd.read_csv("cosmology_supernova.txt", sep=" ", header=None)
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

def gaussian(x, mu, sigma):
    return 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-(x-mu)**2 / (2*sigma**2))

chi2, h, OmegaM, OmegaK = np.loadtxt("results_supernovafitting.txt", skiprows=200, comments="#", unpack=True) # Skip burnin of chains
OmegaLambda = 1 - OmegaM - OmegaK # Assuming OmegaR = OmegaNu = 0

# Find best values
chi2_min_arg = np.argmin(chi2)
chi2_min = chi2[chi2_min_arg]
OmegaM_min = OmegaM[chi2_min_arg]
OmegaK_min = OmegaK[chi2_min_arg]
OmegaLambda_min = OmegaLambda[chi2_min_arg]
h_min = h[chi2_min_arg]

# Good values (1 sigma confidence region)
good_condition = chi2 < chi2_min + 3.53
chi2_good = chi2[good_condition]
OmegaM_good = OmegaM[good_condition]
OmegaK_good = OmegaK[good_condition]
OmegaLambda_good = OmegaLambda[good_condition]
h_good = h[good_condition]

# Standard deviations
chi2_std = np.std(chi2_good)
OmegaM_std = np.std(OmegaM_good)
OmegaK_std = np.std(OmegaK_good)
OmegaLambda_std = np.std(OmegaLambda_good)
h_std = np.std(h_good)

# Mean values
chi2_mean = np.mean(chi2_good)
OmegaM_mean = np.mean(OmegaM_good)
OmegaK_mean = np.mean(OmegaK_good)
OmegaLambda_mean = np.mean(OmegaLambda_good)
h_mean = np.mean(h_good)

# Find values to consider for the Gaussian
n_points = 100
OmegaM_range = np.linspace(np.min(OmegaM_good), np.max(OmegaM_good), n_points)
OmegaK_range = np.linspace(np.min(OmegaK_good), np.max(OmegaK_good), n_points)
OmegaLambda_range = np.linspace(np.min(OmegaLambda_good), np.max(OmegaLambda_good), n_points)
h_range = np.linspace(np.min(h_good), np.max(h_good), n_points)

n_bins = 30 # Number of bins to use for histograms

## OmegaM
plt.figure()
plt.title("Posterior for $\\Omega_M$")
plt.hist(x=OmegaM_good, bins=n_bins, density=True)
plt.plot(OmegaM_range, gaussian(OmegaM_range, OmegaM_mean, OmegaM_std))
plt.axvline(x=OmegaM_min, color="black", linestyle="dashed", label="Planck best-fit value")
plt.xlabel("$\\Omega_M$")
plt.ylabel("Probability density")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_OmegaM.pdf")

## OmegaK
plt.figure()
plt.title("Posterior for $\\Omega_K$")
plt.hist(x=OmegaK_good, bins=n_bins, density=True)
plt.plot(OmegaK_range, gaussian(OmegaK_range, OmegaK_mean, OmegaK_std))
plt.axvline(x=OmegaK_min, color="black", linestyle="dashed", label="Planck best-fit value")
plt.xlabel("$\\Omega_K$")
plt.ylabel("Probability density")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_OmegaK.pdf")

## OmegaLambda
plt.figure()
plt.title("Posterior for $\\Omega_\\Lambda$")
plt.hist(x=OmegaLambda_good, bins=n_bins, density=True)
plt.plot(OmegaLambda_range, gaussian(OmegaLambda_range, OmegaLambda_mean, OmegaLambda_std))
plt.axvline(x=OmegaLambda_min, color="black", linestyle="dashed", label="Planck best-fit value")
plt.xlabel("$\\Omega_\\Lambda$")
plt.ylabel("Probability density")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_OmegaLambda.pdf")

## h
plt.figure()
plt.title("Posterior for $h$")
plt.hist(x=h_good, bins=n_bins, density=True)
plt.plot(h_range, gaussian(h_range, h_mean, h_std))
plt.axvline(x=h_min, color="black", linestyle="dashed", label="Planck best-fit value")
plt.xlabel("$h$")
plt.ylabel("Probability density")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_h.pdf")

## Scatterplot
plt.figure()
sns.scatterplot(x=OmegaM_good, y=OmegaK_good)
plt.xlabel("$\\Omega_M$")
plt.ylabel("$\\Omega_K$")
plt.savefig("Figures/Milestone_1/supernova_scatter.pdf")

# Units (copied from utils.h)
Mpc = 3.08567758e22
H0_over_h = 100 * 1e3 / Mpc
c = 3e8
Gpc = 1000 * Mpc

dL_Gpc = cosmo_df["dL"]/Gpc

## Luminosity distance
plt.figure()
plt.title("Luminosity distance compared to supernova data.")
plt.errorbar(supernova_data["z"], supernova_data["dL"]/supernova_data["z"], yerr=supernova_data["error"]/supernova_data["z"], capsize=3, barsabove=True, fmt=".", label="Supernova data")
plt.plot(cosmo_df["z"], dL_Gpc/cosmo_df["z"], label="Theoretical prediction")
plt.xscale("log")
plt.ylabel("Distance $d_L/z$ (Gpc)")
plt.xlabel("Redshift z")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_dL.pdf")