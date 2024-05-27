import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
# import seaborn as sns

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

# Units (copied from utils.h)
Mpc = 3.08567758e22
H0_over_h = 100 * 1e3 / Mpc
c = 3e8
Gpc = 1000 * Mpc

# Find best values
chi2_min_arg = np.argmin(chi2)
chi2_min = chi2[chi2_min_arg]
OmegaM_min = OmegaM[chi2_min_arg]
OmegaK_min = OmegaK[chi2_min_arg]
OmegaLambda_min = OmegaLambda[chi2_min_arg]
h_min = h[chi2_min_arg]

# 1 sigma confidence values
sigma1_condition = chi2 < chi2_min + 3.53
chi2_sigma1 = chi2[sigma1_condition]
OmegaM_sigma1 = OmegaM[sigma1_condition]
OmegaK_sigma1 = OmegaK[sigma1_condition]
OmegaLambda_sigma1 = OmegaLambda[sigma1_condition]
h_sigma1 = h[sigma1_condition]

# 2 sigma confidence values
sigma2_condition = chi2 < chi2_min + 8.02
chi2_sigma2 = chi2[sigma2_condition]
OmegaM_sigma2 = OmegaM[sigma2_condition]
OmegaK_sigma2 = OmegaK[sigma2_condition]
OmegaLambda_sigma2 = OmegaLambda[sigma2_condition]
h_sigma2 = h[sigma2_condition]

# Standard deviations
chi2_std = np.std(chi2_sigma1)
OmegaM_std = np.std(OmegaM_sigma1)
OmegaK_std = np.std(OmegaK_sigma1)
OmegaLambda_std = np.std(OmegaLambda_sigma1)
h_std = np.std(h_sigma1)

# Mean values
chi2_mean = np.mean(chi2_sigma1)
OmegaM_mean = np.mean(OmegaM_sigma1)
OmegaK_mean = np.mean(OmegaK_sigma1)
OmegaLambda_mean = np.mean(OmegaLambda_sigma1)
h_mean = np.mean(h_sigma1)

# Find values to consider for the Gaussian
n_points = 100
OmegaM_range = np.linspace(np.min(OmegaM_sigma1), np.max(OmegaM_sigma1), n_points)
OmegaK_range = np.linspace(np.min(OmegaK_sigma1), np.max(OmegaK_sigma1), n_points)
OmegaLambda_range = np.linspace(np.min(OmegaLambda_sigma1), np.max(OmegaLambda_sigma1), n_points)
h_range = np.linspace(np.min(h_sigma1), np.max(h_sigma1), n_points)

n_bins = 30 # Number of bins to use for histograms

## OmegaM
plt.figure()
plt.title("Posterior for $\\Omega_M$")
plt.hist(x=OmegaM_sigma1, bins=n_bins, density=True)
plt.plot(OmegaM_range, gaussian(OmegaM_range, OmegaM_mean, OmegaM_std))
plt.axvline(x=OmegaM_min, color="black", linestyle="dashed", label="Planck best-fit value")
plt.xlabel("$\\Omega_M$")
plt.ylabel("Probability density")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_OmegaM.pdf")

## OmegaK
plt.figure()
plt.title("Posterior for $\\Omega_K$")
plt.hist(x=OmegaK_sigma1, bins=n_bins, density=True)
plt.plot(OmegaK_range, gaussian(OmegaK_range, OmegaK_mean, OmegaK_std))
plt.axvline(x=OmegaK_min, color="black", linestyle="dashed", label="Planck best-fit value")
plt.xlabel("$\\Omega_K$")
plt.ylabel("Probability density")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_OmegaK.pdf")

## OmegaLambda
plt.figure()
plt.title("Posterior for $\\Omega_\\Lambda$")
plt.hist(x=OmegaLambda_sigma1, bins=n_bins, density=True)
plt.plot(OmegaLambda_range, gaussian(OmegaLambda_range, OmegaLambda_mean, OmegaLambda_std))
plt.axvline(x=OmegaLambda_min, color="black", linestyle="dashed", label="Planck best-fit value")
plt.xlabel("$\\Omega_\\Lambda$")
plt.ylabel("Probability density")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_OmegaLambda.pdf")

## H0
plt.figure()
plt.title("Posterior PDF for $H_0$")
plt.hist(x=h_sigma1, bins=n_bins, density=True)
plt.plot(h_range, gaussian(h_range, h_mean, h_std))
plt.axvline(x=h_min, color="black", linestyle="dashed", label="Best-fit value")
plt.xlabel("$H_0$ (100km/Mpc)")
plt.ylabel("Probability density")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_H0.pdf")

## Scatterplot ML
plt.figure()
plt.title("Scatterplot of $\\Omega_M$ and $\\Omega_\\Lambda$ from MCMC samples")
plt.scatter(x=OmegaM_sigma2, y=OmegaLambda_sigma2, label="$2\\sigma$ constraint", rasterized=True)
plt.scatter(x=OmegaM_sigma1, y=OmegaLambda_sigma1, label="$1\\sigma$ constraint", rasterized=True)
plt.plot([0.0,1.0], [1.0,0.0], color="black", linestyle="dashed", label="Flat universe")
plt.xlabel("$\\Omega_M$")
plt.ylabel("$\\Omega_\\Lambda$")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_scatter_MLambda.pdf")

## Luminosity distance
dL_Gpc = cosmo_df["dL"]/Gpc
plt.figure()
plt.title("Luminosity distance compared to supernova data.")
plt.errorbar(supernova_data["z"], supernova_data["dL"]/supernova_data["z"], yerr=supernova_data["error"]/supernova_data["z"], capsize=3, barsabove=True, fmt=".", label="Supernova data")
plt.plot(cosmo_df["z"], dL_Gpc/cosmo_df["z"], label="Theoretical prediction")
plt.xscale("log")
plt.ylabel("Distance $d_L/z$ (Gpc)")
plt.xlabel("Redshift z")
plt.legend()
plt.savefig("Figures/Milestone_1/supernova_dL.pdf")