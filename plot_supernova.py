import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

supernova_data = pd.read_csv("data/supernovadata.txt", sep="\s+", comment="#", header=None)
supernova_data.columns = ["z", "dL", "Error"]

# supernova_fit = pd.read_csv("results_supernovafitting.txt", sep="\s+", comment="#", skiprows=range(1,200), header=None)
# supernova_fit.columns = ["chi2", "h", "OmegaM", "OmegaK"]
# chi2_min_ind = np.argmin(supernova_fit["chi2"])
# chi2_min = np.min(supernova_fit["chi2"][chi2_min_ind])
# good_indices = 

# OmegaM = supernova_fit["OmegaM"]
# OmegaK = supernova_fit["OmegaK"]
# print(supernova_fit)

chi2, h, OmegaM, OmegaK = np.loadtxt("results_supernovafitting.txt", skiprows=200, comments="#", unpack=True) # Skip burnin of chains
chi2_min_arg = np.argmin(chi2)
chi2_min = chi2[chi2_min_arg]


good_condition = chi2 < chi2_min + 3.53
chi2_good = chi2[good_condition]
OmegaM_good = OmegaM[good_condition]
OmegaK_good = OmegaK[good_condition]

chi2_std = np.std(chi2_good)

print(chi2_min/31)
print(chi2_std)

plt.figure()
sns.scatterplot(x=OmegaM_good, y=OmegaK_good)
plt.savefig("Figures/supernova_scatter.pdf")