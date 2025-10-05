import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import seaborn as sns
import sys
import pandas as pd

sns.set_theme(
    context="paper",
    style="ticks",
    palette="colorblind",
    font="serif",
    font_scale=1.2,
    color_codes=True,
    rc={"lines.linewidth": 2.5},
)
sns.despine()

# import csv file
if len(sys.argv) != 2:
    raise ValueError("Error! Must pass input file name.")

outfile_name = sys.argv[1]
df = pd.read_csv(outfile_name, header=None)

x_vals = df.iloc[1, :].to_numpy()
flux_vals = df.iloc[0, :].to_numpy()
borders = df.iloc[2, :].to_numpy()

max_flux = max(flux_vals)

# create plot
plt.figure(figsize=(8, 6))
sns.lineplot(x=x_vals, y=flux_vals)
plt.vlines(borders, 0, max_flux, linestyles="dashed", colors="black")
plt.xlabel("x [cm]")
plt.ylabel(r"Flux [cm$^{-2}$-s$^{-1}$-source]")
plt.title("MCSlab 1D Flux Plot")
plt.savefig(f"plot_{outfile_name}.png", dpi=300)
