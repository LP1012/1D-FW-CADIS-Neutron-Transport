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

x_vals = df.iloc[2, :].to_numpy()
flux_pl_vals = df.iloc[0, :].to_numpy()
flux_col_vals = df.iloc[1, :].to_numpy()

borders = df.iloc[3, :].to_numpy()

max_pl_flux = max(flux_pl_vals)
max_col_flux = max(flux_col_vals)

# prepare name
name_without_extension = outfile_name.rsplit("_out.", 1)[0]

# create plot
plt.figure(figsize=(8, 6))
sns.lineplot(x=x_vals, y=flux_pl_vals)
plt.vlines(borders, 0, max_pl_flux, linestyles="dashed", colors="black")
plt.xlabel("x [cm]")
plt.ylabel(r"Flux [cm$^{-2}$-s$^{-1}$-source]")
plt.title("MCSlab 1D Path-Length Flux Plot")
plt.savefig(f"plot_{name_without_extension}_pl.png", dpi=300)
plt.close()

plt.figure(figsize=(8, 6))
sns.lineplot(x=x_vals, y=flux_col_vals)
plt.vlines(borders, 0, max_col_flux, linestyles="dashed", colors="black")
plt.xlabel("x [cm]")
plt.ylabel(r"Flux [cm$^{-2}$-s$^{-1}$-source]")
plt.title("MCSlab 1D Collision Estimate Flux Plot")
plt.savefig(f"plot_{name_without_extension}_col.png", dpi=300)
plt.close()
