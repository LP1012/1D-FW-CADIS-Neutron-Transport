import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sns.set_theme(
    context="paper",
    style="ticks",
    palette="colorblind",
    font="serif",
    font_scale=1.2,
    color_codes=True,
    rc={"lines.linewidth": 3, "figure.figsize": (8, 6), "figure.dpi": 200},
)
sns.despine()

# Input handling
if len(sys.argv) != 2:
    raise ValueError("Usage: python plot_pathlength.py <basename>")

basename = sys.argv[1]

xml_file = basename + ".xml"
filename = basename + "_binned_tallies.csv"

df = pd.read_csv(filename)

col = df["collision_estimate"].to_numpy()
pl = df["pathlength"].to_numpy()
pos = df["position"].to_numpy()
pl_var = df["pathlength_variance"].to_numpy()
pl_std = df["pathlength_std"].to_numpy()
dxs = df["dx"].to_numpy()

relative_error = []
for std, flux in zip(pl_std, pl):
    relative_error.append(std / flux)


integrated_variance = 0
for var, dx in zip(pl_var, dxs):
    integrated_variance += var * dx

upper = pl + pl_std
lower = pl - pl_std
lower = np.clip(lower, 0.0, None)  # don't go below 0

plt.figure()
plt.plot(pos, pl, label="pathlength est.")
plt.plot(pos, col, label="collision est.")
plt.fill_between(pos, lower, upper, alpha=0.3, label=r"$\pm 1\sigma$")
plt.xlabel("x")
plt.ylabel("Normalized Flux")
plt.ylim(bottom=0)
plt.legend()
plt.title(f"Tallied results for {basename}")
plt.savefig(f"{basename}_tallies.png")
plt.close()

plt.figure()
plt.plot(pos, relative_error)
plt.title(f"Relative Error for {basename}")
plt.xlabel("x")
plt.ylabel("Relative Error")
plt.yscale("log")
plt.savefig(f"{basename}_relative_error.png")
plt.close()

print(f"Integrated variance: {integrated_variance}")
