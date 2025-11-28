import numpy as np
import matplotlib.pyplot as plt
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

# ---------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------
if len(sys.argv) != 2:
    raise ValueError("Usage: python plot_pathlength.py <basename>")

basename = sys.argv[1]

xml_file = basename + ".xml"
pathlength_file = basename + "_pl.csv"
regions_file = basename + "_regions.csv"

df = pd.read_csv(pathlength_file, header=0)
regions = pd.read_csv(regions_file, header=0)

# ---------------------------------------------------------------------
# Read XML settings (n_particles, n_generations, n_inactive)
# ---------------------------------------------------------------------
sim_data = pd.read_xml(xml_file, xpath=".//settings")

n_particles = int(sim_data["n_particles"].iloc[0])
n_generations = int(sim_data["n_generations"].iloc[0])
n_inactive = int(sim_data["n_inactive"].iloc[0])

total_active_particles = n_particles * (n_generations - n_inactive)
print("Total active particles =", total_active_particles)

# ---------------------------------------------------------------------
# Region-dependent grid construction
# ---------------------------------------------------------------------
n_cells_per_region = 50

bin_edges_list = []
bin_widths = []

for _, reg in regions.iterrows():
    xmin_r = reg["xmin"]
    xmax_r = reg["xmax"]
    sigma_t = reg["Sigma_t"]

    edges = np.linspace(xmin_r, xmax_r, n_cells_per_region + 1)
    width = (xmax_r - xmin_r) / n_cells_per_region

    bin_edges_list.append(edges)
    bin_widths.extend([width] * n_cells_per_region)

# Concatenate while removing boundary duplicates
bin_edges = np.concatenate(
    [be if i == 0 else be[1:] for i, be in enumerate(bin_edges_list)]
)

bin_widths = np.array(bin_widths)

n_bins = len(bin_widths)

# ---------------------------------------------------------------------
# Path-length binning
# ---------------------------------------------------------------------
pathlength_bins = np.zeros(n_bins)

# Loop through each flight segment
for _, row in df.iterrows():
    x0 = row["start"]
    x1 = row["end"]
    seg_length = row["pathlength"]

    # Ensure ordering (in case start > end for backtracking)
    xmin = min(x0, x1)
    xmax = max(x0, x1)

    # Find all bins the segment overlaps
    bin_i_start = np.digitize([xmin], bin_edges)[0] - 1
    bin_i_end = np.digitize([xmax], bin_edges)[0] - 1

    bin_i_start = max(bin_i_start, 0)
    bin_i_end = min(bin_i_end, n_bins - 1)

    # For each crossed bin, compute overlap fraction
    for b in range(bin_i_start, bin_i_end + 1):
        left = bin_edges[b]
        right = bin_edges[b + 1]

        overlap = max(0.0, min(xmax, right) - max(xmin, left))
        if overlap > 0:
            frac = overlap / (xmax - xmin)
            pathlength_bins[b] += frac * seg_length

# ---------------------------------------------------------------------
# Normalize: pathlength / (N_particles * Δx)
# ---------------------------------------------------------------------
normalized_pathlength = pathlength_bins / (total_active_particles * bin_widths)

# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

plt.figure(figsize=(8, 6))
sns.lineplot(x=bin_centers, y=normalized_pathlength)
plt.xlabel("x")
plt.ylabel("Normalized Path Length Density")
plt.title("Path-Length Flux Approximation")
plt.ylim(bottom=0)
plt.tight_layout()
plt.show()
