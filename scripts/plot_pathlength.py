import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from numba import njit

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

# ---------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------
if len(sys.argv) != 2:
    raise ValueError("Usage: python plot_pathlength.py <basename>")

basename = sys.argv[1]

xml_file = basename + ".xml"
pathlength_file = basename + "_pl.csv"
regions_file = basename + "_regions.csv"

df = pd.read_csv(pathlength_file)
regions = pd.read_csv(regions_file)

# ---------------------------------------------------------------------
# Read XML settings (n_particles, n_generations, n_inactive)
# ---------------------------------------------------------------------
# sim_data = pd.read_xml(xml_file, xpath=".//settings")
with open(xml_file, "rb") as f:
    sim_data = pd.read_xml(f, xpath=".//settings")

n_particles = int(sim_data["n_particles"].iloc[0])
n_generations = int(sim_data["n_generations"].iloc[0])
n_inactive = int(sim_data["n_inactive"].iloc[0])

total_active_particles = n_particles * (n_generations - n_inactive)
print("Total active particles =", total_active_particles)

# ---------------------------------------------------------------------
# Region-dependent grid construction
# ---------------------------------------------------------------------
n_cells_per_region = 25

bin_edges_list = []
bin_widths_list = []

for _, reg in regions.iterrows():
    xmin_r, xmax_r = reg["xmin"], reg["xmax"]
    edges = np.linspace(xmin_r, xmax_r, n_cells_per_region + 1)
    width = (xmax_r - xmin_r) / n_cells_per_region
    bin_edges_list.append(edges)
    bin_widths_list.extend([width] * n_cells_per_region)

# Concatenate while removing boundary duplicates
bin_edges = np.concatenate(
    [be if i == 0 else be[1:] for i, be in enumerate(bin_edges_list)]
)
bin_widths = np.array(bin_widths_list)
n_bins = len(bin_widths)

# ---------------------------------------------------------------------
# Pathlength binning (vectorized + Numba for multi-bin segments)
# ---------------------------------------------------------------------
x0 = df["start"].to_numpy()
x1 = df["end"].to_numpy()
seg_length = df["pathlength"].to_numpy()
weight = df["weight"].to_numpy()

xmin = np.minimum(x0, x1)
xmax = np.maximum(x0, x1)

bin_i_start = np.digitize(xmin, bin_edges) - 1
bin_i_end = np.digitize(xmax, bin_edges) - 1
bin_i_start = np.clip(bin_i_start, 0, n_bins - 1)
bin_i_end = np.clip(bin_i_end, 0, n_bins - 1)

pathlength_bins = np.zeros(n_bins)
squared_sum_bins = np.zeros_like(pathlength_bins)

# Handle single-bin segments fully vectorized
single_mask = bin_i_start == bin_i_end
np.add.at(
    pathlength_bins,
    bin_i_start[single_mask],
    seg_length[single_mask] * weight[single_mask],
)
np.add.at(
    squared_sum_bins,
    bin_i_start[single_mask],
    (seg_length[single_mask] * weight[single_mask]) ** 2,
)

# Multi-bin segments handled with Numba
multi_mask = ~single_mask
multi_indices = np.where(multi_mask)[0]


@njit
def bin_multi_segments(
    pathlength_bins,
    multi_indices,
    bin_i_start,
    bin_i_end,
    xmin,
    xmax,
    seg_length,
    bin_edges,
    weights,
    exp=1,
):
    for i in multi_indices:
        b_start = bin_i_start[i]
        b_end = bin_i_end[i]
        xmin_seg = xmin[i]
        xmax_seg = xmax[i]
        seg_len = seg_length[i]
        weight = weights[i]

        for b in range(b_start, b_end + 1):
            left = bin_edges[b]
            right = bin_edges[b + 1]
            overlap = max(0.0, min(xmax_seg, right) - max(xmin_seg, left))
            if overlap > 0:
                pathlength_bins[b] += (
                    overlap / (xmax_seg - xmin_seg) * (seg_len * weight) ** exp
                )


bin_multi_segments(
    squared_sum_bins,
    multi_indices,
    bin_i_start,
    bin_i_end,
    xmin,
    xmax,
    seg_length,
    bin_edges,
    weight,
    2,
)

# change exponent to 2 to get sum of squares
bin_multi_segments(
    pathlength_bins,
    multi_indices,
    bin_i_start,
    bin_i_end,
    xmin,
    xmax,
    seg_length,
    bin_edges,
    weight,
    1,
)

# ---------------------------------------------------------------------
# Normalize: pathlength / (N_particles * Δx)
# ---------------------------------------------------------------------
normalized_pathlength = pathlength_bins / (total_active_particles * bin_widths)

# ---------------------------------------------------------------------
# Compute variance
# ---------------------------------------------------------------------
normalized_ss = squared_sum_bins / (total_active_particles * bin_widths**2)
variance = (total_active_particles / (total_active_particles - 1)) * (
    normalized_ss - normalized_pathlength**2
)
integrated_variance = np.dot(variance, bin_widths)  # approximate midpoint rule

# ---------------------------------------------------------------------
# Compute relative error
# ---------------------------------------------------------------------
sigma = np.sqrt(variance)
relative_error = sigma / normalized_pathlength

# ---------------------------------------------------------------------
# Plot with uncertainty band
# ---------------------------------------------------------------------
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

lower = normalized_pathlength - sigma
upper = normalized_pathlength + sigma
lower = np.clip(lower, 0.0, None)

plt.figure(figsize=(8, 6))

# Mean flux
plt.plot(
    bin_centers,
    normalized_pathlength,
    label="Flux",
    linewidth=2.5,
)

# Uncertainty band (±1σ)
plt.fill_between(
    bin_centers,
    lower,
    upper,
    alpha=0.3,
    label=r"$\pm 1\sigma$",
)

plt.xlabel("x")
plt.ylabel("Normalized Flux")
plt.title("Path-Length Flux Tally with Statistical Uncertainty")
plt.ylim(bottom=0)
plt.legend()
plt.tight_layout()
plt.savefig(f"{basename}_pathlength_plot.png")
plt.close()

print(f"Integrated variance: {integrated_variance}")
