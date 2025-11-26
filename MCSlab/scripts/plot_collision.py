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
    raise ValueError("Usage: python plot_collision.py <basename>")

basename = sys.argv[1]  # e.g. "run1"
xml_file = basename + ".xml"  # e.g. "input.xml"

# input collision file and region file
collisions_file = basename + "_col.csv"
regions_file = basename + "_regions.csv"

df = pd.read_csv(collisions_file, header=0)
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
# Bin geometry
# ---------------------------------------------------------------------
num_regions = len(regions)
positions = df["position"].to_numpy()

n_cells_per_region = 50  # adjustable
n_bins = n_cells_per_region * num_regions

xmin = regions["xmin"].iloc[0]
xmax = regions["xmax"].iloc[-1]
bin_edges = np.linspace(xmin, xmax, n_bins + 1)

# digitize positions → bin indices in [1, n_bins]
bin_indices = np.digitize(positions, bin_edges) - 1  # shift to [0, n_bins-1]

# initialize bins
flux_bins = np.zeros(n_bins)

# count collisions
for idx in bin_indices:
    if 0 <= idx < n_bins:
        flux_bins[idx] += 1

# ---------------------------------------------------------------------
# Normalize: divide by total Σ_t in the region containing the bin
# ---------------------------------------------------------------------

# Build an array giving Σ_t for each bin
sigma_t_bins = np.zeros(n_bins)

for i, reg in regions.iterrows():
    reg_xmin = reg["xmin"]
    reg_xmax = reg["xmax"]
    sigma_t = reg["Sigma_t"]

    # compute which bin indices belong to this region
    left_bin = int((reg_xmin - xmin) / (xmax - xmin) * n_bins)
    right_bin = int((reg_xmax - xmin) / (xmax - xmin) * n_bins)

    # clip to array limits
    left_bin = max(left_bin, 0)
    right_bin = min(right_bin, n_bins)

    sigma_t_bins[left_bin:right_bin] = sigma_t

# normalized flux = collisions / (Σ_t * total active particles)
normalized_flux = flux_bins / (sigma_t_bins * total_active_particles)

# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

plt.figure(figsize=(8, 6))
sns.lineplot(x=bin_centers, y=normalized_flux)
plt.xlabel("x")
plt.ylabel("Normalized Flux")
plt.title("Collision Density / Σ_t")
plt.tight_layout()
plt.show()
