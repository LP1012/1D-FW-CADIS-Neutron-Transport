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
    rc={"lines.linewidth": 3, "figure.figsize": (8, 6), "figure.dpi": 200},
)
sns.despine()

# ---------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------
if len(sys.argv) != 2:
    raise ValueError("Usage: python plot_collision.py <basename>")

basename = sys.argv[1]  # e.g. "run1"

xml_file = basename + ".xml"
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
# Region-dependent grid construction
# ---------------------------------------------------------------------
n_cells_per_region = 25  # adjustable

bin_edges_list = []
sigma_t_bins = []
bin_widths = []

for _, reg in regions.iterrows():
    xmin_r = reg["xmin"]
    xmax_r = reg["xmax"]
    sigma_t = reg["Sigma_t"]

    # region-specific equal spacing
    edges = np.linspace(xmin_r, xmax_r, n_cells_per_region + 1)
    width = (xmax_r - xmin_r) / n_cells_per_region

    bin_edges_list.append(edges)
    bin_widths.extend([width] * n_cells_per_region)
    sigma_t_bins.extend([sigma_t] * n_cells_per_region)

# concatenate sub-edges while removing duplicates at boundaries
bin_edges = np.concatenate(
    [be if i == 0 else be[1:] for i, be in enumerate(bin_edges_list)]
)

bin_widths = np.array(bin_widths)
sigma_t_bins = np.array(sigma_t_bins)

n_bins = len(sigma_t_bins)
positions = df["position"].to_numpy()
weights = df["weight"].to_numpy()

# ---------------------------------------------------------------------
# Digitize collision positions into the region-aware bins
# ---------------------------------------------------------------------
bin_indices = np.digitize(positions, bin_edges) - 1

flux_bins = np.zeros(n_bins)
for idx in bin_indices:
    if 0 <= idx < n_bins:
        flux_bins[idx] += 1

# ---------------------------------------------------------------------
# Normalize flux: collisions / (total active particles * Δx)
# ---------------------------------------------------------------------
normalized_flux = np.zeros_like(flux_bins, dtype=float)
nonzero_mask = sigma_t_bins > 0
normalized_flux[nonzero_mask] = flux_bins[nonzero_mask] / (
    sigma_t_bins[nonzero_mask] * total_active_particles * bin_widths[nonzero_mask]
)
# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

plt.figure(figsize=(8, 6))
sns.lineplot(x=bin_centers, y=normalized_flux)
plt.xlabel("x")
plt.ylabel("Normalized Flux")
plt.title("Collision Estimator Flux Tally")
plt.ylim(bottom=0)
plt.tight_layout()
plt.savefig(f"{basename}_collision_plot.png")
