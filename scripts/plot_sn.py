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

forward_file = "SN_output_forward_flux_" + basename + ".csv"
adjoint_file = "SN_output_adjoint_flux_" + basename + ".csv"

dff = pd.read_csv(forward_file)
dfa = pd.read_csv(adjoint_file)

forward_pos = dff["position"].to_numpy()
adjoint_pos = dfa["position"].to_numpy()

foward_flux = dff["scalar_flux"].to_numpy()
adjoint_flux = dfa["scalar_flux"].to_numpy()
importance_map = 1 / adjoint_flux
importance_map = importance_map / max(importance_map)

plt.figure(figsize=(8, 6))

plt.plot(
    forward_pos,
    foward_flux,
    label="Forward Flux",
    linewidth=2.5,
)
plt.plot(
    adjoint_pos,
    adjoint_flux,
    label="Adjoint Flux",
    linewidth=2.5,
)
plt.plot(
    adjoint_pos,
    importance_map,
    label="Normalized Importance Map",
    linewidth=2.5,
)

plt.xlabel("x")
plt.ylabel("Normalized Flux")
plt.title(f"Discrete Ordinates Results for {basename}")
plt.ylim(bottom=0)
plt.legend()
plt.tight_layout()
plt.savefig(f"{basename}_SN_plot.png")
plt.close()
