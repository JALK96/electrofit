import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Polygon

sns.set_context("talk", font_scale=0.7)
# Parameters from the book excerpt
mu = -0.1  # mean of ΔU (in kJ/mol)
sigma2 = 0.18**2  # variance
sigma = np.sqrt(sigma2)
beta = 3.0  # in kJ/mol^-1 units, doesn't matter for relative plot

# Define ΔU range
x = np.linspace(-1.0, 0.6, 500)

# P₀(ΔU) ~ Gaussian (not normalized to area = 1, just peak height)
P0 = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-((x - mu) ** 2) / (2 * sigma2))
P0 /= np.trapz(P0, x)

# exp(-β ΔU)
exp_beta_dU = np.exp(-beta * x)
exp_beta_dU_norm = exp_beta_dU / np.trapz(exp_beta_dU, x)

# Product = P₀(ΔU) * exp(-β ΔU), normalized to area = 1
product = P0 * exp_beta_dU_norm * 3
# product /= np.trapz(product, x)

new_x = np.linspace(-1.3, 0.7, 500)
# Create the plot
fig, ax = plt.subplots(figsize=(6, 4))

# Plot the curves
ax.plot(x, P0, "k--", label=r"$P_0(\Delta U)$")
ax.plot(
    x, exp_beta_dU_norm, "0.6", linestyle="dashed", label=r"$\exp(-\beta \Delta U)$"
)
ax.plot(new_x, product, "k-", label=r"$P(\Delta U) \times \exp(-\beta \Delta U)$")

# Highlight the low-ΔU region under the product curve
x_fill = new_x[new_x < -0.6]
y_fill = product[new_x < -0.6]
verts = [(x_fill[0], 0)] + list(zip(x_fill, y_fill)) + [(x_fill[-1], 0)]
poly = Polygon(verts, facecolor="none", hatch="///", edgecolor="k")
ax.add_patch(poly)

# Add labels near each curve
ax.annotate(r"$P(\Delta U_{ij})$", xy=(0.05, 1.7), fontsize=14)
ax.annotate(r"$\exp(-\beta \Delta U_{ij})$", xy=(-0.8, 2.0), fontsize=14)
ax.annotate(
    r"$P(\Delta U) \times$" + "\n" + r"$\exp(-\beta \Delta U_{ij})$",
    xy=(-0.98, 0.7),
    fontsize=14,
)

# Axes labels
ax.set_xlabel(r"$\Delta U$")
ax.set_ylabel(r"$P(\Delta U)$")
ax.set_xlim([-1.0, 0.6])
ax.set_ylim([0.0, 2.4])

plt.tight_layout()
plt.savefig("plot.pdf", transparent=True)
