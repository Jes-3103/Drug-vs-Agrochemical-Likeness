# figures.py
# Draws the two figures from results/screening_results.csv and saves them
# as PNGs in results/. Uses matplotlib only, so it needs no RDKit.

import os
import matplotlib
matplotlib.use("Agg")          # write to a file, don't open a window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use("dark_background")

os.makedirs("results", exist_ok=True)
df = pd.read_csv("results/screening_results.csv")


# --- Figure 1: chemical space --------------------------------------------

fig, ax = plt.subplots(figsize=(7, 5))
for cat, marker in [("drug", "o"), ("agrochemical", "s")]:
    sub = df[df["category"] == cat]
    ax.scatter(sub["MW"], sub["LogP"], marker=marker, label=cat, s=60)
ax.axvline(500, ls="--", c="grey")   # Ro5 MW limit
ax.axhline(5, ls="--", c="grey")     # Ro5 logP limit
ax.set_xlabel("Molecular weight (g/mol)")
ax.set_ylabel("cLogP (Crippen)")
ax.set_title("Chemical space: drugs vs. agrochemicals")
ax.legend()
fig.tight_layout()
fig.savefig("results/mw_vs_logp.png", dpi=200, facecolor=fig.get_facecolor())
print("saved results/mw_vs_logp.png")


# --- Figure 2: pass rates by class and rule set ---------------------------

rates = df.groupby("category")[["Lipinski_pass", "Tice_herb_pass",
                                "Tice_ins_pass"]].mean()

rules = ["Lipinski_pass", "Tice_herb_pass", "Tice_ins_pass"]
pretty = ["Lipinski\n(drug rule)", "Tice herbicide\n(agro rule)",
          "Tice insecticide\n(agro rule)"]

x = np.arange(len(rules)) * 1.3   # group centres
w = 0.3                           # bar width
gap = 0.03                        # space between the two bars in a pair

fig, ax = plt.subplots(figsize=(9, 5.5))
for i, (cat, col) in enumerate(zip(["drug", "agrochemical"],
                                   ["#4da6ff", "#5fd38d"])):
    vals = [rates.loc[cat, r] * 100 for r in rules]
    n = (df["category"] == cat).sum()
    bars = ax.bar(x + (i - 0.5) * (w + gap), vals, width=w,
                  label=f"{cat}s (n={n})", color=col, edgecolor="white", lw=0.6)
    ax.bar_label(bars, fmt="%.0f%%", padding=3, fontsize=9)

ax.set_xticks(x)
ax.set_xticklabels(pretty)
ax.set_ylim(0, 115)
ax.set_ylabel("Compounds passing (%)")
ax.set_title("Rule pass rates by compound class")
ax.grid(axis="y", ls=":", alpha=0.3)
ax.set_axisbelow(True)
ax.legend()
fig.tight_layout()
fig.savefig("results/pass_rates.png", dpi=200, facecolor=fig.get_facecolor())
print("saved results/pass_rates.png")
