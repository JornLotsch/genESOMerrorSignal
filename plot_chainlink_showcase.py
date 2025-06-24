#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 12:18:34 2025

@author: joern
"""

# %% imports

import os
os.chdir("/home/joern/Aktuell/GenerativeESOM/08AnalyseProgramme/R/genESOMerrorSignal")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd

# %% paths
pfad_o = "/home/joern/Aktuell/GenerativeESOM/"
pfad_u1 = "08AnalyseProgramme/R/genESOMerrorSignal/"
pfad_u2 = "08AnalyseProgramme/Python/"
pfad_u3 = "09Originale/"
pfad_umx = "04Umatrix/"
pfad_umx3 = "03Umatrix/"

# %% Functions
def annotate_axes(ax, text, fontsize=18):
    ax.text(-.08, 1.03, text, transform=ax.transAxes,
            ha="center", va="center", fontsize=fontsize, color="black")

def annotate_axes_3d(ax, text, fontsize=18):
    ax.text2D(0.05, 0.9, text, transform=ax.transAxes,
              ha="left", va="top", fontsize=fontsize, color="black")

colors = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

sns.set_palette(sns.color_palette(colors))

d = dict(zip(range(1, 20), colors))

# %% Load Data
df = pd.read_csv(pfad_o + pfad_u1 + "df_analysis.csv")
df["class"] = df["class"].replace(d)

data_levels = ["Original", "n_generated_1", "n_generated_10", "n_generated_50", "n_generated_100"]

with sns.axes_style("whitegrid"):
    fig = plt.figure(figsize=(30, 12))
    gs = gridspec.GridSpec(2, 5, figure=fig, wspace=0.2, hspace=0.15)

    for row, algo in enumerate(["genESOM", "genGMM"]):
        for col, datalevel in enumerate(data_levels):
            ax = fig.add_subplot(gs[row, col], projection='3d')
            subset = df[(df["GenerationAlgorithm"] == algo) & (df["Data"] == datalevel)]
            ax.scatter(subset["X1"], subset["X2"], subset["X3"], c=subset["class"], s=15, alpha=0.7)
            ax.set_xlim(-2.4, 2.4)
            ax.set_ylim(-2.4, 2.4)
            ax.set_zlim(-2.4, 2.4)
            ax.set_xlabel("X1")
            ax.set_ylabel("X2")
            ax.set_zlabel("X3")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_zticks([])
            if row == 0:
                ax.set_title(datalevel, fontsize=16)
            if col == 0:
                annotate_axes_3d(ax, algo, fontsize=16)

    fig.suptitle(
        "Chainlink: generative ESOM (top) versus generative GMM (bottom) data generation",
        fontsize=24,
        y=0.97,
        x=0.12,
        ha='left'
    )
    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.show()
    fig.savefig("Chainlink_ESOM_vs_GMM_grid.png", format="png", bbox_inches='tight')
