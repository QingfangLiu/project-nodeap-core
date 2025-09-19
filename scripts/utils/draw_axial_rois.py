#!/usr/bin/env python3
"""
draw_axial_rois.py

This script plots spherical ROIs on an MNI152 template using Nilearn.
The resulting figure is saved to the Figs_paper folder of the project.

Requirements:
- nilearn
- numpy
- nibabel
"""

import os
import numpy as np
import nibabel as nib
from nilearn import plotting, image, datasets

# ---------------------------------------------------------------------
# Project folder path (match R Setup.R)
# ---------------------------------------------------------------------
project_folder = "/Users/liuq13/project-nodeap-core"

# Define output folder and file
use_folder = os.path.join(project_folder, "Figs_paper")
os.makedirs(use_folder, exist_ok=True)
file_path = os.path.join(use_folder, "roi_plot.png")

# ---------------------------------------------------------------------
# Load background image (e.g., MNI152 template)
# ---------------------------------------------------------------------
bg_img = datasets.load_mni152_template()

# ---------------------------------------------------------------------
# Define coordinates, radii, and colors
# ---------------------------------------------------------------------
coords = [[34, 54, -14], [28, 38, -16]]  # example MNI coordinates
radii = [6, 6]  # radius in mm (currently not used for markers)
colors = ['#ED7014', '#B8008B']  # hex color codes

# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------
display = plotting.plot_anat(
    bg_img,
    display_mode='z',
    cut_coords=[-16],
    title='',
    annotate=False,
    draw_cross=False
)

for coord, color in zip(coords, colors):
    display.add_markers([coord], marker_color=color, marker_size=100)

plotting.show()
display.savefig(file_path, dpi=1200)
print(f"ROI plot saved to: {file_path}")
