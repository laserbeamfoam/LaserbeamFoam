'''
License
    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published 
    by the Free Software Foundation, either version 3 of the License, 
    or (at your option) any later version.
    
    This program is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
    
    See the GNU General Public License for more details. You should have 
    received a copy of the GNU General Public License along with this 
    program. If not, see <https://www.gnu.org/licenses/>. 
    

Purpose:
  Postprocessing script for finished simulation cases using ParaView 
  (pvpython):
  - Performs a mesh-aware continuity check of the single track.
  - If the track is continuous, computes per cross-section metrics:
    W (width), H (track height), H_depth (depth), and Porosity.

Assumptions:
  - ParaView (with pvpython) is installed and accessible in PATH.
  - Case outputs are compatible with laserbeamFoam (OpenFOAM v2412).
  - The computational mesh is regular and its spacing is known, allowing 
  expected cross-sections along the track to be enumerated.
  - Within each cross-section, cell volumes are uniform (or accounted for), so 
  that porosity computed from cell counts is equivalent to volume fraction.

Method:
  1) Continuity check (mesh-aware):
     - All expected cross-sections are determined from the known mesh spacing.
     - For each cross-section, it is verified whether at least one cell is 
     classified as material.
     - If any expected cross-section contains no material cells, the track is 
     considered non-continuous and no further analysis is performed.
     - The result is saved as a boolean file: "continuous.joblib".

  2) Per cross-section analysis (performed only if the track is continuous) – 
    two-stage workflow:

     Stage A — Row-level processing:
       - Each cross-section (x–z plane at fixed y) is processed row by row 
       (along the z-direction).
       - For each row, the total number of material cells, pore cells, and the 
       horizontal extent of material regions are recorded, among other 
       variables used in the subsequent analysis.
       - This information is collected for all rows across all cross-sections.

     Stage B — Metric computation:
       - Pores are classified as either:
         • Surface pores: located at the topmost layer of the cross-section.
         • Internal pores: located below the top surface layer.

       - **Width (W):**  
         Maximum horizontal extent of the melt pool within the cross-section.  
         If multiple disconnected material regions exist, the largest 
         horizontal span is used.

       - **Track height (H):**  
         Vertical extent of the melt pool, defined based on pore distribution:  
         • If surface pores are present: height = difference between the lowest
         material row and the row where the surface pore begins.  
         • If only internal pores are present (or none): height = difference 
         between the highest and lowest material rows.

       - **Depth (H_depth):**  
         Vertical distance between the lowest material row and the row where 
         the maximum horizontal width occurs (i.e., where the track width is
                                              largest).

       - **Porosity:**  
         Fraction of pore volume within the cross-section:  
         Porosity = number of pore cells / (number of pore cells + number of 
                                            material cells).  
         As cell volumes are uniform, this equals (pore volume / 
                                                   total cross-section volume).

     - Metrics are computed for all cross-sections and exported as a CSV file.

Outputs:
  - continuous.joblib          → Boolean indicating track continuity 
                                 (mesh-aware check).
  - cross_section_metrics.csv  → Per cross-section values of W, H, H_depth, and
  Porosity  (only if the track is continuous).

Authors
    Simon A. Rodriguez, University College Dublin (UCD)
    Petar Cosic, University College Dublin (UCD)
    Tom Flint, University of Manchester
    Philip Cardiff, University College Dublin (UCD)

'''

import os
import sys

sys.path.insert(0, os.getcwd())

from functions import terminal, calculate_geometry_full_meltpool, plotResults
from input_data import *

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


terminal(
    f'bash -c "source {OF_LOCATION} && pvpython {os.path.join(SCRIPT_DIR, "extract_meltpool.py")}"'
)
calculate_geometry_full_meltpool(CSV_3D = "./meltpool.csv")

# Plot regardless of continuity; require metrics CSV only
if PLOT_GEOMETRY_VS_Y_LOCATION and os.path.exists("./cross_sections_statistics.csv"):
    plotResults()
else:
    print("Skipping plotting: metrics CSV missing")

# Always try to produce a slice preview if available
from functions import plot_slice_preview
plot_slice_preview(slice_csv="meltpool_slice_xmid.csv", output_png="SlicePreview.png")
terminal("mkdir -p results_plots && mv *.png results_plots")

print("Geometry measurement finished.")
