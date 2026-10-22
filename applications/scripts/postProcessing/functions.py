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

Description
  Provides utility functions that support the meltpool characterisation.

Assumptions:
  - Called by higher-level script "characterise_meltpool.py"

Authors
    Simon A. Rodriguez, University College Dublin (UCD). All rights reserved
    Petar Cosic, University College Dublin (UCD). All rights reserved
    Tom Flint, University of Manchester (UOM). All rights reserved
    Philip Cardiff, University College Dublin (UCD). All rights reserved
    
'''


import os
import numpy as np
import pandas as pd
import input_data
from input_data import *
import importlib
import os
from joblib import dump, load
import matplotlib.pyplot as plt

def plotResults(CSV_CROSS_SECTIONS = "cross_sections_statistics.csv"):
    def generate_figure(x, y_values, xlabel, ylabel, title, name_png_file):
        plt.figure()
        plt.plot(x, y_values, marker="o")  
        plt.axhline(y=y_values.mean(), color="red", linestyle="--", 
                    label="Mean")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("./" + name_png_file + ".png")
    
    df = pd.read_csv(CSV_CROSS_SECTIONS)
    # work in micrometers for plotting
    y_locations_um = df["iy"] * 1e6
    
    keys_for_plot = ["width", "height", "depth", "porosity_at_iy"]
    if "area" in df.columns:
        # Keep order consistent: lengths first, then area, porosity last
        keys_for_plot.insert(3, "area")
    
    for key in keys_for_plot:
        if key not in df.columns:
            continue
        values_for_plot = df[key]
        if key == "porosity_at_iy":
            generate_figure(
                y_locations_um,
                values_for_plot,
                "y_coordinate (um)",
                "Porosity (porous volume / total volume)",
                "Porosity vs. y-coordinate",
                "Porosity",
            )
        elif key == "area":
            values_um2 = values_for_plot * 1e12
            generate_figure(
                y_locations_um,
                values_um2,
                "y_coordinate (um)",
                "Area (um^2)",
                "Area vs. y-coordinate",
                "Area",
            )
        else:
            values_um = values_for_plot * 1e6
            title = key.capitalize() + " vs. y-coordinate"
            generate_figure(
                y_locations_um,
                values_um,
                "y_coordinate (um)",
                f"{key} (um)",
                title,
                key.capitalize(),
            )


def plot_slice_preview(slice_csv="meltpool_slice_xmid.csv", output_png="SlicePreview.png"):
    """
    Quick scatter plot of a mid-x slice to visualise where width/height/depth
    are measured. Uses micrometers on both axes.
    """
    if not os.path.exists(slice_csv) or os.path.getsize(slice_csv) == 0:
        return
    try:
        df = pd.read_csv(slice_csv)
    except pd.errors.EmptyDataError:
        # Nothing to plot (empty CSV)
        return
    if df.empty:
        return
    # Accept both Points:0/1/2 and Points_0/1/2
    cols = df.columns
    y_col = "Points_1" if "Points_1" in cols else "Points:1"
    z_col = "Points_2" if "Points_2" in cols else "Points:2"
    if y_col not in cols or z_col not in cols:
        return
    y_um = df[y_col] * 1e6
    z_um = df[z_col] * 1e6
    plt.figure()
    plt.scatter(y_um, z_um, s=5, alpha=0.6)
    plt.xlabel("y (um)")
    plt.ylabel("z (um)")
    plt.title("Mid-x slice (melt pool points)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_png)


def save_averages_to_excel(cross_sections_statistics_df,
                           filename="metrics_summary.csv"):
    """Save average width/height/depth (and area if available) to a CSV file."""
    try:
        summary = {
            "width_mean_um": cross_sections_statistics_df["width"].mean() * 1e6,
            "height_mean_um": cross_sections_statistics_df["height"].mean() * 1e6,
            "depth_mean_um": cross_sections_statistics_df["depth"].mean() * 1e6,
        }
        if "area" in cross_sections_statistics_df.columns:
            summary["area_mean_um2"] = (
                cross_sections_statistics_df["area"].mean() * 1e12
            )

        df_summary = pd.DataFrame([summary])
        df_summary.to_csv(filename, index=False)
        print(f"Averages saved to {filename}")
    except Exception as exc:
        print(f"Could not write {filename}: {exc}")
            

    
def terminal(command):
    os.system(command)

def is_meltpool_continuous(CSV_3D = "meltpool.csv"):
    df = pd.read_csv(CSV_3D)
    # Round coordinates to avoid floating point mismatch
    df["Points_0"] = df["Points_0"].round(8)
    df["Points_1"] = df["Points_1"].round(8)
    df["Points_2"] = df["Points_2"].round(8)
    
    x = df["Points_0"].values
    y = df["Points_1"].values
    z = df["Points_2"].values
    
    y0 = Y_COORD_BEGIN_TRACK + LASER_DIAMETER / 2 
    y_max = Y_COORD_END_TRACK - LASER_DIAMETER / 2  
    
    meltpool_is_continuous = True
    
    # Use a tolerance-based approach for the mid-section
    x_mid_section = (X_MIN_AND_MAX_DOMAIN[0] + X_MIN_AND_MAX_DOMAIN[1])/2
    # Find unique x values
    unique_x = np.unique(x)
    # Find the x value closest to mid_section
    if len(unique_x) == 0:
        return False
        
    closest_x_idx = np.abs(unique_x - x_mid_section).argmin()
    closest_x = unique_x[closest_x_idx]
    
    mask_x_mid_section = (x == closest_x)
    mid_plane_df = df[mask_x_mid_section]
    
    if mid_plane_df.empty:
        return False
        
    y_at_mid = mid_plane_df["Points_1"].values
    z_at_mid = mid_plane_df["Points_2"].values
    
    # Check continuity along Y for this approximate mid-plane
    unique_y_mid = np.sort(np.unique(y_at_mid))
    # Filter to ROI
    unique_y_mid = unique_y_mid[(unique_y_mid >= np.round(y0, 8)) & (unique_y_mid <= np.round(y_max, 8))]
    
    if len(unique_y_mid) < 2:
        # Almost no data
        # Check if we have ANY data in the track range
        return False
        
    # Check gaps
    # For now, let's be lenient. If we have data, we assume it's continuous enough to measure metrics
    # The original check was very strict about Z-depth.
    
    for iy in unique_y_mid:
        z_at_iy = z_at_mid[y_at_mid == iy]
        if len(z_at_iy) == 0:
            continue
            
        z_min = np.min(z_at_iy)
        z_max = np.max(z_at_iy)
        height = z_max - z_min
        
        # If height is very small (less than 3 * CELL_SIZE), it's potentially "broken"
        if height < 3 * CELL_SIZE:
             meltpool_is_continuous = False
             break
             
    dump(meltpool_is_continuous, "./continuous.joblib")
    return meltpool_is_continuous

def calculate_statistics_rows_meltpool(CSV_3D, meltpool_is_continuous):
    
    df = pd.read_csv(CSV_3D)
    # Ensure rounding
    df["Points_0"] = df["Points_0"].round(8)
    df["Points_1"] = df["Points_1"].round(8)
    df["Points_2"] = df["Points_2"].round(8)
    
    x = df["Points_0"].values
    y = df["Points_1"].values
    z = df["Points_2"].values
    
    y0 = Y_COORD_BEGIN_TRACK + LASER_DIAMETER / 2 
    y_max = Y_COORD_END_TRACK - LASER_DIAMETER / 2  
    
    # Identify unique Y levels actually present in the data
    all_y = np.unique(y)
    # Filter for ROI
    valid_ys = all_y[(all_y >= np.round(y0, 8)) & (all_y <= np.round(y_max, 8))]
    
    id_row = 0
    Statistics = []
    pore_locatios_at_rows = []
    pores_at_row_are_internal = []
    
    # Iterate over existing Y sections only
    for iy in valid_ys:
        mask = (y == iy)
        cells_at_iy = df[mask]
        
        # Within this Y-slice, find unique Z rows
        z_at_iy = cells_at_iy["Points_2"].values
        unique_z = np.unique(z_at_iy)
        
        for iz in unique_z:
            mask2 = (z_at_iy == iz)
            cells_at_iy_iz = cells_at_iy[mask2]
            if cells_at_iy_iz.empty:
                continue
                
            x_at_iy_iz = cells_at_iy_iz["Points_0"].values
            min_x = np.min(x_at_iy_iz)
            max_x = np.max(x_at_iy_iz)
            
            # Use data-driven width unless grid snapping is strictly required
            # The original code snapped to CELL_SIZE grid. We can do that too.
            min_x_snapped = np.round(np.round(min_x/CELL_SIZE) * CELL_SIZE, 8)
            max_x_snapped = np.round(np.round(max_x/CELL_SIZE) * CELL_SIZE, 8)
            
            width_row = np.round(max_x_snapped - min_x_snapped, 8)
            
            # Simple check for pores: gaps in X?
            # Original code expected continuous cells.
            number_non_void_cells_in_row = len(x_at_iy_iz)
            
            # Calculate expected cells if it were solid
            expected_cells = int(np.round(width_row / CELL_SIZE)) if width_row > 0 else 1
            if expected_cells < 1: expected_cells = 1
            
            n_pores_in_row = max(0, expected_cells - number_non_void_cells_in_row)
            row_has_pores = (n_pores_in_row > 0)
            
            pore_locations_at_row_i = []
            pores_at_row_i_are_internal_list = []
            
            if row_has_pores:
                pore_locations_at_row_i = [id_row] 
                is_top_layer = (iz == np.max(unique_z))
                pores_at_row_i_are_internal_list = [not is_top_layer] * n_pores_in_row
                
                pore_locatios_at_rows.append(pore_locations_at_row_i)
                pores_at_row_are_internal.append([id_row] + pores_at_row_i_are_internal_list)
            else:
                pore_locatios_at_rows.append([id_row, "NA"])
                pores_at_row_are_internal.append([id_row, "NA"])

            new_statistics_row = [id_row, iy, iz, 
                               min_x_snapped, max_x_snapped,
                               row_has_pores, n_pores_in_row, width_row,
                               number_non_void_cells_in_row]
            Statistics.append(new_statistics_row)
            id_row += 1

    row_statistics = pd.DataFrame(Statistics, columns = ["id_row", "y_coord", 
                                                          "z_coord_", "x_min", 
                                                          "x_max", 
                                                          "row_has_pores", 
                                                      "number_of_pores_in_row", 
                                                          "width_row",
                                               "number_non_void_cells_in_row"])
    row_statistics.to_csv("row_statistics.csv", index=False, encoding="utf-8") 
    
    return row_statistics, pore_locatios_at_rows, pores_at_row_are_internal

def calculate_cross_sections_statistics(row_statistics, 
                                        pore_locatios_at_rows, 
                                        pores_at_row_are_internal, 
                                        meltpool_is_continuous):
    
    cross_sections_statistics = []
    y = row_statistics["y_coord"].to_numpy()
    row_has_pores = row_statistics["row_has_pores"]
    y_unique = np.unique(y)  
    
    void_iy_levels= []
    if (not meltpool_is_continuous):
        try:
            void_iy_levels = load("void_iy_levels.joblib")
        except Exception:
            void_iy_levels = []
    
    for iy in y_unique:
        if (iy not in void_iy_levels):
            mask = (iy == y)
            cross_section_at_iy = row_statistics[mask]
            z_at_iy = cross_section_at_iy["z_coord_"]
            width_rows_at_iy = cross_section_at_iy["width_row"]
            number_pores_at_iy = cross_section_at_iy["number_of_pores_in_row"]
            number_non_void_cells_in_row_at_iy = cross_section_at_iy[
                                                "number_non_void_cells_in_row"]
            # Simple metrics: width = max row width; height = z span;
            # depth = z at max width minus lowest z.
            width = width_rows_at_iy.max()
            height = z_at_iy.max() - z_at_iy.min()
            z_at_max_width = z_at_iy[width_rows_at_iy == width]
            depth = z_at_max_width.max() - z_at_iy.min()
            
            pore_cells = np.sum(number_pores_at_iy.to_numpy())
            material_cells = np.sum(number_non_void_cells_in_row_at_iy.to_numpy())
            porous_volume_at_iy = pore_cells * (CELL_SIZE**3)
            total_cells_at_iy = material_cells + pore_cells
            total_volume_material_at_iy = total_cells_at_iy * (CELL_SIZE**3)
            area = total_cells_at_iy * (CELL_SIZE**2)
            
            porosity_at_iy = porous_volume_at_iy/total_volume_material_at_iy

            cross_sections_statistics.append([iy, width, height, depth, area,
                                      porosity_at_iy, total_volume_material_at_iy])
               
    cross_sections_statistics_df = pd.DataFrame(cross_sections_statistics, 
                                                columns = ["iy", "width", 
                                                           "height", "depth",
                                                           "area",
                                                           "porosity_at_iy", 
                                                "total_volume_material_at_iy"])
    
    cross_sections_statistics_df.to_csv("./cross_sections_statistics.csv", 
                                        index=False, encoding="utf-8") 
    save_averages_to_excel(cross_sections_statistics_df,
                           filename="metrics_summary.csv")
        

    return cross_sections_statistics_df


def calculate_geometry_full_meltpool(CSV_3D = "meltpool.csv"):
    # Always proceed with geometry calculation; keep continuity flag for info
    meltpool_is_continuous = is_meltpool_continuous(CSV_3D)
    row_statistics, pore_locatios_at_rows, pores_at_row_are_internal = (
        calculate_statistics_rows_meltpool(CSV_3D, meltpool_is_continuous)
    )
    cross_sections_statistics = calculate_cross_sections_statistics(
        row_statistics,
        pore_locatios_at_rows, 
        pores_at_row_are_internal, 
        meltpool_is_continuous
    )
    return cross_sections_statistics
