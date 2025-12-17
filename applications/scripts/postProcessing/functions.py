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
                           filename="metrics_summary.xlsx"):
    """Save average width/height/depth (and area if available) to an Excel file."""
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
        df_summary.to_excel(filename, index=False)
        print(f"Averages saved to {filename}")
    except Exception as exc:
        print(f"Could not write {filename}: {exc}")
            

    
def terminal(command):
    os.system(command)

def is_meltpool_continuous(CSV_3D = "meltpool.csv"):
    df = pd.read_csv(CSV_3D)
    x = df["Points_0"].to_numpy()
    y = df["Points_1"].to_numpy()
    z = df["Points_2"].to_numpy()
    y0 = Y_COORD_BEGIN_TRACK + LASER_DIAMETER / 2 
    y_max = Y_COORD_END_TRACK - LASER_DIAMETER / 2  
    factor = 1000
    TOL = CELL_SIZE/factor
    meltpool_is_continuous = True
    # Build the y-levels
    y_levels = []
    y_level = np.round(y0, 8)
    # First, check if a y-z slice at x = x_mid_section is continuous
    x_mid_section = (X_MIN_AND_MAX_DOMAIN[0] + X_MIN_AND_MAX_DOMAIN[1])/2
    mask_x_mid_section = (x == x_mid_section)
    mid_plane_x = df[mask_x_mid_section]
    y_at_mid_plane_x = mid_plane_x["Points_1"]
    z_at_mid_plane_x = mid_plane_x["Points_2"]
    # This loop answers the question: are there cells with less than 3 cells
    # (or 4 points) aligned along the z-axis for every y-location in the 
    # x_mid_plane?
    while (y_level <= y_max):
        # Check if there are cells at every y_level
        mask_y_section_at_y_level_and_x_mid_plane = np.isclose(np.round(
                                                 y_level, 8), y_at_mid_plane_x)
        
        z_at_y_level_and_x_mid_plane = z_at_mid_plane_x[
                                     mask_y_section_at_y_level_and_x_mid_plane]
        
        
        z_min_at_y_level_and_x_mid_plane = np.min(z_at_y_level_and_x_mid_plane)
        z_max_at_y_level_and_x_mid_plane = np.max(z_at_y_level_and_x_mid_plane)
        z_min_at_iy_that_is_in_original_mesh = np.round(np.round(
                    z_min_at_y_level_and_x_mid_plane/CELL_SIZE) * CELL_SIZE, 8)
        z_max_at_iy_that_is_in_original_mesh = np.round(np.round(
                    z_max_at_y_level_and_x_mid_plane/CELL_SIZE) * CELL_SIZE, 8)
    
    
        z0_at_y_level = z_min_at_iy_that_is_in_original_mesh.copy()
        expected_levels_count = np.round((z_max_at_iy_that_is_in_original_mesh 
                           - z_min_at_iy_that_is_in_original_mesh) / CELL_SIZE)
        
        levels_count = 0
        while (z0_at_y_level <= z_max_at_iy_that_is_in_original_mesh):
            mask = (np.round(z0_at_y_level, 8) == np.round(
                                              z_at_y_level_and_x_mid_plane, 8))    
            if (np.sum(mask) >= 1):
                levels_count = levels_count + 1            
            z0_at_y_level = np.round(z0_at_y_level + CELL_SIZE, 8)
        
        if (levels_count < 4): # 4
            meltpool_is_continuous = False
            break
            
        if (meltpool_is_continuous):            
            y_level = np.round(y_level + CELL_SIZE, 8)
    
        elif (y_level < y_max):
            y_level =  2 * y_max
            
    dump(meltpool_is_continuous, "./continuous.joblib")

    if (not meltpool_is_continuous): 
     # Build the y-levels
        y_levels = []
        y_level = np.round(y0, 8)
        # This loop answers the question: "Are there y_levels that are  
        # completely void. If so, where are they?
        while (y_level <= y_max):
            y_levels.append(y_level)
            y_level = np.round(y_level + np.round(CELL_SIZE, 8), 8)
        
        are_there_material_cells_at_iy = True
        void_iy_levels = []
        
        for iy in y_levels:
            mask = (iy == np.round(y, 8))
            if (np.sum(mask) == 0):
                are_there_material_cells_at_iy = False
                void_iy_levels.append(iy)
    
        if (len(void_iy_levels) > 0):
            dump(void_iy_levels, "./void_iy_levels.joblib")
    
        if (len(void_iy_levels) > 0):
            dump(void_iy_levels, "./void_iy_levels.joblib")
    
    return meltpool_is_continuous

def calculate_statistics_rows_meltpool(CSV_3D, meltpool_is_continuous):
    
    # This function takes the .csv file tht represents the meltpool and 
    # calculates several metrics for every row in the meltpool. 
    # The resulting objects are: 
    # 1. Statistics = [id_row, y_coord, z_coord_ x_min, x_max, 
    # row_has_pores, number_of_pores_in_row, width_row, 
    # number_non_void_cells_in_row]
    # 2. pore_locatios_at_rows = A file with "NA" if the row has no pores. 
    # Otherwise, it has list with the x_coord of every pore at the row
    # 3. pores_at_row_are_internal = A file with "NA" if the row has no pores.
    # Otherwise, it has "True" or "False" for every pore in the row. True if 
    # the pore is internal, Flase, otherwise
    
    df = pd.read_csv(CSV_3D)
    x = df["Points_0"].to_numpy()
    y = df["Points_1"].to_numpy()
    z = df["Points_2"].to_numpy()
    
    void_iy_levels= []
    if (not meltpool_is_continuous):
        try:
            void_iy_levels = load("void_iy_levels.joblib")
        except Exception:
            void_iy_levels = []
    
    y0 = Y_COORD_BEGIN_TRACK + LASER_DIAMETER / 2 
    y_max = Y_COORD_END_TRACK - LASER_DIAMETER / 2  
    factor = 1000
    TOL = CELL_SIZE/factor
    iy = y0
    id_row = 0
    Statistics = []
    pore_locatios_at_rows = []
    pores_at_row_are_internal = []
    
    # Iterate over all the y-sections
    while (iy <= y_max):
                
        if(iy not in void_iy_levels):
            mask = (iy == np.round(y, 8))
            cells_at_iy = df[mask]
            if cells_at_iy.empty:
                # no cells at this y, skip safely
                iy = np.round(iy + CELL_SIZE, 8)
                continue
            x_at_iy = cells_at_iy["Points_0"].to_numpy()
            y_at_iy = cells_at_iy["Points_1"].to_numpy()
            z_at_iy = cells_at_iy["Points_2"].to_numpy()
            
            z_min_at_iy = np.min(z_at_iy)
            z_max_at_iy = np.max(z_at_iy)
            x_min_at_iy = np.min(x_at_iy) 
            x_max_at_iy = np.max(x_at_iy) 
            z_min_at_iy_that_is_in_original_mesh = np.round(np.round(
                                         z_min_at_iy/CELL_SIZE) * CELL_SIZE, 8)
            z_max_at_iy_that_is_in_original_mesh = np.round(np.round(
                                         z_max_at_iy/CELL_SIZE) * CELL_SIZE, 8)
            x_min_at_iy_that_is_in_original_mesh = np.round(np.round(
                                         x_min_at_iy/CELL_SIZE) * CELL_SIZE, 8)
            x_max_at_iy_that_is_in_original_mesh = np.round(np.round(
                                         x_max_at_iy/CELL_SIZE) * CELL_SIZE, 8)
            
            iz = z_min_at_iy_that_is_in_original_mesh
            
            while (iz <= z_max_at_iy_that_is_in_original_mesh):
                mask2 = ((iz == np.round(z_at_iy, 8)))
                cells_at_iy_iz = cells_at_iy[mask2]
                x_at_iy_iz = cells_at_iy_iz["Points_0"].to_numpy()
                if (x_at_iy_iz.shape[0] == 0):
                    iz = z_max_at_iy_that_is_in_original_mesh
                    iz = np.round(iz + CELL_SIZE, 8) 
                              
                else:
                    min_x_at_iy_iz = np.min(x_at_iy_iz)
                    max_x_at_iy_iz = np.max(x_at_iy_iz)
                    
                    
                    min_x_at_iy_iz_that_is_in_original_mesh = np.round(
                             np.round(min_x_at_iy_iz/CELL_SIZE) * CELL_SIZE, 8)
                    max_x_at_iy_iz_that_is_in_original_mesh = np.round(
                             np.round(max_x_at_iy_iz/CELL_SIZE) * CELL_SIZE, 8)
                    
                    distance_minx_max_at_zlevel = np.round(
                        max_x_at_iy_iz_that_is_in_original_mesh - 
                        min_x_at_iy_iz_that_is_in_original_mesh, 8)
                    expected_number_cells_at_iy_iz = int(
                                         distance_minx_max_at_zlevel/CELL_SIZE)
                   
                    ix = min_x_at_iy_iz_that_is_in_original_mesh
                    init_in = ix
                    number_non_void_cells_in_row = 0
                    
                    row_has_pores = False
                    n_pores_in_row = 0
                    width_row = np.round(
                                      max_x_at_iy_iz_that_is_in_original_mesh - 
                                       min_x_at_iy_iz_that_is_in_original_mesh, 
                                         8)
                    
                    pore_locations_at_row_i = []
                    pores_at_row_i_are_internal = []
                    if (expected_number_cells_at_iy_iz > 1):
                        while (ix < max_x_at_iy_iz_that_is_in_original_mesh):
                            if (np.sum([ix == np.round(x_at_iy_iz, 8)]) > 0):
                                cell_is_a_pore = False
                                number_non_void_cells_in_row = (
                                                  number_non_void_cells_in_row 
                                                  + 1)
                            else:
                                cell_is_a_pore = True
                                n_pores_in_row = n_pores_in_row + 1
                                row_has_pores = True
                                pore_locations_at_row_i.append(ix)
                                mask3 = (ix == np.round(x_at_iy, 8))
                                cells_at_ix_iy = cells_at_iy[mask3]
                                z_at_ix_iy = cells_at_ix_iy["Points_2"]
                                pores_at_row_i_are_internal.append(
                                                   np.sum(iz < z_at_ix_iy) > 0)
                                                        
                            ix = np.round(ix + CELL_SIZE, 8)
                            
                    if (expected_number_cells_at_iy_iz == 1):
                        number_non_void_cells_in_row = 1
                    
                    if (n_pores_in_row > 0):
                        pore_locations_at_row_i.insert(0, id_row)
                        pore_locatios_at_rows.append(pore_locations_at_row_i)
                        pores_at_row_i_are_internal.insert(0, id_row)
                        pores_at_row_are_internal.append(
                                                   pores_at_row_i_are_internal)
                    else:
                        pore_locatios_at_rows.append([id_row, "NA"])
                        pores_at_row_are_internal.append([id_row, "NA"])
                    
                    new_statistics_row = [id_row, iy, iz, 
                                       min_x_at_iy_iz_that_is_in_original_mesh, 
                                       max_x_at_iy_iz_that_is_in_original_mesh, 
                                      row_has_pores, n_pores_in_row, width_row,
                                                  number_non_void_cells_in_row]
                    
                    Statistics.append(new_statistics_row)
   
                    iz = np.round(iz + CELL_SIZE, 8)
                    id_row = id_row + 1

        iy = np.round(iy + CELL_SIZE, 8)
            
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
                           filename="metrics_summary.xlsx")
        

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
