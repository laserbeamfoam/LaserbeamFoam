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
  - Central configuration file for the workflow.
  - Acts as a single source of truth for other scripts

Assumptions:
  - Must be edited before running to reflect the target system setup

Authors
    Simon A. Rodriguez, University College Dublin (UCD). All rights reserved
    Petar Cosic, University College Dublin (UCD). All rights reserved
    Tom Flint, University of Manchester (UOM). All rights reserved
    Philip Cardiff, University College Dublin (UCD). All rights reserved

'''

LASER_DIAMETER = 10e-6 # In meters
OF_LOCATION = "$HOME/OpenFOAM/OpenFOAM-v2412/etc/bashrc"
CELL_SIZE = 2.5e-06
X_MIN_AND_MAX_DOMAIN = [0, 0.0001]
Y_COORD_BEGIN_TRACK = 40e-6
Y_COORD_END_TRACK = 110e-6
PLOT_GEOMETRY_VS_Y_LOCATION = True