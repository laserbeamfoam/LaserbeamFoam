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
  pvpython script that extracts the meltpool and save it as a csv file in the 
  current directory

Assumptions:
  - Called by higher-level script "characterise_meltpool.py"

Authors
    Simon A. Rodriguez, University College Dublin (UCD). All rights reserved
    Petar Cosic, University College Dublin (UCD). All rights reserved
    Tom Flint, University of Manchester (UOM). All rights reserved
    Philip Cardiff, University College Dublin (UCD). All rights reserved
    
'''


# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

from pathlib import Path
import numpy as np
from input_data import *

current_case = Path.cwd().name

# get active source.
test_case_1foam = GetActiveSource()

# Properties modified on test_case_1foam
test_case_1foam = OpenFOAMReader(FileName="./main.foam")
test_case_1foam.SkipZeroTime = 0
test_case_1foam.MeshRegions = ['back', 'bottomWall', 'front', 'internalMesh', 
                               'leftWall', 'rightWall', 'topWall']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1980, 907]

# show data in view
test_case_1foamDisplay = Show(test_case_1foam, renderView1)

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')
pLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
pLUT.InterpretValuesAsCategories = 0
pLUT.AnnotationsInitialized = 0
pLUT.ShowCategoricalColorsinDataRangeOnly = 0
pLUT.RescaleOnVisibilityChange = 0
pLUT.EnableOpacityMapping = 0
pLUT.RGBPoints = [-79534.890625, 0.231373, 0.298039, 0.752941, 490587.5546875, 0.865003, 0.865003, 0.865003, 1060710.0, 0.705882, 0.0156863, 0.14902]
pLUT.UseLogScale = 0
pLUT.ColorSpace = 'Diverging'
pLUT.UseBelowRangeColor = 0
pLUT.BelowRangeColor = [0.0, 0.0, 0.0]
pLUT.UseAboveRangeColor = 0
pLUT.AboveRangeColor = [0.5, 0.5, 0.5]
pLUT.NanColor = [1.0, 1.0, 0.0]
pLUT.NanOpacity = 1.0
pLUT.Discretize = 1
pLUT.NumberOfTableValues = 256
pLUT.ScalarRangeInitialized = 1.0
pLUT.HSVWrap = 0
pLUT.VectorComponent = 0
pLUT.VectorMode = 'Magnitude'
pLUT.AllowDuplicateScalars = 1
pLUT.Annotations = []
pLUT.ActiveAnnotatedValues = []
pLUT.IndexedColors = []
pLUT.IndexedOpacities = []

# trace defaults for the display properties.
test_case_1foamDisplay.Representation = 'Surface'
test_case_1foamDisplay.AmbientColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.ColorArrayName = ['POINTS', 'p']
test_case_1foamDisplay.DiffuseColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.LookupTable = pLUT
test_case_1foamDisplay.MapScalars = 1
test_case_1foamDisplay.MultiComponentsMapping = 0
test_case_1foamDisplay.InterpolateScalarsBeforeMapping = 1
test_case_1foamDisplay.Opacity = 1.0
test_case_1foamDisplay.PointSize = 2.0
test_case_1foamDisplay.LineWidth = 1.0
test_case_1foamDisplay.RenderLinesAsTubes = 0
test_case_1foamDisplay.RenderPointsAsSpheres = 0
test_case_1foamDisplay.Interpolation = 'Gouraud'
test_case_1foamDisplay.Specular = 0.0
test_case_1foamDisplay.SpecularColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.SpecularPower = 100.0
test_case_1foamDisplay.Luminosity = 0.0
test_case_1foamDisplay.Ambient = 0.0
test_case_1foamDisplay.Diffuse = 1.0
test_case_1foamDisplay.EdgeColor = [0.0, 0.0, 0.5]
test_case_1foamDisplay.BackfaceRepresentation = 'Follow Frontface'
test_case_1foamDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.BackfaceOpacity = 1.0
test_case_1foamDisplay.Position = [0.0, 0.0, 0.0]
test_case_1foamDisplay.Scale = [1.0, 1.0, 1.0]
test_case_1foamDisplay.Orientation = [0.0, 0.0, 0.0]
test_case_1foamDisplay.Origin = [0.0, 0.0, 0.0]
test_case_1foamDisplay.Pickable = 1
test_case_1foamDisplay.Texture = None
test_case_1foamDisplay.Triangulate = 0
test_case_1foamDisplay.UseShaderReplacements = 0
test_case_1foamDisplay.ShaderReplacements = ''
test_case_1foamDisplay.NonlinearSubdivisionLevel = 1
test_case_1foamDisplay.UseDataPartitions = 0
test_case_1foamDisplay.OSPRayUseScaleArray = 0
test_case_1foamDisplay.OSPRayScaleArray = 'p'
test_case_1foamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
test_case_1foamDisplay.OSPRayMaterial = 'None'
test_case_1foamDisplay.Orient = 0
test_case_1foamDisplay.OrientationMode = 'Direction'
test_case_1foamDisplay.SelectOrientationVectors = 'U'
test_case_1foamDisplay.Scaling = 0
test_case_1foamDisplay.ScaleMode = 'No Data Scaling Off'
test_case_1foamDisplay.ScaleFactor = 7.999999797903001e-05
test_case_1foamDisplay.SelectScaleArray = 'p'
test_case_1foamDisplay.GlyphType = 'Arrow'
test_case_1foamDisplay.UseGlyphTable = 0
test_case_1foamDisplay.GlyphTableIndexArray = 'p'
test_case_1foamDisplay.UseCompositeGlyphTable = 0
test_case_1foamDisplay.UseGlyphCullingAndLOD = 0
test_case_1foamDisplay.LODValues = []
test_case_1foamDisplay.ColorByLODIndex = 0
test_case_1foamDisplay.GaussianRadius = 3.999999898951501e-06
test_case_1foamDisplay.ShaderPreset = 'Sphere'
test_case_1foamDisplay.CustomTriangleScale = 3
test_case_1foamDisplay.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
test_case_1foamDisplay.Emissive = 0
test_case_1foamDisplay.ScaleByArray = 0
test_case_1foamDisplay.SetScaleArray = ['POINTS', 'p']
test_case_1foamDisplay.ScaleArrayComponent = ''
test_case_1foamDisplay.UseScaleFunction = 1
test_case_1foamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
test_case_1foamDisplay.OpacityByArray = 0
test_case_1foamDisplay.OpacityArray = ['POINTS', 'p']
test_case_1foamDisplay.OpacityArrayComponent = ''
test_case_1foamDisplay.OpacityTransferFunction = 'PiecewiseFunction'
test_case_1foamDisplay.DataAxesGrid = 'GridAxesRepresentation'
test_case_1foamDisplay.SelectionCellLabelBold = 0
test_case_1foamDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
test_case_1foamDisplay.SelectionCellLabelFontFamily = 'Arial'
test_case_1foamDisplay.SelectionCellLabelFontFile = ''
test_case_1foamDisplay.SelectionCellLabelFontSize = 18
test_case_1foamDisplay.SelectionCellLabelItalic = 0
test_case_1foamDisplay.SelectionCellLabelJustification = 'Left'
test_case_1foamDisplay.SelectionCellLabelOpacity = 1.0
test_case_1foamDisplay.SelectionCellLabelShadow = 0
test_case_1foamDisplay.SelectionPointLabelBold = 0
test_case_1foamDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
test_case_1foamDisplay.SelectionPointLabelFontFamily = 'Arial'
test_case_1foamDisplay.SelectionPointLabelFontFile = ''
test_case_1foamDisplay.SelectionPointLabelFontSize = 18
test_case_1foamDisplay.SelectionPointLabelItalic = 0
test_case_1foamDisplay.SelectionPointLabelJustification = 'Left'
test_case_1foamDisplay.SelectionPointLabelOpacity = 1.0
test_case_1foamDisplay.SelectionPointLabelShadow = 0
test_case_1foamDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
test_case_1foamDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
test_case_1foamDisplay.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
test_case_1foamDisplay.GlyphType.TipResolution = 6
test_case_1foamDisplay.GlyphType.TipRadius = 0.1
test_case_1foamDisplay.GlyphType.TipLength = 0.35
test_case_1foamDisplay.GlyphType.ShaftResolution = 6
test_case_1foamDisplay.GlyphType.ShaftRadius = 0.03
test_case_1foamDisplay.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
test_case_1foamDisplay.ScaleTransferFunction.Points = [-79534.890625, 0.0, 0.5, 0.0, 1060710.0, 1.0, 0.5, 0.0]
test_case_1foamDisplay.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
test_case_1foamDisplay.OpacityTransferFunction.Points = [-79534.890625, 0.0, 0.5, 0.0, 1060710.0, 1.0, 0.5, 0.0]
test_case_1foamDisplay.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
test_case_1foamDisplay.DataAxesGrid.XTitle = 'X Axis'
test_case_1foamDisplay.DataAxesGrid.YTitle = 'Y Axis'
test_case_1foamDisplay.DataAxesGrid.ZTitle = 'Z Axis'
test_case_1foamDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
test_case_1foamDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
test_case_1foamDisplay.DataAxesGrid.XTitleFontFile = ''
test_case_1foamDisplay.DataAxesGrid.XTitleBold = 0
test_case_1foamDisplay.DataAxesGrid.XTitleItalic = 0
test_case_1foamDisplay.DataAxesGrid.XTitleFontSize = 12
test_case_1foamDisplay.DataAxesGrid.XTitleShadow = 0
test_case_1foamDisplay.DataAxesGrid.XTitleOpacity = 1.0
test_case_1foamDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
test_case_1foamDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
test_case_1foamDisplay.DataAxesGrid.YTitleFontFile = ''
test_case_1foamDisplay.DataAxesGrid.YTitleBold = 0
test_case_1foamDisplay.DataAxesGrid.YTitleItalic = 0
test_case_1foamDisplay.DataAxesGrid.YTitleFontSize = 12
test_case_1foamDisplay.DataAxesGrid.YTitleShadow = 0
test_case_1foamDisplay.DataAxesGrid.YTitleOpacity = 1.0
test_case_1foamDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
test_case_1foamDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
test_case_1foamDisplay.DataAxesGrid.ZTitleFontFile = ''
test_case_1foamDisplay.DataAxesGrid.ZTitleBold = 0
test_case_1foamDisplay.DataAxesGrid.ZTitleItalic = 0
test_case_1foamDisplay.DataAxesGrid.ZTitleFontSize = 12
test_case_1foamDisplay.DataAxesGrid.ZTitleShadow = 0
test_case_1foamDisplay.DataAxesGrid.ZTitleOpacity = 1.0
test_case_1foamDisplay.DataAxesGrid.FacesToRender = 63
test_case_1foamDisplay.DataAxesGrid.CullBackface = 0
test_case_1foamDisplay.DataAxesGrid.CullFrontface = 1
test_case_1foamDisplay.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.DataAxesGrid.ShowGrid = 0
test_case_1foamDisplay.DataAxesGrid.ShowEdges = 1
test_case_1foamDisplay.DataAxesGrid.ShowTicks = 1
test_case_1foamDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
test_case_1foamDisplay.DataAxesGrid.AxesToLabel = 63
test_case_1foamDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
test_case_1foamDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
test_case_1foamDisplay.DataAxesGrid.XLabelFontFile = ''
test_case_1foamDisplay.DataAxesGrid.XLabelBold = 0
test_case_1foamDisplay.DataAxesGrid.XLabelItalic = 0
test_case_1foamDisplay.DataAxesGrid.XLabelFontSize = 12
test_case_1foamDisplay.DataAxesGrid.XLabelShadow = 0
test_case_1foamDisplay.DataAxesGrid.XLabelOpacity = 1.0
test_case_1foamDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
test_case_1foamDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
test_case_1foamDisplay.DataAxesGrid.YLabelFontFile = ''
test_case_1foamDisplay.DataAxesGrid.YLabelBold = 0
test_case_1foamDisplay.DataAxesGrid.YLabelItalic = 0
test_case_1foamDisplay.DataAxesGrid.YLabelFontSize = 12
test_case_1foamDisplay.DataAxesGrid.YLabelShadow = 0
test_case_1foamDisplay.DataAxesGrid.YLabelOpacity = 1.0
test_case_1foamDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
test_case_1foamDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
test_case_1foamDisplay.DataAxesGrid.ZLabelFontFile = ''
test_case_1foamDisplay.DataAxesGrid.ZLabelBold = 0
test_case_1foamDisplay.DataAxesGrid.ZLabelItalic = 0
test_case_1foamDisplay.DataAxesGrid.ZLabelFontSize = 12
test_case_1foamDisplay.DataAxesGrid.ZLabelShadow = 0
test_case_1foamDisplay.DataAxesGrid.ZLabelOpacity = 1.0
test_case_1foamDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
test_case_1foamDisplay.DataAxesGrid.XAxisPrecision = 2
test_case_1foamDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
test_case_1foamDisplay.DataAxesGrid.XAxisLabels = []
test_case_1foamDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
test_case_1foamDisplay.DataAxesGrid.YAxisPrecision = 2
test_case_1foamDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
test_case_1foamDisplay.DataAxesGrid.YAxisLabels = []
test_case_1foamDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
test_case_1foamDisplay.DataAxesGrid.ZAxisPrecision = 2
test_case_1foamDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
test_case_1foamDisplay.DataAxesGrid.ZAxisLabels = []
test_case_1foamDisplay.DataAxesGrid.UseCustomBounds = 0
test_case_1foamDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
test_case_1foamDisplay.PolarAxes.Visibility = 0
test_case_1foamDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
test_case_1foamDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
test_case_1foamDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
test_case_1foamDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
test_case_1foamDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
test_case_1foamDisplay.PolarAxes.EnableCustomRange = 0
test_case_1foamDisplay.PolarAxes.CustomRange = [0.0, 1.0]
test_case_1foamDisplay.PolarAxes.PolarAxisVisibility = 1
test_case_1foamDisplay.PolarAxes.RadialAxesVisibility = 1
test_case_1foamDisplay.PolarAxes.DrawRadialGridlines = 1
test_case_1foamDisplay.PolarAxes.PolarArcsVisibility = 1
test_case_1foamDisplay.PolarAxes.DrawPolarArcsGridlines = 1
test_case_1foamDisplay.PolarAxes.NumberOfRadialAxes = 0
test_case_1foamDisplay.PolarAxes.AutoSubdividePolarAxis = 1
test_case_1foamDisplay.PolarAxes.NumberOfPolarAxis = 0
test_case_1foamDisplay.PolarAxes.MinimumRadius = 0.0
test_case_1foamDisplay.PolarAxes.MinimumAngle = 0.0
test_case_1foamDisplay.PolarAxes.MaximumAngle = 90.0
test_case_1foamDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
test_case_1foamDisplay.PolarAxes.Ratio = 1.0
test_case_1foamDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
test_case_1foamDisplay.PolarAxes.PolarAxisTitleVisibility = 1
test_case_1foamDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
test_case_1foamDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
test_case_1foamDisplay.PolarAxes.PolarLabelVisibility = 1
test_case_1foamDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
test_case_1foamDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
test_case_1foamDisplay.PolarAxes.RadialLabelVisibility = 1
test_case_1foamDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
test_case_1foamDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
test_case_1foamDisplay.PolarAxes.RadialUnitsVisibility = 1
test_case_1foamDisplay.PolarAxes.ScreenSize = 10.0
test_case_1foamDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
test_case_1foamDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
test_case_1foamDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
test_case_1foamDisplay.PolarAxes.PolarAxisTitleFontFile = ''
test_case_1foamDisplay.PolarAxes.PolarAxisTitleBold = 0
test_case_1foamDisplay.PolarAxes.PolarAxisTitleItalic = 0
test_case_1foamDisplay.PolarAxes.PolarAxisTitleShadow = 0
test_case_1foamDisplay.PolarAxes.PolarAxisTitleFontSize = 12
test_case_1foamDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
test_case_1foamDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
test_case_1foamDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
test_case_1foamDisplay.PolarAxes.PolarAxisLabelFontFile = ''
test_case_1foamDisplay.PolarAxes.PolarAxisLabelBold = 0
test_case_1foamDisplay.PolarAxes.PolarAxisLabelItalic = 0
test_case_1foamDisplay.PolarAxes.PolarAxisLabelShadow = 0
test_case_1foamDisplay.PolarAxes.PolarAxisLabelFontSize = 12
test_case_1foamDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
test_case_1foamDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
test_case_1foamDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
test_case_1foamDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
test_case_1foamDisplay.PolarAxes.LastRadialAxisTextBold = 0
test_case_1foamDisplay.PolarAxes.LastRadialAxisTextItalic = 0
test_case_1foamDisplay.PolarAxes.LastRadialAxisTextShadow = 0
test_case_1foamDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
test_case_1foamDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
test_case_1foamDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
test_case_1foamDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
test_case_1foamDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
test_case_1foamDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
test_case_1foamDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
test_case_1foamDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
test_case_1foamDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
test_case_1foamDisplay.PolarAxes.EnableDistanceLOD = 1
test_case_1foamDisplay.PolarAxes.DistanceLODThreshold = 0.7
test_case_1foamDisplay.PolarAxes.EnableViewAngleLOD = 1
test_case_1foamDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
test_case_1foamDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
test_case_1foamDisplay.PolarAxes.PolarTicksVisibility = 1
test_case_1foamDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
test_case_1foamDisplay.PolarAxes.TickLocation = 'Both'
test_case_1foamDisplay.PolarAxes.AxisTickVisibility = 1
test_case_1foamDisplay.PolarAxes.AxisMinorTickVisibility = 0
test_case_1foamDisplay.PolarAxes.ArcTickVisibility = 1
test_case_1foamDisplay.PolarAxes.ArcMinorTickVisibility = 0
test_case_1foamDisplay.PolarAxes.DeltaAngleMajor = 10.0
test_case_1foamDisplay.PolarAxes.DeltaAngleMinor = 5.0
test_case_1foamDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
test_case_1foamDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
test_case_1foamDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
test_case_1foamDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
test_case_1foamDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
test_case_1foamDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
test_case_1foamDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
test_case_1foamDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
test_case_1foamDisplay.PolarAxes.ArcMajorTickSize = 0.0
test_case_1foamDisplay.PolarAxes.ArcTickRatioSize = 0.3
test_case_1foamDisplay.PolarAxes.ArcMajorTickThickness = 1.0
test_case_1foamDisplay.PolarAxes.ArcTickRatioThickness = 0.5
test_case_1foamDisplay.PolarAxes.Use2DMode = 0
test_case_1foamDisplay.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
test_case_1foamDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')
pPWF.Points = [-79534.890625, 0.0, 0.5, 0.0, 1060710.0, 1.0, 0.5, 0.0]
pPWF.AllowDuplicateScalars = 1
pPWF.UseLogScale = 0
pPWF.ScalarRangeInitialized = 1

# set scalar coloring
ColorBy(test_case_1foamDisplay, ('POINTS', 'T'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
test_case_1foamDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
test_case_1foamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'T'
tLUT = GetColorTransferFunction('T')
tLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
tLUT.InterpretValuesAsCategories = 0
tLUT.AnnotationsInitialized = 0
tLUT.ShowCategoricalColorsinDataRangeOnly = 0
tLUT.RescaleOnVisibilityChange = 0
tLUT.EnableOpacityMapping = 0
tLUT.RGBPoints = [298.0, 0.231373, 0.298039, 0.752941, 1856.8681640625, 0.865003, 0.865003, 0.865003, 3415.736328125, 0.705882, 0.0156863, 0.14902]
tLUT.UseLogScale = 0
tLUT.ColorSpace = 'Diverging'
tLUT.UseBelowRangeColor = 0
tLUT.BelowRangeColor = [0.0, 0.0, 0.0]
tLUT.UseAboveRangeColor = 0
tLUT.AboveRangeColor = [0.5, 0.5, 0.5]
tLUT.NanColor = [1.0, 1.0, 0.0]
tLUT.NanOpacity = 1.0
tLUT.Discretize = 1
tLUT.NumberOfTableValues = 256
tLUT.ScalarRangeInitialized = 1.0
tLUT.HSVWrap = 0
tLUT.VectorComponent = 0
tLUT.VectorMode = 'Magnitude'
tLUT.AllowDuplicateScalars = 1
tLUT.Annotations = []
tLUT.ActiveAnnotatedValues = []
tLUT.IndexedColors = []
tLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'T'
tPWF = GetOpacityTransferFunction('T')
tPWF.Points = [298.0, 0.0, 0.5, 0.0, 3415.736328125, 1.0, 0.5, 0.0]
tPWF.AllowDuplicateScalars = 1
tPWF.UseLogScale = 0
tPWF.ScalarRangeInitialized = 1

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

animationScene1.GoToLast()

# create a new 'Clip'
clip1 = Clip(Input=test_case_1foam)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'p']
clip1.Value = 474724.91796875
clip1.Invert = 1
clip1.Crinkleclip = 0
clip1.Exact = 0

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [9.999999747378752e-05, 0.00039999998989515007, 9.999999747378752e-05]
clip1.ClipType.Normal = [1.0, 0.0, 0.0]
clip1.ClipType.Offset = 0.0

# Properties modified on clip1
clip1.ClipType = 'Scalar'
# clip1.Scalars = ['POINTS', 'alpha.material']
clip1.Scalars = ['POINTS', 'alpha.metal']
clip1.Value = 0.5
clip1.Invert = 0

# show data in view
clip1Display = Show(clip1, renderView1)

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.AmbientColor = [1.0, 1.0, 1.0]
clip1Display.ColorArrayName = ['POINTS', 'p']
clip1Display.DiffuseColor = [1.0, 1.0, 1.0]
clip1Display.LookupTable = pLUT
clip1Display.MapScalars = 1
clip1Display.MultiComponentsMapping = 0
clip1Display.InterpolateScalarsBeforeMapping = 1
clip1Display.Opacity = 1.0
clip1Display.PointSize = 2.0
clip1Display.LineWidth = 1.0
clip1Display.RenderLinesAsTubes = 0
clip1Display.RenderPointsAsSpheres = 0
clip1Display.Interpolation = 'Gouraud'
clip1Display.Specular = 0.0
clip1Display.SpecularColor = [1.0, 1.0, 1.0]
clip1Display.SpecularPower = 100.0
clip1Display.Luminosity = 0.0
clip1Display.Ambient = 0.0
clip1Display.Diffuse = 1.0
clip1Display.EdgeColor = [0.0, 0.0, 0.5]
clip1Display.BackfaceRepresentation = 'Follow Frontface'
clip1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
clip1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
clip1Display.BackfaceOpacity = 1.0
clip1Display.Position = [0.0, 0.0, 0.0]
clip1Display.Scale = [1.0, 1.0, 1.0]
clip1Display.Orientation = [0.0, 0.0, 0.0]
clip1Display.Origin = [0.0, 0.0, 0.0]
clip1Display.Pickable = 1
clip1Display.Texture = None
clip1Display.Triangulate = 0
clip1Display.UseShaderReplacements = 0
clip1Display.ShaderReplacements = ''
clip1Display.NonlinearSubdivisionLevel = 1
clip1Display.UseDataPartitions = 0
clip1Display.OSPRayUseScaleArray = 0
clip1Display.OSPRayScaleArray = 'p'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.OSPRayMaterial = 'None'
clip1Display.Orient = 0
clip1Display.OrientationMode = 'Direction'
clip1Display.SelectOrientationVectors = 'U'
clip1Display.Scaling = 0
clip1Display.ScaleMode = 'No Data Scaling Off'
clip1Display.ScaleFactor = 7.999999797903001e-05
clip1Display.SelectScaleArray = 'p'
clip1Display.GlyphType = 'Arrow'
clip1Display.UseGlyphTable = 0
clip1Display.GlyphTableIndexArray = 'p'
clip1Display.UseCompositeGlyphTable = 0
clip1Display.UseGlyphCullingAndLOD = 0
clip1Display.LODValues = []
clip1Display.ColorByLODIndex = 0
clip1Display.GaussianRadius = 3.999999898951501e-06
clip1Display.ShaderPreset = 'Sphere'
clip1Display.CustomTriangleScale = 3
clip1Display.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
clip1Display.Emissive = 0
clip1Display.ScaleByArray = 0
clip1Display.SetScaleArray = ['POINTS', 'p']
clip1Display.ScaleArrayComponent = ''
clip1Display.UseScaleFunction = 1
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityByArray = 0
clip1Display.OpacityArray = ['POINTS', 'p']
clip1Display.OpacityArrayComponent = ''
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.SelectionCellLabelBold = 0
clip1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
clip1Display.SelectionCellLabelFontFamily = 'Arial'
clip1Display.SelectionCellLabelFontFile = ''
clip1Display.SelectionCellLabelFontSize = 18
clip1Display.SelectionCellLabelItalic = 0
clip1Display.SelectionCellLabelJustification = 'Left'
clip1Display.SelectionCellLabelOpacity = 1.0
clip1Display.SelectionCellLabelShadow = 0
clip1Display.SelectionPointLabelBold = 0
clip1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
clip1Display.SelectionPointLabelFontFamily = 'Arial'
clip1Display.SelectionPointLabelFontFile = ''
clip1Display.SelectionPointLabelFontSize = 18
clip1Display.SelectionPointLabelItalic = 0
clip1Display.SelectionPointLabelJustification = 'Left'
clip1Display.SelectionPointLabelOpacity = 1.0
clip1Display.SelectionPointLabelShadow = 0
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = pPWF
clip1Display.ScalarOpacityUnitDistance = 7.030666005480218e-06
clip1Display.ExtractedBlockIndex = 1
clip1Display.SelectMapper = 'Projected tetra'
clip1Display.SamplingDimensions = [128, 128, 128]
clip1Display.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
clip1Display.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
clip1Display.GlyphType.TipResolution = 6
clip1Display.GlyphType.TipRadius = 0.1
clip1Display.GlyphType.TipLength = 0.35
clip1Display.GlyphType.ShaftResolution = 6
clip1Display.GlyphType.ShaftRadius = 0.03
clip1Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [-111260.1640625, 0.0, 0.5, 0.0, 1060710.0, 1.0, 0.5, 0.0]
clip1Display.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [-111260.1640625, 0.0, 0.5, 0.0, 1060710.0, 1.0, 0.5, 0.0]
clip1Display.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip1Display.DataAxesGrid.XTitle = 'X Axis'
clip1Display.DataAxesGrid.YTitle = 'Y Axis'
clip1Display.DataAxesGrid.ZTitle = 'Z Axis'
clip1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.XTitleFontFamily = 'Arial'
clip1Display.DataAxesGrid.XTitleFontFile = ''
clip1Display.DataAxesGrid.XTitleBold = 0
clip1Display.DataAxesGrid.XTitleItalic = 0
clip1Display.DataAxesGrid.XTitleFontSize = 12
clip1Display.DataAxesGrid.XTitleShadow = 0
clip1Display.DataAxesGrid.XTitleOpacity = 1.0
clip1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.YTitleFontFamily = 'Arial'
clip1Display.DataAxesGrid.YTitleFontFile = ''
clip1Display.DataAxesGrid.YTitleBold = 0
clip1Display.DataAxesGrid.YTitleItalic = 0
clip1Display.DataAxesGrid.YTitleFontSize = 12
clip1Display.DataAxesGrid.YTitleShadow = 0
clip1Display.DataAxesGrid.YTitleOpacity = 1.0
clip1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.ZTitleFontFamily = 'Arial'
clip1Display.DataAxesGrid.ZTitleFontFile = ''
clip1Display.DataAxesGrid.ZTitleBold = 0
clip1Display.DataAxesGrid.ZTitleItalic = 0
clip1Display.DataAxesGrid.ZTitleFontSize = 12
clip1Display.DataAxesGrid.ZTitleShadow = 0
clip1Display.DataAxesGrid.ZTitleOpacity = 1.0
clip1Display.DataAxesGrid.FacesToRender = 63
clip1Display.DataAxesGrid.CullBackface = 0
clip1Display.DataAxesGrid.CullFrontface = 1
clip1Display.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
clip1Display.DataAxesGrid.ShowGrid = 0
clip1Display.DataAxesGrid.ShowEdges = 1
clip1Display.DataAxesGrid.ShowTicks = 1
clip1Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
clip1Display.DataAxesGrid.AxesToLabel = 63
clip1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.XLabelFontFamily = 'Arial'
clip1Display.DataAxesGrid.XLabelFontFile = ''
clip1Display.DataAxesGrid.XLabelBold = 0
clip1Display.DataAxesGrid.XLabelItalic = 0
clip1Display.DataAxesGrid.XLabelFontSize = 12
clip1Display.DataAxesGrid.XLabelShadow = 0
clip1Display.DataAxesGrid.XLabelOpacity = 1.0
clip1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.YLabelFontFamily = 'Arial'
clip1Display.DataAxesGrid.YLabelFontFile = ''
clip1Display.DataAxesGrid.YLabelBold = 0
clip1Display.DataAxesGrid.YLabelItalic = 0
clip1Display.DataAxesGrid.YLabelFontSize = 12
clip1Display.DataAxesGrid.YLabelShadow = 0
clip1Display.DataAxesGrid.YLabelOpacity = 1.0
clip1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
clip1Display.DataAxesGrid.ZLabelFontFamily = 'Arial'
clip1Display.DataAxesGrid.ZLabelFontFile = ''
clip1Display.DataAxesGrid.ZLabelBold = 0
clip1Display.DataAxesGrid.ZLabelItalic = 0
clip1Display.DataAxesGrid.ZLabelFontSize = 12
clip1Display.DataAxesGrid.ZLabelShadow = 0
clip1Display.DataAxesGrid.ZLabelOpacity = 1.0
clip1Display.DataAxesGrid.XAxisNotation = 'Mixed'
clip1Display.DataAxesGrid.XAxisPrecision = 2
clip1Display.DataAxesGrid.XAxisUseCustomLabels = 0
clip1Display.DataAxesGrid.XAxisLabels = []
clip1Display.DataAxesGrid.YAxisNotation = 'Mixed'
clip1Display.DataAxesGrid.YAxisPrecision = 2
clip1Display.DataAxesGrid.YAxisUseCustomLabels = 0
clip1Display.DataAxesGrid.YAxisLabels = []
clip1Display.DataAxesGrid.ZAxisNotation = 'Mixed'
clip1Display.DataAxesGrid.ZAxisPrecision = 2
clip1Display.DataAxesGrid.ZAxisUseCustomLabels = 0
clip1Display.DataAxesGrid.ZAxisLabels = []
clip1Display.DataAxesGrid.UseCustomBounds = 0
clip1Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip1Display.PolarAxes.Visibility = 0
clip1Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
clip1Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
clip1Display.PolarAxes.EnableCustomRange = 0
clip1Display.PolarAxes.CustomRange = [0.0, 1.0]
clip1Display.PolarAxes.PolarAxisVisibility = 1
clip1Display.PolarAxes.RadialAxesVisibility = 1
clip1Display.PolarAxes.DrawRadialGridlines = 1
clip1Display.PolarAxes.PolarArcsVisibility = 1
clip1Display.PolarAxes.DrawPolarArcsGridlines = 1
clip1Display.PolarAxes.NumberOfRadialAxes = 0
clip1Display.PolarAxes.AutoSubdividePolarAxis = 1
clip1Display.PolarAxes.NumberOfPolarAxis = 0
clip1Display.PolarAxes.MinimumRadius = 0.0
clip1Display.PolarAxes.MinimumAngle = 0.0
clip1Display.PolarAxes.MaximumAngle = 90.0
clip1Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
clip1Display.PolarAxes.Ratio = 1.0
clip1Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.PolarAxisTitleVisibility = 1
clip1Display.PolarAxes.PolarAxisTitle = 'Radial Distance'
clip1Display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
clip1Display.PolarAxes.PolarLabelVisibility = 1
clip1Display.PolarAxes.PolarLabelFormat = '%-#6.3g'
clip1Display.PolarAxes.PolarLabelExponentLocation = 'Labels'
clip1Display.PolarAxes.RadialLabelVisibility = 1
clip1Display.PolarAxes.RadialLabelFormat = '%-#3.1f'
clip1Display.PolarAxes.RadialLabelLocation = 'Bottom'
clip1Display.PolarAxes.RadialUnitsVisibility = 1
clip1Display.PolarAxes.ScreenSize = 10.0
clip1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.PolarAxisTitleOpacity = 1.0
clip1Display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
clip1Display.PolarAxes.PolarAxisTitleFontFile = ''
clip1Display.PolarAxes.PolarAxisTitleBold = 0
clip1Display.PolarAxes.PolarAxisTitleItalic = 0
clip1Display.PolarAxes.PolarAxisTitleShadow = 0
clip1Display.PolarAxes.PolarAxisTitleFontSize = 12
clip1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.PolarAxisLabelOpacity = 1.0
clip1Display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
clip1Display.PolarAxes.PolarAxisLabelFontFile = ''
clip1Display.PolarAxes.PolarAxisLabelBold = 0
clip1Display.PolarAxes.PolarAxisLabelItalic = 0
clip1Display.PolarAxes.PolarAxisLabelShadow = 0
clip1Display.PolarAxes.PolarAxisLabelFontSize = 12
clip1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
clip1Display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
clip1Display.PolarAxes.LastRadialAxisTextFontFile = ''
clip1Display.PolarAxes.LastRadialAxisTextBold = 0
clip1Display.PolarAxes.LastRadialAxisTextItalic = 0
clip1Display.PolarAxes.LastRadialAxisTextShadow = 0
clip1Display.PolarAxes.LastRadialAxisTextFontSize = 12
clip1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
clip1Display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
clip1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
clip1Display.PolarAxes.SecondaryRadialAxesTextBold = 0
clip1Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
clip1Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
clip1Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
clip1Display.PolarAxes.EnableDistanceLOD = 1
clip1Display.PolarAxes.DistanceLODThreshold = 0.7
clip1Display.PolarAxes.EnableViewAngleLOD = 1
clip1Display.PolarAxes.ViewAngleLODThreshold = 0.7
clip1Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
clip1Display.PolarAxes.PolarTicksVisibility = 1
clip1Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
clip1Display.PolarAxes.TickLocation = 'Both'
clip1Display.PolarAxes.AxisTickVisibility = 1
clip1Display.PolarAxes.AxisMinorTickVisibility = 0
clip1Display.PolarAxes.ArcTickVisibility = 1
clip1Display.PolarAxes.ArcMinorTickVisibility = 0
clip1Display.PolarAxes.DeltaAngleMajor = 10.0
clip1Display.PolarAxes.DeltaAngleMinor = 5.0
clip1Display.PolarAxes.PolarAxisMajorTickSize = 0.0
clip1Display.PolarAxes.PolarAxisTickRatioSize = 0.3
clip1Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
clip1Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
clip1Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
clip1Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
clip1Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
clip1Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
clip1Display.PolarAxes.ArcMajorTickSize = 0.0
clip1Display.PolarAxes.ArcTickRatioSize = 0.3
clip1Display.PolarAxes.ArcMajorTickThickness = 1.0
clip1Display.PolarAxes.ArcTickRatioThickness = 0.5
clip1Display.PolarAxes.Use2DMode = 0
clip1Display.PolarAxes.UseLogAxis = 0

# hide data in view
Hide(test_case_1foam, renderView1)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'T'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# current camera placement for renderView1
renderView1.CameraPosition = [9.999999747378752e-05, 0.00039999998989515007, 0.0017392304406046689]
renderView1.CameraFocalPoint = [9.999999747378752e-05, 0.00039999998989515007, 9.999999747378752e-05]
renderView1.CameraParallelScale = 0.00042426405799411666

renderView1.ResetCamera()
renderView1.StillRender()


# get layout
layout1 = GetLayout()

# save screenshot
# SaveScreenshot('./Clip1.png', layout1, ImageResolution=[1980, 907],
#     FontScaling='Scale fonts proportionally',
#     OverrideColorPalette='',
#     StereoMode='No change',
#     TransparentBackground=0, 
#     # PNG options
#     CompressionLevel='5')

# create a new 'Clip'
clip2 = Clip(Input=clip1)
clip2.ClipType = 'Plane'
clip2.Scalars = ['POINTS', 'p']
clip2.Value = 474724.91796875
clip2.Invert = 1
clip2.Crinkleclip = 0
clip2.Exact = 0

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [9.999999747378752e-05, 0.00039999998989515007, 7.499999628635123e-05]
clip2.ClipType.Normal = [1.0, 0.0, 0.0]
clip2.ClipType.Offset = 0.0

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip2.ClipType)

# Properties modified on clip2
clip2.Invert = 0

# Properties modified on clip2.ClipType
clip2.ClipType.Origin = [9.999999747378752e-05, Y_COORD_BEGIN_TRACK + LASER_DIAMETER/2, 7.499999628635123e-05]
clip2.ClipType.Normal = [0.0, 1.0, 0.0]

# show data in view
clip2Display = Show(clip2, renderView1)

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.AmbientColor = [1.0, 1.0, 1.0]
clip2Display.ColorArrayName = ['POINTS', 'p']
clip2Display.DiffuseColor = [1.0, 1.0, 1.0]
clip2Display.LookupTable = pLUT
clip2Display.MapScalars = 1
clip2Display.MultiComponentsMapping = 0
clip2Display.InterpolateScalarsBeforeMapping = 1
clip2Display.Opacity = 1.0
clip2Display.PointSize = 2.0
clip2Display.LineWidth = 1.0
clip2Display.RenderLinesAsTubes = 0
clip2Display.RenderPointsAsSpheres = 0
clip2Display.Interpolation = 'Gouraud'
clip2Display.Specular = 0.0
clip2Display.SpecularColor = [1.0, 1.0, 1.0]
clip2Display.SpecularPower = 100.0
clip2Display.Luminosity = 0.0
clip2Display.Ambient = 0.0
clip2Display.Diffuse = 1.0
clip2Display.EdgeColor = [0.0, 0.0, 0.5]
clip2Display.BackfaceRepresentation = 'Follow Frontface'
clip2Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
clip2Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
clip2Display.BackfaceOpacity = 1.0
clip2Display.Position = [0.0, 0.0, 0.0]
clip2Display.Scale = [1.0, 1.0, 1.0]
clip2Display.Orientation = [0.0, 0.0, 0.0]
clip2Display.Origin = [0.0, 0.0, 0.0]
clip2Display.Pickable = 1
clip2Display.Texture = None
clip2Display.Triangulate = 0
clip2Display.UseShaderReplacements = 0
clip2Display.ShaderReplacements = ''
clip2Display.NonlinearSubdivisionLevel = 1
clip2Display.UseDataPartitions = 0
clip2Display.OSPRayUseScaleArray = 0
clip2Display.OSPRayScaleArray = 'p'
clip2Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip2Display.OSPRayMaterial = 'None'
clip2Display.Orient = 0
clip2Display.OrientationMode = 'Direction'
clip2Display.SelectOrientationVectors = 'U'
clip2Display.Scaling = 0
clip2Display.ScaleMode = 'No Data Scaling Off'
clip2Display.ScaleFactor = 6.5999997605104e-05
clip2Display.SelectScaleArray = 'p'
clip2Display.GlyphType = 'Arrow'
clip2Display.UseGlyphTable = 0
clip2Display.GlyphTableIndexArray = 'p'
clip2Display.UseCompositeGlyphTable = 0
clip2Display.UseGlyphCullingAndLOD = 0
clip2Display.LODValues = []
clip2Display.ColorByLODIndex = 0
clip2Display.GaussianRadius = 3.2999998802552e-06
clip2Display.ShaderPreset = 'Sphere'
clip2Display.CustomTriangleScale = 3
clip2Display.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
clip2Display.Emissive = 0
clip2Display.ScaleByArray = 0
clip2Display.SetScaleArray = ['POINTS', 'p']
clip2Display.ScaleArrayComponent = ''
clip2Display.UseScaleFunction = 1
clip2Display.ScaleTransferFunction = 'PiecewiseFunction'
clip2Display.OpacityByArray = 0
clip2Display.OpacityArray = ['POINTS', 'p']
clip2Display.OpacityArrayComponent = ''
clip2Display.OpacityTransferFunction = 'PiecewiseFunction'
clip2Display.DataAxesGrid = 'GridAxesRepresentation'
clip2Display.SelectionCellLabelBold = 0
clip2Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
clip2Display.SelectionCellLabelFontFamily = 'Arial'
clip2Display.SelectionCellLabelFontFile = ''
clip2Display.SelectionCellLabelFontSize = 18
clip2Display.SelectionCellLabelItalic = 0
clip2Display.SelectionCellLabelJustification = 'Left'
clip2Display.SelectionCellLabelOpacity = 1.0
clip2Display.SelectionCellLabelShadow = 0
clip2Display.SelectionPointLabelBold = 0
clip2Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
clip2Display.SelectionPointLabelFontFamily = 'Arial'
clip2Display.SelectionPointLabelFontFile = ''
clip2Display.SelectionPointLabelFontSize = 18
clip2Display.SelectionPointLabelItalic = 0
clip2Display.SelectionPointLabelJustification = 'Left'
clip2Display.SelectionPointLabelOpacity = 1.0
clip2Display.SelectionPointLabelShadow = 0
clip2Display.PolarAxes = 'PolarAxesRepresentation'
clip2Display.ScalarOpacityFunction = pPWF
clip2Display.ScalarOpacityUnitDistance = 6.33561131050142e-06
clip2Display.ExtractedBlockIndex = 1
clip2Display.SelectMapper = 'Projected tetra'
clip2Display.SamplingDimensions = [128, 128, 128]
clip2Display.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
clip2Display.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
clip2Display.GlyphType.TipResolution = 6
clip2Display.GlyphType.TipRadius = 0.1
clip2Display.GlyphType.TipLength = 0.35
clip2Display.GlyphType.ShaftResolution = 6
clip2Display.GlyphType.ShaftRadius = 0.03
clip2Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip2Display.ScaleTransferFunction.Points = [-100115.2890625, 0.0, 0.5, 0.0, 1060710.0, 1.0, 0.5, 0.0]
clip2Display.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip2Display.OpacityTransferFunction.Points = [-100115.2890625, 0.0, 0.5, 0.0, 1060710.0, 1.0, 0.5, 0.0]
clip2Display.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip2Display.DataAxesGrid.XTitle = 'X Axis'
clip2Display.DataAxesGrid.YTitle = 'Y Axis'
clip2Display.DataAxesGrid.ZTitle = 'Z Axis'
clip2Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
clip2Display.DataAxesGrid.XTitleFontFamily = 'Arial'
clip2Display.DataAxesGrid.XTitleFontFile = ''
clip2Display.DataAxesGrid.XTitleBold = 0
clip2Display.DataAxesGrid.XTitleItalic = 0
clip2Display.DataAxesGrid.XTitleFontSize = 12
clip2Display.DataAxesGrid.XTitleShadow = 0
clip2Display.DataAxesGrid.XTitleOpacity = 1.0
clip2Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
clip2Display.DataAxesGrid.YTitleFontFamily = 'Arial'
clip2Display.DataAxesGrid.YTitleFontFile = ''
clip2Display.DataAxesGrid.YTitleBold = 0
clip2Display.DataAxesGrid.YTitleItalic = 0
clip2Display.DataAxesGrid.YTitleFontSize = 12
clip2Display.DataAxesGrid.YTitleShadow = 0
clip2Display.DataAxesGrid.YTitleOpacity = 1.0
clip2Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
clip2Display.DataAxesGrid.ZTitleFontFamily = 'Arial'
clip2Display.DataAxesGrid.ZTitleFontFile = ''
clip2Display.DataAxesGrid.ZTitleBold = 0
clip2Display.DataAxesGrid.ZTitleItalic = 0
clip2Display.DataAxesGrid.ZTitleFontSize = 12
clip2Display.DataAxesGrid.ZTitleShadow = 0
clip2Display.DataAxesGrid.ZTitleOpacity = 1.0
clip2Display.DataAxesGrid.FacesToRender = 63
clip2Display.DataAxesGrid.CullBackface = 0
clip2Display.DataAxesGrid.CullFrontface = 1
clip2Display.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
clip2Display.DataAxesGrid.ShowGrid = 0
clip2Display.DataAxesGrid.ShowEdges = 1
clip2Display.DataAxesGrid.ShowTicks = 1
clip2Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
clip2Display.DataAxesGrid.AxesToLabel = 63
clip2Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
clip2Display.DataAxesGrid.XLabelFontFamily = 'Arial'
clip2Display.DataAxesGrid.XLabelFontFile = ''
clip2Display.DataAxesGrid.XLabelBold = 0
clip2Display.DataAxesGrid.XLabelItalic = 0
clip2Display.DataAxesGrid.XLabelFontSize = 12
clip2Display.DataAxesGrid.XLabelShadow = 0
clip2Display.DataAxesGrid.XLabelOpacity = 1.0
clip2Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
clip2Display.DataAxesGrid.YLabelFontFamily = 'Arial'
clip2Display.DataAxesGrid.YLabelFontFile = ''
clip2Display.DataAxesGrid.YLabelBold = 0
clip2Display.DataAxesGrid.YLabelItalic = 0
clip2Display.DataAxesGrid.YLabelFontSize = 12
clip2Display.DataAxesGrid.YLabelShadow = 0
clip2Display.DataAxesGrid.YLabelOpacity = 1.0
clip2Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
clip2Display.DataAxesGrid.ZLabelFontFamily = 'Arial'
clip2Display.DataAxesGrid.ZLabelFontFile = ''
clip2Display.DataAxesGrid.ZLabelBold = 0
clip2Display.DataAxesGrid.ZLabelItalic = 0
clip2Display.DataAxesGrid.ZLabelFontSize = 12
clip2Display.DataAxesGrid.ZLabelShadow = 0
clip2Display.DataAxesGrid.ZLabelOpacity = 1.0
clip2Display.DataAxesGrid.XAxisNotation = 'Mixed'
clip2Display.DataAxesGrid.XAxisPrecision = 2
clip2Display.DataAxesGrid.XAxisUseCustomLabels = 0
clip2Display.DataAxesGrid.XAxisLabels = []
clip2Display.DataAxesGrid.YAxisNotation = 'Mixed'
clip2Display.DataAxesGrid.YAxisPrecision = 2
clip2Display.DataAxesGrid.YAxisUseCustomLabels = 0
clip2Display.DataAxesGrid.YAxisLabels = []
clip2Display.DataAxesGrid.ZAxisNotation = 'Mixed'
clip2Display.DataAxesGrid.ZAxisPrecision = 2
clip2Display.DataAxesGrid.ZAxisUseCustomLabels = 0
clip2Display.DataAxesGrid.ZAxisLabels = []
clip2Display.DataAxesGrid.UseCustomBounds = 0
clip2Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip2Display.PolarAxes.Visibility = 0
clip2Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
clip2Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
clip2Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
clip2Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
clip2Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
clip2Display.PolarAxes.EnableCustomRange = 0
clip2Display.PolarAxes.CustomRange = [0.0, 1.0]
clip2Display.PolarAxes.PolarAxisVisibility = 1
clip2Display.PolarAxes.RadialAxesVisibility = 1
clip2Display.PolarAxes.DrawRadialGridlines = 1
clip2Display.PolarAxes.PolarArcsVisibility = 1
clip2Display.PolarAxes.DrawPolarArcsGridlines = 1
clip2Display.PolarAxes.NumberOfRadialAxes = 0
clip2Display.PolarAxes.AutoSubdividePolarAxis = 1
clip2Display.PolarAxes.NumberOfPolarAxis = 0
clip2Display.PolarAxes.MinimumRadius = 0.0
clip2Display.PolarAxes.MinimumAngle = 0.0
clip2Display.PolarAxes.MaximumAngle = 90.0
clip2Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
clip2Display.PolarAxes.Ratio = 1.0
clip2Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
clip2Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
clip2Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
clip2Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
clip2Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
clip2Display.PolarAxes.PolarAxisTitleVisibility = 1
clip2Display.PolarAxes.PolarAxisTitle = 'Radial Distance'
clip2Display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
clip2Display.PolarAxes.PolarLabelVisibility = 1
clip2Display.PolarAxes.PolarLabelFormat = '%-#6.3g'
clip2Display.PolarAxes.PolarLabelExponentLocation = 'Labels'
clip2Display.PolarAxes.RadialLabelVisibility = 1
clip2Display.PolarAxes.RadialLabelFormat = '%-#3.1f'
clip2Display.PolarAxes.RadialLabelLocation = 'Bottom'
clip2Display.PolarAxes.RadialUnitsVisibility = 1
clip2Display.PolarAxes.ScreenSize = 10.0
clip2Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
clip2Display.PolarAxes.PolarAxisTitleOpacity = 1.0
clip2Display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
clip2Display.PolarAxes.PolarAxisTitleFontFile = ''
clip2Display.PolarAxes.PolarAxisTitleBold = 0
clip2Display.PolarAxes.PolarAxisTitleItalic = 0
clip2Display.PolarAxes.PolarAxisTitleShadow = 0
clip2Display.PolarAxes.PolarAxisTitleFontSize = 12
clip2Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
clip2Display.PolarAxes.PolarAxisLabelOpacity = 1.0
clip2Display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
clip2Display.PolarAxes.PolarAxisLabelFontFile = ''
clip2Display.PolarAxes.PolarAxisLabelBold = 0
clip2Display.PolarAxes.PolarAxisLabelItalic = 0
clip2Display.PolarAxes.PolarAxisLabelShadow = 0
clip2Display.PolarAxes.PolarAxisLabelFontSize = 12
clip2Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
clip2Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
clip2Display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
clip2Display.PolarAxes.LastRadialAxisTextFontFile = ''
clip2Display.PolarAxes.LastRadialAxisTextBold = 0
clip2Display.PolarAxes.LastRadialAxisTextItalic = 0
clip2Display.PolarAxes.LastRadialAxisTextShadow = 0
clip2Display.PolarAxes.LastRadialAxisTextFontSize = 12
clip2Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
clip2Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
clip2Display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
clip2Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
clip2Display.PolarAxes.SecondaryRadialAxesTextBold = 0
clip2Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
clip2Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
clip2Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
clip2Display.PolarAxes.EnableDistanceLOD = 1
clip2Display.PolarAxes.DistanceLODThreshold = 0.7
clip2Display.PolarAxes.EnableViewAngleLOD = 1
clip2Display.PolarAxes.ViewAngleLODThreshold = 0.7
clip2Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
clip2Display.PolarAxes.PolarTicksVisibility = 1
clip2Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
clip2Display.PolarAxes.TickLocation = 'Both'
clip2Display.PolarAxes.AxisTickVisibility = 1
clip2Display.PolarAxes.AxisMinorTickVisibility = 0
clip2Display.PolarAxes.ArcTickVisibility = 1
clip2Display.PolarAxes.ArcMinorTickVisibility = 0
clip2Display.PolarAxes.DeltaAngleMajor = 10.0
clip2Display.PolarAxes.DeltaAngleMinor = 5.0
clip2Display.PolarAxes.PolarAxisMajorTickSize = 0.0
clip2Display.PolarAxes.PolarAxisTickRatioSize = 0.3
clip2Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
clip2Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
clip2Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
clip2Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
clip2Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
clip2Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
clip2Display.PolarAxes.ArcMajorTickSize = 0.0
clip2Display.PolarAxes.ArcTickRatioSize = 0.3
clip2Display.PolarAxes.ArcMajorTickThickness = 1.0
clip2Display.PolarAxes.ArcTickRatioThickness = 0.5
clip2Display.PolarAxes.Use2DMode = 0
clip2Display.PolarAxes.UseLogAxis = 0

# hide data in view
Hide(clip1, renderView1)

# show color bar/color legend
clip2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip2Display, ('POINTS', 'T'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip2Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip2Display.SetScalarBarVisibility(renderView1, True)

# current camera placement for renderView1
renderView1.CameraPosition = [9.999999747378752e-05, 0.00039999998989515007, 
                              0.0017392304406046689]
renderView1.CameraFocalPoint = [9.999999747378752e-05, 0.00039999998989515007, 
                                9.999999747378752e-05]
renderView1.CameraParallelScale = 0.00042426405799411666

renderView1.ResetCamera()
renderView1.StillRender()

# save screenshot
# SaveScreenshot('./Clip2.png', layout1, ImageResolution=[1980, 907],
#     FontScaling='Scale fonts proportionally',
#     OverrideColorPalette='',
#     StereoMode='No change',
#     TransparentBackground=0, 
#     # PNG options
#     CompressionLevel='5')

# create a new 'Clip'
clip3 = Clip(Input=clip2)
clip3.ClipType = 'Plane'
clip3.Scalars = ['POINTS', 'p']
clip3.Value = 480297.35546875
clip3.Invert = 1
clip3.Crinkleclip = 0
clip3.Exact = 0

# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Origin = [9.999999747378752e-05, 0.00046999999176478013, 
                         7.375000132014975e-05]
clip3.ClipType.Normal = [1.0, 0.0, 0.0]
clip3.ClipType.Offset = 0.0

# Properties modified on clip3.ClipType
clip3.ClipType.Origin = [9.999999747378752e-05, Y_COORD_END_TRACK - 
                                       LASER_DIAMETER/2, 7.375000132014975e-05]
clip3.ClipType.Normal = [0.0, 1.0, 0.0]

# show data in view
clip3Display = Show(clip3, renderView1)

# trace defaults for the display properties.
clip3Display.Representation = 'Surface'
clip3Display.AmbientColor = [1.0, 1.0, 1.0]
clip3Display.ColorArrayName = ['POINTS', 'p']
clip3Display.DiffuseColor = [1.0, 1.0, 1.0]
clip3Display.LookupTable = pLUT
clip3Display.MapScalars = 1
clip3Display.MultiComponentsMapping = 0
clip3Display.InterpolateScalarsBeforeMapping = 1
clip3Display.Opacity = 1.0
clip3Display.PointSize = 2.0
clip3Display.LineWidth = 1.0
clip3Display.RenderLinesAsTubes = 0
clip3Display.RenderPointsAsSpheres = 0
clip3Display.Interpolation = 'Gouraud'
clip3Display.Specular = 0.0
clip3Display.SpecularColor = [1.0, 1.0, 1.0]
clip3Display.SpecularPower = 100.0
clip3Display.Luminosity = 0.0
clip3Display.Ambient = 0.0
clip3Display.Diffuse = 1.0
clip3Display.EdgeColor = [0.0, 0.0, 0.5]
clip3Display.BackfaceRepresentation = 'Follow Frontface'
clip3Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
clip3Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
clip3Display.BackfaceOpacity = 1.0
clip3Display.Position = [0.0, 0.0, 0.0]
clip3Display.Scale = [1.0, 1.0, 1.0]
clip3Display.Orientation = [0.0, 0.0, 0.0]
clip3Display.Origin = [0.0, 0.0, 0.0]
clip3Display.Pickable = 1
clip3Display.Texture = None
clip3Display.Triangulate = 0
clip3Display.UseShaderReplacements = 0
clip3Display.ShaderReplacements = ''
clip3Display.NonlinearSubdivisionLevel = 1
clip3Display.UseDataPartitions = 0
clip3Display.OSPRayUseScaleArray = 0
clip3Display.OSPRayScaleArray = 'p'
clip3Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip3Display.OSPRayMaterial = 'None'
clip3Display.Orient = 0
clip3Display.OrientationMode = 'Direction'
clip3Display.SelectOrientationVectors = 'U'
clip3Display.Scaling = 0
clip3Display.ScaleMode = 'No Data Scaling Off'
clip3Display.ScaleFactor = 5.200000159675256e-05
clip3Display.SelectScaleArray = 'p'
clip3Display.GlyphType = 'Arrow'
clip3Display.UseGlyphTable = 0
clip3Display.GlyphTableIndexArray = 'p'
clip3Display.UseCompositeGlyphTable = 0
clip3Display.UseGlyphCullingAndLOD = 0
clip3Display.LODValues = []
clip3Display.ColorByLODIndex = 0
clip3Display.GaussianRadius = 2.6000000798376276e-06
clip3Display.ShaderPreset = 'Sphere'
clip3Display.CustomTriangleScale = 3
clip3Display.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
clip3Display.Emissive = 0
clip3Display.ScaleByArray = 0
clip3Display.SetScaleArray = ['POINTS', 'p']
clip3Display.ScaleArrayComponent = ''
clip3Display.UseScaleFunction = 1
clip3Display.ScaleTransferFunction = 'PiecewiseFunction'
clip3Display.OpacityByArray = 0
clip3Display.OpacityArray = ['POINTS', 'p']
clip3Display.OpacityArrayComponent = ''
clip3Display.OpacityTransferFunction = 'PiecewiseFunction'
clip3Display.DataAxesGrid = 'GridAxesRepresentation'
clip3Display.SelectionCellLabelBold = 0
clip3Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
clip3Display.SelectionCellLabelFontFamily = 'Arial'
clip3Display.SelectionCellLabelFontFile = ''
clip3Display.SelectionCellLabelFontSize = 18
clip3Display.SelectionCellLabelItalic = 0
clip3Display.SelectionCellLabelJustification = 'Left'
clip3Display.SelectionCellLabelOpacity = 1.0
clip3Display.SelectionCellLabelShadow = 0
clip3Display.SelectionPointLabelBold = 0
clip3Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
clip3Display.SelectionPointLabelFontFamily = 'Arial'
clip3Display.SelectionPointLabelFontFile = ''
clip3Display.SelectionPointLabelFontSize = 18
clip3Display.SelectionPointLabelItalic = 0
clip3Display.SelectionPointLabelJustification = 'Left'
clip3Display.SelectionPointLabelOpacity = 1.0
clip3Display.SelectionPointLabelShadow = 0
clip3Display.PolarAxes = 'PolarAxesRepresentation'
clip3Display.ScalarOpacityFunction = pPWF
clip3Display.ScalarOpacityUnitDistance = 5.623239319539663e-06
clip3Display.ExtractedBlockIndex = 1
clip3Display.SelectMapper = 'Projected tetra'
clip3Display.SamplingDimensions = [128, 128, 128]
clip3Display.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip3Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
clip3Display.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
clip3Display.GlyphType.TipResolution = 6
clip3Display.GlyphType.TipRadius = 0.1
clip3Display.GlyphType.TipLength = 0.35
clip3Display.GlyphType.ShaftResolution = 6
clip3Display.GlyphType.ShaftRadius = 0.03
clip3Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip3Display.ScaleTransferFunction.Points = [-100115.2890625, 0.0, 0.5, 0.0, 1060710.0, 1.0, 0.5, 0.0]
clip3Display.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip3Display.OpacityTransferFunction.Points = [-100115.2890625, 0.0, 0.5, 0.0, 1060710.0, 1.0, 0.5, 0.0]
clip3Display.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip3Display.DataAxesGrid.XTitle = 'X Axis'
clip3Display.DataAxesGrid.YTitle = 'Y Axis'
clip3Display.DataAxesGrid.ZTitle = 'Z Axis'
clip3Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
clip3Display.DataAxesGrid.XTitleFontFamily = 'Arial'
clip3Display.DataAxesGrid.XTitleFontFile = ''
clip3Display.DataAxesGrid.XTitleBold = 0
clip3Display.DataAxesGrid.XTitleItalic = 0
clip3Display.DataAxesGrid.XTitleFontSize = 12
clip3Display.DataAxesGrid.XTitleShadow = 0
clip3Display.DataAxesGrid.XTitleOpacity = 1.0
clip3Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
clip3Display.DataAxesGrid.YTitleFontFamily = 'Arial'
clip3Display.DataAxesGrid.YTitleFontFile = ''
clip3Display.DataAxesGrid.YTitleBold = 0
clip3Display.DataAxesGrid.YTitleItalic = 0
clip3Display.DataAxesGrid.YTitleFontSize = 12
clip3Display.DataAxesGrid.YTitleShadow = 0
clip3Display.DataAxesGrid.YTitleOpacity = 1.0
clip3Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
clip3Display.DataAxesGrid.ZTitleFontFamily = 'Arial'
clip3Display.DataAxesGrid.ZTitleFontFile = ''
clip3Display.DataAxesGrid.ZTitleBold = 0
clip3Display.DataAxesGrid.ZTitleItalic = 0
clip3Display.DataAxesGrid.ZTitleFontSize = 12
clip3Display.DataAxesGrid.ZTitleShadow = 0
clip3Display.DataAxesGrid.ZTitleOpacity = 1.0
clip3Display.DataAxesGrid.FacesToRender = 63
clip3Display.DataAxesGrid.CullBackface = 0
clip3Display.DataAxesGrid.CullFrontface = 1
clip3Display.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
clip3Display.DataAxesGrid.ShowGrid = 0
clip3Display.DataAxesGrid.ShowEdges = 1
clip3Display.DataAxesGrid.ShowTicks = 1
clip3Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
clip3Display.DataAxesGrid.AxesToLabel = 63
clip3Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
clip3Display.DataAxesGrid.XLabelFontFamily = 'Arial'
clip3Display.DataAxesGrid.XLabelFontFile = ''
clip3Display.DataAxesGrid.XLabelBold = 0
clip3Display.DataAxesGrid.XLabelItalic = 0
clip3Display.DataAxesGrid.XLabelFontSize = 12
clip3Display.DataAxesGrid.XLabelShadow = 0
clip3Display.DataAxesGrid.XLabelOpacity = 1.0
clip3Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
clip3Display.DataAxesGrid.YLabelFontFamily = 'Arial'
clip3Display.DataAxesGrid.YLabelFontFile = ''
clip3Display.DataAxesGrid.YLabelBold = 0
clip3Display.DataAxesGrid.YLabelItalic = 0
clip3Display.DataAxesGrid.YLabelFontSize = 12
clip3Display.DataAxesGrid.YLabelShadow = 0
clip3Display.DataAxesGrid.YLabelOpacity = 1.0
clip3Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
clip3Display.DataAxesGrid.ZLabelFontFamily = 'Arial'
clip3Display.DataAxesGrid.ZLabelFontFile = ''
clip3Display.DataAxesGrid.ZLabelBold = 0
clip3Display.DataAxesGrid.ZLabelItalic = 0
clip3Display.DataAxesGrid.ZLabelFontSize = 12
clip3Display.DataAxesGrid.ZLabelShadow = 0
clip3Display.DataAxesGrid.ZLabelOpacity = 1.0
clip3Display.DataAxesGrid.XAxisNotation = 'Mixed'
clip3Display.DataAxesGrid.XAxisPrecision = 2
clip3Display.DataAxesGrid.XAxisUseCustomLabels = 0
clip3Display.DataAxesGrid.XAxisLabels = []
clip3Display.DataAxesGrid.YAxisNotation = 'Mixed'
clip3Display.DataAxesGrid.YAxisPrecision = 2
clip3Display.DataAxesGrid.YAxisUseCustomLabels = 0
clip3Display.DataAxesGrid.YAxisLabels = []
clip3Display.DataAxesGrid.ZAxisNotation = 'Mixed'
clip3Display.DataAxesGrid.ZAxisPrecision = 2
clip3Display.DataAxesGrid.ZAxisUseCustomLabels = 0
clip3Display.DataAxesGrid.ZAxisLabels = []
clip3Display.DataAxesGrid.UseCustomBounds = 0
clip3Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip3Display.PolarAxes.Visibility = 0
clip3Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
clip3Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
clip3Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
clip3Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
clip3Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
clip3Display.PolarAxes.EnableCustomRange = 0
clip3Display.PolarAxes.CustomRange = [0.0, 1.0]
clip3Display.PolarAxes.PolarAxisVisibility = 1
clip3Display.PolarAxes.RadialAxesVisibility = 1
clip3Display.PolarAxes.DrawRadialGridlines = 1
clip3Display.PolarAxes.PolarArcsVisibility = 1
clip3Display.PolarAxes.DrawPolarArcsGridlines = 1
clip3Display.PolarAxes.NumberOfRadialAxes = 0
clip3Display.PolarAxes.AutoSubdividePolarAxis = 1
clip3Display.PolarAxes.NumberOfPolarAxis = 0
clip3Display.PolarAxes.MinimumRadius = 0.0
clip3Display.PolarAxes.MinimumAngle = 0.0
clip3Display.PolarAxes.MaximumAngle = 90.0
clip3Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
clip3Display.PolarAxes.Ratio = 1.0
clip3Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
clip3Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
clip3Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
clip3Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
clip3Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
clip3Display.PolarAxes.PolarAxisTitleVisibility = 1
clip3Display.PolarAxes.PolarAxisTitle = 'Radial Distance'
clip3Display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
clip3Display.PolarAxes.PolarLabelVisibility = 1
clip3Display.PolarAxes.PolarLabelFormat = '%-#6.3g'
clip3Display.PolarAxes.PolarLabelExponentLocation = 'Labels'
clip3Display.PolarAxes.RadialLabelVisibility = 1
clip3Display.PolarAxes.RadialLabelFormat = '%-#3.1f'
clip3Display.PolarAxes.RadialLabelLocation = 'Bottom'
clip3Display.PolarAxes.RadialUnitsVisibility = 1
clip3Display.PolarAxes.ScreenSize = 10.0
clip3Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
clip3Display.PolarAxes.PolarAxisTitleOpacity = 1.0
clip3Display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
clip3Display.PolarAxes.PolarAxisTitleFontFile = ''
clip3Display.PolarAxes.PolarAxisTitleBold = 0
clip3Display.PolarAxes.PolarAxisTitleItalic = 0
clip3Display.PolarAxes.PolarAxisTitleShadow = 0
clip3Display.PolarAxes.PolarAxisTitleFontSize = 12
clip3Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
clip3Display.PolarAxes.PolarAxisLabelOpacity = 1.0
clip3Display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
clip3Display.PolarAxes.PolarAxisLabelFontFile = ''
clip3Display.PolarAxes.PolarAxisLabelBold = 0
clip3Display.PolarAxes.PolarAxisLabelItalic = 0
clip3Display.PolarAxes.PolarAxisLabelShadow = 0
clip3Display.PolarAxes.PolarAxisLabelFontSize = 12
clip3Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
clip3Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
clip3Display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
clip3Display.PolarAxes.LastRadialAxisTextFontFile = ''
clip3Display.PolarAxes.LastRadialAxisTextBold = 0
clip3Display.PolarAxes.LastRadialAxisTextItalic = 0
clip3Display.PolarAxes.LastRadialAxisTextShadow = 0
clip3Display.PolarAxes.LastRadialAxisTextFontSize = 12
clip3Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
clip3Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
clip3Display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
clip3Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
clip3Display.PolarAxes.SecondaryRadialAxesTextBold = 0
clip3Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
clip3Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
clip3Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
clip3Display.PolarAxes.EnableDistanceLOD = 1
clip3Display.PolarAxes.DistanceLODThreshold = 0.7
clip3Display.PolarAxes.EnableViewAngleLOD = 1
clip3Display.PolarAxes.ViewAngleLODThreshold = 0.7
clip3Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
clip3Display.PolarAxes.PolarTicksVisibility = 1
clip3Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
clip3Display.PolarAxes.TickLocation = 'Both'
clip3Display.PolarAxes.AxisTickVisibility = 1
clip3Display.PolarAxes.AxisMinorTickVisibility = 0
clip3Display.PolarAxes.ArcTickVisibility = 1
clip3Display.PolarAxes.ArcMinorTickVisibility = 0
clip3Display.PolarAxes.DeltaAngleMajor = 10.0
clip3Display.PolarAxes.DeltaAngleMinor = 5.0
clip3Display.PolarAxes.PolarAxisMajorTickSize = 0.0
clip3Display.PolarAxes.PolarAxisTickRatioSize = 0.3
clip3Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
clip3Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
clip3Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
clip3Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
clip3Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
clip3Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
clip3Display.PolarAxes.ArcMajorTickSize = 0.0
clip3Display.PolarAxes.ArcTickRatioSize = 0.3
clip3Display.PolarAxes.ArcMajorTickThickness = 1.0
clip3Display.PolarAxes.ArcTickRatioThickness = 0.5
clip3Display.PolarAxes.Use2DMode = 0
clip3Display.PolarAxes.UseLogAxis = 0

# hide data in view
Hide(clip2, renderView1)

# show color bar/color legend
clip3Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip3Display, ('POINTS', 'T'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip3Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip3Display.SetScalarBarVisibility(renderView1, True)

# current camera placement for renderView1
renderView1.CameraPosition = [9.999999747378752e-05, 0.00039999998989515007, 
                                                         0.0017392304406046689]
renderView1.CameraFocalPoint = [9.999999747378752e-05, 0.00039999998989515007, 
                                                         9.999999747378752e-05]
renderView1.CameraParallelScale = 0.00042426405799411666


renderView1.ResetCamera()
renderView1.StillRender()

# save screenshot
# SaveScreenshot('./Clip3.png', layout1, ImageResolution=[1980, 907],
#     FontScaling='Scale fonts proportionally',
#     OverrideColorPalette='',
#     StereoMode='No change',
#     TransparentBackground=0, 
#     # PNG options
#     CompressionLevel='5')

# create a new 'Clip'
clip4 = Clip(Input=clip3)
clip4.ClipType = 'Plane'
clip4.Scalars = ['POINTS', 'p']
clip4.Value = 480297.35546875
clip4.Invert = 1
clip4.Crinkleclip = 0
clip4.Exact = 0

# init the 'Plane' selected for 'ClipType'
clip4.ClipType.Origin = [9.999999747378752e-05, 0.0004000000117230229, 
                                                         7.375000132014975e-05]
clip4.ClipType.Normal = [1.0, 0.0, 0.0]
clip4.ClipType.Offset = 0.0

# Properties modified on clip4
clip4.ClipType = 'Scalar'
clip4.Scalars = ['POINTS', 'meltHistory']
clip4.Value = 0.1
clip4.Invert = 0

# show data in view
clip4Display = Show(clip4, renderView1)

# trace defaults for the display properties.
clip4Display.Representation = 'Surface'
clip4Display.AmbientColor = [1.0, 1.0, 1.0]
clip4Display.ColorArrayName = ['POINTS', 'p']
clip4Display.DiffuseColor = [1.0, 1.0, 1.0]
clip4Display.LookupTable = pLUT
clip4Display.MapScalars = 1
clip4Display.MultiComponentsMapping = 0
clip4Display.InterpolateScalarsBeforeMapping = 1
clip4Display.Opacity = 1.0
clip4Display.PointSize = 2.0
clip4Display.LineWidth = 1.0
clip4Display.RenderLinesAsTubes = 0
clip4Display.RenderPointsAsSpheres = 0
clip4Display.Interpolation = 'Gouraud'
clip4Display.Specular = 0.0
clip4Display.SpecularColor = [1.0, 1.0, 1.0]
clip4Display.SpecularPower = 100.0
clip4Display.Luminosity = 0.0
clip4Display.Ambient = 0.0
clip4Display.Diffuse = 1.0
clip4Display.EdgeColor = [0.0, 0.0, 0.5]
clip4Display.BackfaceRepresentation = 'Follow Frontface'
clip4Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
clip4Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
clip4Display.BackfaceOpacity = 1.0
clip4Display.Position = [0.0, 0.0, 0.0]
clip4Display.Scale = [1.0, 1.0, 1.0]
clip4Display.Orientation = [0.0, 0.0, 0.0]
clip4Display.Origin = [0.0, 0.0, 0.0]
clip4Display.Pickable = 1
clip4Display.Texture = None
clip4Display.Triangulate = 0
clip4Display.UseShaderReplacements = 0
clip4Display.ShaderReplacements = ''
clip4Display.NonlinearSubdivisionLevel = 1
clip4Display.UseDataPartitions = 0
clip4Display.OSPRayUseScaleArray = 0
clip4Display.OSPRayScaleArray = 'p'
clip4Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip4Display.OSPRayMaterial = 'None'
clip4Display.Orient = 0
clip4Display.OrientationMode = 'Direction'
clip4Display.SelectOrientationVectors = 'U'
clip4Display.Scaling = 0
clip4Display.ScaleMode = 'No Data Scaling Off'
clip4Display.ScaleFactor = 5.200000159675256e-05
clip4Display.SelectScaleArray = 'p'
clip4Display.GlyphType = 'Arrow'
clip4Display.UseGlyphTable = 0
clip4Display.GlyphTableIndexArray = 'p'
clip4Display.UseCompositeGlyphTable = 0
clip4Display.UseGlyphCullingAndLOD = 0
clip4Display.LODValues = []
clip4Display.ColorByLODIndex = 0
clip4Display.GaussianRadius = 2.6000000798376276e-06
clip4Display.ShaderPreset = 'Sphere'
clip4Display.CustomTriangleScale = 3
clip4Display.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
clip4Display.Emissive = 0
clip4Display.ScaleByArray = 0
clip4Display.SetScaleArray = ['POINTS', 'p']
clip4Display.ScaleArrayComponent = ''
clip4Display.UseScaleFunction = 1
clip4Display.ScaleTransferFunction = 'PiecewiseFunction'
clip4Display.OpacityByArray = 0
clip4Display.OpacityArray = ['POINTS', 'p']
clip4Display.OpacityArrayComponent = ''
clip4Display.OpacityTransferFunction = 'PiecewiseFunction'
clip4Display.DataAxesGrid = 'GridAxesRepresentation'
clip4Display.SelectionCellLabelBold = 0
clip4Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
clip4Display.SelectionCellLabelFontFamily = 'Arial'
clip4Display.SelectionCellLabelFontFile = ''
clip4Display.SelectionCellLabelFontSize = 18
clip4Display.SelectionCellLabelItalic = 0
clip4Display.SelectionCellLabelJustification = 'Left'
clip4Display.SelectionCellLabelOpacity = 1.0
clip4Display.SelectionCellLabelShadow = 0
clip4Display.SelectionPointLabelBold = 0
clip4Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
clip4Display.SelectionPointLabelFontFamily = 'Arial'
clip4Display.SelectionPointLabelFontFile = ''
clip4Display.SelectionPointLabelFontSize = 18
clip4Display.SelectionPointLabelItalic = 0
clip4Display.SelectionPointLabelJustification = 'Left'
clip4Display.SelectionPointLabelOpacity = 1.0
clip4Display.SelectionPointLabelShadow = 0
clip4Display.PolarAxes = 'PolarAxesRepresentation'
clip4Display.ScalarOpacityFunction = pPWF
clip4Display.ScalarOpacityUnitDistance = 1.062582986007341e-05
clip4Display.ExtractedBlockIndex = 1
clip4Display.SelectMapper = 'Projected tetra'
clip4Display.SamplingDimensions = [128, 128, 128]
clip4Display.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip4Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 
                                                                           0.0]
clip4Display.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
clip4Display.GlyphType.TipResolution = 6
clip4Display.GlyphType.TipRadius = 0.1
clip4Display.GlyphType.TipLength = 0.35
clip4Display.GlyphType.ShaftResolution = 6
clip4Display.GlyphType.ShaftRadius = 0.03
clip4Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip4Display.ScaleTransferFunction.Points = [-100115.2890625, 0.0, 0.5, 0.0, 
                                                       477681.0, 1.0, 0.5, 0.0]
clip4Display.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip4Display.OpacityTransferFunction.Points = [-100115.2890625, 0.0, 0.5, 0.0, 
                                                       477681.0, 1.0, 0.5, 0.0]
clip4Display.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip4Display.DataAxesGrid.XTitle = 'X Axis'
clip4Display.DataAxesGrid.YTitle = 'Y Axis'
clip4Display.DataAxesGrid.ZTitle = 'Z Axis'
clip4Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
clip4Display.DataAxesGrid.XTitleFontFamily = 'Arial'
clip4Display.DataAxesGrid.XTitleFontFile = ''
clip4Display.DataAxesGrid.XTitleBold = 0
clip4Display.DataAxesGrid.XTitleItalic = 0
clip4Display.DataAxesGrid.XTitleFontSize = 12
clip4Display.DataAxesGrid.XTitleShadow = 0
clip4Display.DataAxesGrid.XTitleOpacity = 1.0
clip4Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
clip4Display.DataAxesGrid.YTitleFontFamily = 'Arial'
clip4Display.DataAxesGrid.YTitleFontFile = ''
clip4Display.DataAxesGrid.YTitleBold = 0
clip4Display.DataAxesGrid.YTitleItalic = 0
clip4Display.DataAxesGrid.YTitleFontSize = 12
clip4Display.DataAxesGrid.YTitleShadow = 0
clip4Display.DataAxesGrid.YTitleOpacity = 1.0
clip4Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
clip4Display.DataAxesGrid.ZTitleFontFamily = 'Arial'
clip4Display.DataAxesGrid.ZTitleFontFile = ''
clip4Display.DataAxesGrid.ZTitleBold = 0
clip4Display.DataAxesGrid.ZTitleItalic = 0
clip4Display.DataAxesGrid.ZTitleFontSize = 12
clip4Display.DataAxesGrid.ZTitleShadow = 0
clip4Display.DataAxesGrid.ZTitleOpacity = 1.0
clip4Display.DataAxesGrid.FacesToRender = 63
clip4Display.DataAxesGrid.CullBackface = 0
clip4Display.DataAxesGrid.CullFrontface = 1
clip4Display.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
clip4Display.DataAxesGrid.ShowGrid = 0
clip4Display.DataAxesGrid.ShowEdges = 1
clip4Display.DataAxesGrid.ShowTicks = 1
clip4Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
clip4Display.DataAxesGrid.AxesToLabel = 63
clip4Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
clip4Display.DataAxesGrid.XLabelFontFamily = 'Arial'
clip4Display.DataAxesGrid.XLabelFontFile = ''
clip4Display.DataAxesGrid.XLabelBold = 0
clip4Display.DataAxesGrid.XLabelItalic = 0
clip4Display.DataAxesGrid.XLabelFontSize = 12
clip4Display.DataAxesGrid.XLabelShadow = 0
clip4Display.DataAxesGrid.XLabelOpacity = 1.0
clip4Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
clip4Display.DataAxesGrid.YLabelFontFamily = 'Arial'
clip4Display.DataAxesGrid.YLabelFontFile = ''
clip4Display.DataAxesGrid.YLabelBold = 0
clip4Display.DataAxesGrid.YLabelItalic = 0
clip4Display.DataAxesGrid.YLabelFontSize = 12
clip4Display.DataAxesGrid.YLabelShadow = 0
clip4Display.DataAxesGrid.YLabelOpacity = 1.0
clip4Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
clip4Display.DataAxesGrid.ZLabelFontFamily = 'Arial'
clip4Display.DataAxesGrid.ZLabelFontFile = ''
clip4Display.DataAxesGrid.ZLabelBold = 0
clip4Display.DataAxesGrid.ZLabelItalic = 0
clip4Display.DataAxesGrid.ZLabelFontSize = 12
clip4Display.DataAxesGrid.ZLabelShadow = 0
clip4Display.DataAxesGrid.ZLabelOpacity = 1.0
clip4Display.DataAxesGrid.XAxisNotation = 'Mixed'
clip4Display.DataAxesGrid.XAxisPrecision = 2
clip4Display.DataAxesGrid.XAxisUseCustomLabels = 0
clip4Display.DataAxesGrid.XAxisLabels = []
clip4Display.DataAxesGrid.YAxisNotation = 'Mixed'
clip4Display.DataAxesGrid.YAxisPrecision = 2
clip4Display.DataAxesGrid.YAxisUseCustomLabels = 0
clip4Display.DataAxesGrid.YAxisLabels = []
clip4Display.DataAxesGrid.ZAxisNotation = 'Mixed'
clip4Display.DataAxesGrid.ZAxisPrecision = 2
clip4Display.DataAxesGrid.ZAxisUseCustomLabels = 0
clip4Display.DataAxesGrid.ZAxisLabels = []
clip4Display.DataAxesGrid.UseCustomBounds = 0
clip4Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip4Display.PolarAxes.Visibility = 0
clip4Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
clip4Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
clip4Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
clip4Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
clip4Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
clip4Display.PolarAxes.EnableCustomRange = 0
clip4Display.PolarAxes.CustomRange = [0.0, 1.0]
clip4Display.PolarAxes.PolarAxisVisibility = 1
clip4Display.PolarAxes.RadialAxesVisibility = 1
clip4Display.PolarAxes.DrawRadialGridlines = 1
clip4Display.PolarAxes.PolarArcsVisibility = 1
clip4Display.PolarAxes.DrawPolarArcsGridlines = 1
clip4Display.PolarAxes.NumberOfRadialAxes = 0
clip4Display.PolarAxes.AutoSubdividePolarAxis = 1
clip4Display.PolarAxes.NumberOfPolarAxis = 0
clip4Display.PolarAxes.MinimumRadius = 0.0
clip4Display.PolarAxes.MinimumAngle = 0.0
clip4Display.PolarAxes.MaximumAngle = 90.0
clip4Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
clip4Display.PolarAxes.Ratio = 1.0
clip4Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
clip4Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
clip4Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
clip4Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
clip4Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
clip4Display.PolarAxes.PolarAxisTitleVisibility = 1
clip4Display.PolarAxes.PolarAxisTitle = 'Radial Distance'
clip4Display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
clip4Display.PolarAxes.PolarLabelVisibility = 1
clip4Display.PolarAxes.PolarLabelFormat = '%-#6.3g'
clip4Display.PolarAxes.PolarLabelExponentLocation = 'Labels'
clip4Display.PolarAxes.RadialLabelVisibility = 1
clip4Display.PolarAxes.RadialLabelFormat = '%-#3.1f'
clip4Display.PolarAxes.RadialLabelLocation = 'Bottom'
clip4Display.PolarAxes.RadialUnitsVisibility = 1
clip4Display.PolarAxes.ScreenSize = 10.0
clip4Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
clip4Display.PolarAxes.PolarAxisTitleOpacity = 1.0
clip4Display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
clip4Display.PolarAxes.PolarAxisTitleFontFile = ''
clip4Display.PolarAxes.PolarAxisTitleBold = 0
clip4Display.PolarAxes.PolarAxisTitleItalic = 0
clip4Display.PolarAxes.PolarAxisTitleShadow = 0
clip4Display.PolarAxes.PolarAxisTitleFontSize = 12
clip4Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
clip4Display.PolarAxes.PolarAxisLabelOpacity = 1.0
clip4Display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
clip4Display.PolarAxes.PolarAxisLabelFontFile = ''
clip4Display.PolarAxes.PolarAxisLabelBold = 0
clip4Display.PolarAxes.PolarAxisLabelItalic = 0
clip4Display.PolarAxes.PolarAxisLabelShadow = 0
clip4Display.PolarAxes.PolarAxisLabelFontSize = 12
clip4Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
clip4Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
clip4Display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
clip4Display.PolarAxes.LastRadialAxisTextFontFile = ''
clip4Display.PolarAxes.LastRadialAxisTextBold = 0
clip4Display.PolarAxes.LastRadialAxisTextItalic = 0
clip4Display.PolarAxes.LastRadialAxisTextShadow = 0
clip4Display.PolarAxes.LastRadialAxisTextFontSize = 12
clip4Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
clip4Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
clip4Display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
clip4Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
clip4Display.PolarAxes.SecondaryRadialAxesTextBold = 0
clip4Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
clip4Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
clip4Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
clip4Display.PolarAxes.EnableDistanceLOD = 1
clip4Display.PolarAxes.DistanceLODThreshold = 0.7
clip4Display.PolarAxes.EnableViewAngleLOD = 1
clip4Display.PolarAxes.ViewAngleLODThreshold = 0.7
clip4Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
clip4Display.PolarAxes.PolarTicksVisibility = 1
clip4Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
clip4Display.PolarAxes.TickLocation = 'Both'
clip4Display.PolarAxes.AxisTickVisibility = 1
clip4Display.PolarAxes.AxisMinorTickVisibility = 0
clip4Display.PolarAxes.ArcTickVisibility = 1
clip4Display.PolarAxes.ArcMinorTickVisibility = 0
clip4Display.PolarAxes.DeltaAngleMajor = 10.0
clip4Display.PolarAxes.DeltaAngleMinor = 5.0
clip4Display.PolarAxes.PolarAxisMajorTickSize = 0.0
clip4Display.PolarAxes.PolarAxisTickRatioSize = 0.3
clip4Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
clip4Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
clip4Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
clip4Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
clip4Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
clip4Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
clip4Display.PolarAxes.ArcMajorTickSize = 0.0
clip4Display.PolarAxes.ArcTickRatioSize = 0.3
clip4Display.PolarAxes.ArcMajorTickThickness = 1.0
clip4Display.PolarAxes.ArcTickRatioThickness = 0.5
clip4Display.PolarAxes.Use2DMode = 0
clip4Display.PolarAxes.UseLogAxis = 0

# hide data in view
Hide(clip3, renderView1)

# show color bar/color legend
clip4Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip4Display, ('POINTS', 'T'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip4Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip4Display.SetScalarBarVisibility(renderView1, True)

# current camera placement for renderView1
renderView1.CameraPosition = [9.999999747378752e-05, 0.00039999998989515007, 
                                                         0.0017392304406046689]
renderView1.CameraFocalPoint = [9.999999747378752e-05, 0.00039999998989515007, 
                                                         9.999999747378752e-05]
renderView1.CameraParallelScale = 0.00042426405799411666


renderView1.ResetCamera()
renderView1.StillRender()

# save screenshot
# SaveScreenshot('./Clip4.png', layout1, ImageResolution=[1980, 907],
#     FontScaling='Scale fonts proportionally',
#     OverrideColorPalette='',
#     StereoMode='No change',
#     TransparentBackground=0, 
#     # PNG options
#     CompressionLevel='5')

# create a new 'Clip'
clip5 = Clip(Input=clip4)
clip5.ClipType = 'Plane'
clip5.Scalars = ['POINTS', 'p']
clip5.Value = 188782.85546875
clip5.Invert = 1
clip5.Crinkleclip = 0
clip5.Exact = 0

# init the 'Plane' selected for 'ClipType'
clip5.ClipType.Origin = [9.625130587664898e-05, 0.0004000000117230229, 
                                                        0.00011250020179431885]
clip5.ClipType.Normal = [1.0, 0.0, 0.0]
clip5.ClipType.Offset = 0.0

# Properties modified on clip5
clip5.ClipType = 'Scalar'
clip5.Scalars = ['POINTS', 'alpha.metal']
clip5.Value = 0.9999
clip5.Invert = 0

# show data in view
clip5Display = Show(clip5, renderView1)

# trace defaults for the display properties.
clip5Display.Representation = 'Surface'
clip5Display.AmbientColor = [1.0, 1.0, 1.0]
clip5Display.ColorArrayName = ['POINTS', 'p']
clip5Display.DiffuseColor = [1.0, 1.0, 1.0]
clip5Display.LookupTable = pLUT
clip5Display.MapScalars = 1
clip5Display.MultiComponentsMapping = 0
clip5Display.InterpolateScalarsBeforeMapping = 1
clip5Display.Opacity = 1.0
clip5Display.PointSize = 2.0
clip5Display.LineWidth = 1.0
clip5Display.RenderLinesAsTubes = 0
clip5Display.RenderPointsAsSpheres = 0
clip5Display.Interpolation = 'Gouraud'
clip5Display.Specular = 0.0
clip5Display.SpecularColor = [1.0, 1.0, 1.0]
clip5Display.SpecularPower = 100.0
clip5Display.Luminosity = 0.0
clip5Display.Ambient = 0.0
clip5Display.Diffuse = 1.0
clip5Display.EdgeColor = [0.0, 0.0, 0.5]
clip5Display.BackfaceRepresentation = 'Follow Frontface'
clip5Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
clip5Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
clip5Display.BackfaceOpacity = 1.0
clip5Display.Position = [0.0, 0.0, 0.0]
clip5Display.Scale = [1.0, 1.0, 1.0]
clip5Display.Orientation = [0.0, 0.0, 0.0]
clip5Display.Origin = [0.0, 0.0, 0.0]
clip5Display.Pickable = 1
clip5Display.Texture = None
clip5Display.Triangulate = 0
clip5Display.UseShaderReplacements = 0
clip5Display.ShaderReplacements = ''
clip5Display.NonlinearSubdivisionLevel = 1
clip5Display.UseDataPartitions = 0
clip5Display.OSPRayUseScaleArray = 0
clip5Display.OSPRayScaleArray = 'p'
clip5Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip5Display.OSPRayMaterial = 'None'
clip5Display.Orient = 0
clip5Display.OrientationMode = 'Direction'
clip5Display.SelectOrientationVectors = 'U'
clip5Display.Scaling = 0
clip5Display.ScaleMode = 'No Data Scaling Off'
clip5Display.ScaleFactor = 5.200000159675256e-05
clip5Display.SelectScaleArray = 'p'
clip5Display.GlyphType = 'Arrow'
clip5Display.UseGlyphTable = 0
clip5Display.GlyphTableIndexArray = 'p'
clip5Display.UseCompositeGlyphTable = 0
clip5Display.UseGlyphCullingAndLOD = 0
clip5Display.LODValues = []
clip5Display.ColorByLODIndex = 0
clip5Display.GaussianRadius = 2.6000000798376276e-06
clip5Display.ShaderPreset = 'Sphere'
clip5Display.CustomTriangleScale = 3
clip5Display.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
clip5Display.Emissive = 0
clip5Display.ScaleByArray = 0
clip5Display.SetScaleArray = ['POINTS', 'p']
clip5Display.ScaleArrayComponent = ''
clip5Display.UseScaleFunction = 1
clip5Display.ScaleTransferFunction = 'PiecewiseFunction'
clip5Display.OpacityByArray = 0
clip5Display.OpacityArray = ['POINTS', 'p']
clip5Display.OpacityArrayComponent = ''
clip5Display.OpacityTransferFunction = 'PiecewiseFunction'
clip5Display.DataAxesGrid = 'GridAxesRepresentation'
clip5Display.SelectionCellLabelBold = 0
clip5Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
clip5Display.SelectionCellLabelFontFamily = 'Arial'
clip5Display.SelectionCellLabelFontFile = ''
clip5Display.SelectionCellLabelFontSize = 18
clip5Display.SelectionCellLabelItalic = 0
clip5Display.SelectionCellLabelJustification = 'Left'
clip5Display.SelectionCellLabelOpacity = 1.0
clip5Display.SelectionCellLabelShadow = 0
clip5Display.SelectionPointLabelBold = 0
clip5Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
clip5Display.SelectionPointLabelFontFamily = 'Arial'
clip5Display.SelectionPointLabelFontFile = ''
clip5Display.SelectionPointLabelFontSize = 18
clip5Display.SelectionPointLabelItalic = 0
clip5Display.SelectionPointLabelJustification = 'Left'
clip5Display.SelectionPointLabelOpacity = 1.0
clip5Display.SelectionPointLabelShadow = 0
clip5Display.PolarAxes = 'PolarAxesRepresentation'
clip5Display.ScalarOpacityFunction = pPWF
clip5Display.ScalarOpacityUnitDistance = 1.0906548326314667e-05
clip5Display.ExtractedBlockIndex = 1
clip5Display.SelectMapper = 'Projected tetra'
clip5Display.SamplingDimensions = [128, 128, 128]
clip5Display.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip5Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 
                                                                           0.0]
clip5Display.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
clip5Display.GlyphType.TipResolution = 6
clip5Display.GlyphType.TipRadius = 0.1
clip5Display.GlyphType.TipLength = 0.35
clip5Display.GlyphType.ShaftResolution = 6
clip5Display.GlyphType.ShaftRadius = 0.03
clip5Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip5Display.ScaleTransferFunction.Points = [-67124.1875, 0.0, 0.5, 0.0, 
                                                    294047.9375, 1.0, 0.5, 0.0]
clip5Display.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip5Display.OpacityTransferFunction.Points = [-67124.1875, 0.0, 0.5, 0.0, 
                                                    294047.9375, 1.0, 0.5, 0.0]
clip5Display.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip5Display.DataAxesGrid.XTitle = 'X Axis'
clip5Display.DataAxesGrid.YTitle = 'Y Axis'
clip5Display.DataAxesGrid.ZTitle = 'Z Axis'
clip5Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
clip5Display.DataAxesGrid.XTitleFontFamily = 'Arial'
clip5Display.DataAxesGrid.XTitleFontFile = ''
clip5Display.DataAxesGrid.XTitleBold = 0
clip5Display.DataAxesGrid.XTitleItalic = 0
clip5Display.DataAxesGrid.XTitleFontSize = 12
clip5Display.DataAxesGrid.XTitleShadow = 0
clip5Display.DataAxesGrid.XTitleOpacity = 1.0
clip5Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
clip5Display.DataAxesGrid.YTitleFontFamily = 'Arial'
clip5Display.DataAxesGrid.YTitleFontFile = ''
clip5Display.DataAxesGrid.YTitleBold = 0
clip5Display.DataAxesGrid.YTitleItalic = 0
clip5Display.DataAxesGrid.YTitleFontSize = 12
clip5Display.DataAxesGrid.YTitleShadow = 0
clip5Display.DataAxesGrid.YTitleOpacity = 1.0
clip5Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
clip5Display.DataAxesGrid.ZTitleFontFamily = 'Arial'
clip5Display.DataAxesGrid.ZTitleFontFile = ''
clip5Display.DataAxesGrid.ZTitleBold = 0
clip5Display.DataAxesGrid.ZTitleItalic = 0
clip5Display.DataAxesGrid.ZTitleFontSize = 12
clip5Display.DataAxesGrid.ZTitleShadow = 0
clip5Display.DataAxesGrid.ZTitleOpacity = 1.0
clip5Display.DataAxesGrid.FacesToRender = 63
clip5Display.DataAxesGrid.CullBackface = 0
clip5Display.DataAxesGrid.CullFrontface = 1
clip5Display.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
clip5Display.DataAxesGrid.ShowGrid = 0
clip5Display.DataAxesGrid.ShowEdges = 1
clip5Display.DataAxesGrid.ShowTicks = 1
clip5Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
clip5Display.DataAxesGrid.AxesToLabel = 63
clip5Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
clip5Display.DataAxesGrid.XLabelFontFamily = 'Arial'
clip5Display.DataAxesGrid.XLabelFontFile = ''
clip5Display.DataAxesGrid.XLabelBold = 0
clip5Display.DataAxesGrid.XLabelItalic = 0
clip5Display.DataAxesGrid.XLabelFontSize = 12
clip5Display.DataAxesGrid.XLabelShadow = 0
clip5Display.DataAxesGrid.XLabelOpacity = 1.0
clip5Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
clip5Display.DataAxesGrid.YLabelFontFamily = 'Arial'
clip5Display.DataAxesGrid.YLabelFontFile = ''
clip5Display.DataAxesGrid.YLabelBold = 0
clip5Display.DataAxesGrid.YLabelItalic = 0
clip5Display.DataAxesGrid.YLabelFontSize = 12
clip5Display.DataAxesGrid.YLabelShadow = 0
clip5Display.DataAxesGrid.YLabelOpacity = 1.0
clip5Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
clip5Display.DataAxesGrid.ZLabelFontFamily = 'Arial'
clip5Display.DataAxesGrid.ZLabelFontFile = ''
clip5Display.DataAxesGrid.ZLabelBold = 0
clip5Display.DataAxesGrid.ZLabelItalic = 0
clip5Display.DataAxesGrid.ZLabelFontSize = 12
clip5Display.DataAxesGrid.ZLabelShadow = 0
clip5Display.DataAxesGrid.ZLabelOpacity = 1.0
clip5Display.DataAxesGrid.XAxisNotation = 'Mixed'
clip5Display.DataAxesGrid.XAxisPrecision = 2
clip5Display.DataAxesGrid.XAxisUseCustomLabels = 0
clip5Display.DataAxesGrid.XAxisLabels = []
clip5Display.DataAxesGrid.YAxisNotation = 'Mixed'
clip5Display.DataAxesGrid.YAxisPrecision = 2
clip5Display.DataAxesGrid.YAxisUseCustomLabels = 0
clip5Display.DataAxesGrid.YAxisLabels = []
clip5Display.DataAxesGrid.ZAxisNotation = 'Mixed'
clip5Display.DataAxesGrid.ZAxisPrecision = 2
clip5Display.DataAxesGrid.ZAxisUseCustomLabels = 0
clip5Display.DataAxesGrid.ZAxisLabels = []
clip5Display.DataAxesGrid.UseCustomBounds = 0
clip5Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip5Display.PolarAxes.Visibility = 0
clip5Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
clip5Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
clip5Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
clip5Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
clip5Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
clip5Display.PolarAxes.EnableCustomRange = 0
clip5Display.PolarAxes.CustomRange = [0.0, 1.0]
clip5Display.PolarAxes.PolarAxisVisibility = 1
clip5Display.PolarAxes.RadialAxesVisibility = 1
clip5Display.PolarAxes.DrawRadialGridlines = 1
clip5Display.PolarAxes.PolarArcsVisibility = 1
clip5Display.PolarAxes.DrawPolarArcsGridlines = 1
clip5Display.PolarAxes.NumberOfRadialAxes = 0
clip5Display.PolarAxes.AutoSubdividePolarAxis = 1
clip5Display.PolarAxes.NumberOfPolarAxis = 0
clip5Display.PolarAxes.MinimumRadius = 0.0
clip5Display.PolarAxes.MinimumAngle = 0.0
clip5Display.PolarAxes.MaximumAngle = 90.0
clip5Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
clip5Display.PolarAxes.Ratio = 1.0
clip5Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
clip5Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
clip5Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
clip5Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
clip5Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
clip5Display.PolarAxes.PolarAxisTitleVisibility = 1
clip5Display.PolarAxes.PolarAxisTitle = 'Radial Distance'
clip5Display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
clip5Display.PolarAxes.PolarLabelVisibility = 1
clip5Display.PolarAxes.PolarLabelFormat = '%-#6.3g'
clip5Display.PolarAxes.PolarLabelExponentLocation = 'Labels'
clip5Display.PolarAxes.RadialLabelVisibility = 1
clip5Display.PolarAxes.RadialLabelFormat = '%-#3.1f'
clip5Display.PolarAxes.RadialLabelLocation = 'Bottom'
clip5Display.PolarAxes.RadialUnitsVisibility = 1
clip5Display.PolarAxes.ScreenSize = 10.0
clip5Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
clip5Display.PolarAxes.PolarAxisTitleOpacity = 1.0
clip5Display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
clip5Display.PolarAxes.PolarAxisTitleFontFile = ''
clip5Display.PolarAxes.PolarAxisTitleBold = 0
clip5Display.PolarAxes.PolarAxisTitleItalic = 0
clip5Display.PolarAxes.PolarAxisTitleShadow = 0
clip5Display.PolarAxes.PolarAxisTitleFontSize = 12
clip5Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
clip5Display.PolarAxes.PolarAxisLabelOpacity = 1.0
clip5Display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
clip5Display.PolarAxes.PolarAxisLabelFontFile = ''
clip5Display.PolarAxes.PolarAxisLabelBold = 0
clip5Display.PolarAxes.PolarAxisLabelItalic = 0
clip5Display.PolarAxes.PolarAxisLabelShadow = 0
clip5Display.PolarAxes.PolarAxisLabelFontSize = 12
clip5Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
clip5Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
clip5Display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
clip5Display.PolarAxes.LastRadialAxisTextFontFile = ''
clip5Display.PolarAxes.LastRadialAxisTextBold = 0
clip5Display.PolarAxes.LastRadialAxisTextItalic = 0
clip5Display.PolarAxes.LastRadialAxisTextShadow = 0
clip5Display.PolarAxes.LastRadialAxisTextFontSize = 12
clip5Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
clip5Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
clip5Display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
clip5Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
clip5Display.PolarAxes.SecondaryRadialAxesTextBold = 0
clip5Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
clip5Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
clip5Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
clip5Display.PolarAxes.EnableDistanceLOD = 1
clip5Display.PolarAxes.DistanceLODThreshold = 0.7
clip5Display.PolarAxes.EnableViewAngleLOD = 1
clip5Display.PolarAxes.ViewAngleLODThreshold = 0.7
clip5Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
clip5Display.PolarAxes.PolarTicksVisibility = 1
clip5Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
clip5Display.PolarAxes.TickLocation = 'Both'
clip5Display.PolarAxes.AxisTickVisibility = 1
clip5Display.PolarAxes.AxisMinorTickVisibility = 0
clip5Display.PolarAxes.ArcTickVisibility = 1
clip5Display.PolarAxes.ArcMinorTickVisibility = 0
clip5Display.PolarAxes.DeltaAngleMajor = 10.0
clip5Display.PolarAxes.DeltaAngleMinor = 5.0
clip5Display.PolarAxes.PolarAxisMajorTickSize = 0.0
clip5Display.PolarAxes.PolarAxisTickRatioSize = 0.3
clip5Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
clip5Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
clip5Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
clip5Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
clip5Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
clip5Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
clip5Display.PolarAxes.ArcMajorTickSize = 0.0
clip5Display.PolarAxes.ArcTickRatioSize = 0.3
clip5Display.PolarAxes.ArcMajorTickThickness = 1.0
clip5Display.PolarAxes.ArcTickRatioThickness = 0.5
clip5Display.PolarAxes.Use2DMode = 0
clip5Display.PolarAxes.UseLogAxis = 0

# hide data in view
Hide(clip4, renderView1)

# show color bar/color legend
clip5Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip5Display, ('POINTS', 'T'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip5Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip5Display.SetScalarBarVisibility(renderView1, True)

# current camera placement for renderView1
renderView1.CameraPosition = [9.999999747378752e-05, 0.00039999998989515007, 
                                                         0.0017392304406046689]
renderView1.CameraFocalPoint = [9.999999747378752e-05, 0.00039999998989515007, 
                                                         9.999999747378752e-05]
renderView1.CameraParallelScale = 0.00042426405799411666


renderView1.ResetCamera()
renderView1.StillRender()

# # save screenshot
# SaveScreenshot('./Clip5.png', layout1, ImageResolution=[1980, 907],
#     FontScaling='Scale fonts proportionally',
#     OverrideColorPalette='',
#     StereoMode='No change',
#     TransparentBackground=0, 
#     # PNG options
#     CompressionLevel='5')

# destroy renderView1
Delete(renderView1)
del renderView1

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.UseCache = 0
spreadSheetView1.ViewSize = [400, 400]
spreadSheetView1.CellFontSize = 9
spreadSheetView1.HeaderFontSize = 9
spreadSheetView1.SelectionOnly = 0
spreadSheetView1.GenerateCellConnectivity = 0
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.InvertOrder = 0
spreadSheetView1.BlockSize = 1024
spreadSheetView1.HiddenColumnLabels = []
spreadSheetView1.FieldAssociation = 'Point Data'

# show data in view
clip5Display = Show(clip5, spreadSheetView1)

# trace defaults for the display properties.
clip5Display.CompositeDataSetIndex = [0]

# assign view to a particular cell in the layout
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=0)

# export view
ExportView('./meltpool.csv', view=spreadSheetView1)

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).