# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XDMF Reader'
visuxdmf = XDMFReader(registrationName='visu.xdmf', FileNames=['arquivos/visu.xdmf'])

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# Adjust camera

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(visuxdmf)

# show data in view
visuxdmfDisplay = Show(visuxdmf, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
visuxdmfDisplay.Representation = 'Outline'
visuxdmfDisplay.ColorArrayName = ['CELLS', '']
visuxdmfDisplay.SelectTCoordArray = 'None'
visuxdmfDisplay.SelectNormalArray = 'None'
visuxdmfDisplay.SelectTangentArray = 'None'
visuxdmfDisplay.OSPRayScaleArray = 'V'
visuxdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
visuxdmfDisplay.SelectOrientationVectors = 'V'
visuxdmfDisplay.ScaleFactor = 0.5
visuxdmfDisplay.SelectScaleArray = 'level_set'
visuxdmfDisplay.GlyphType = 'Arrow'
visuxdmfDisplay.GlyphTableIndexArray = 'level_set'
visuxdmfDisplay.GaussianRadius = 0.025
visuxdmfDisplay.SetScaleArray = ['POINTS', 'V']
visuxdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
visuxdmfDisplay.OpacityArray = ['POINTS', 'V']
visuxdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
visuxdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
visuxdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
visuxdmfDisplay.ScalarOpacityUnitDistance = 0.13093483631660657
visuxdmfDisplay.SelectInputVectors = ['POINTS', 'V']
visuxdmfDisplay.WriteLog = ''

# get the material library
materialLibrary1 = GetMaterialLibrary()



# show data in view
visuxdmfDisplay = Show(visuxdmf, renderView1, 'StructuredGridRepresentation')

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=visuxdmf)
cellDatatoPointData1.CellDataArraytoprocess = ['level_set', 'obstaculo', 'pressure', 'visc', 'vorticity']

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Outline'
cellDatatoPointData1Display.ColorArrayName = ['POINTS', '']
cellDatatoPointData1Display.SelectTCoordArray = 'None'
cellDatatoPointData1Display.SelectNormalArray = 'None'
cellDatatoPointData1Display.SelectTangentArray = 'None'
cellDatatoPointData1Display.OSPRayScaleArray = 'level_set'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.SelectOrientationVectors = 'vorticity'
cellDatatoPointData1Display.ScaleFactor = 0.5
cellDatatoPointData1Display.SelectScaleArray = 'level_set'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'level_set'
cellDatatoPointData1Display.GaussianRadius = 0.025
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'level_set']
cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'level_set']
cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'
cellDatatoPointData1Display.ScalarOpacityUnitDistance = 0.13093483631660657
cellDatatoPointData1Display.SelectInputVectors = ['POINTS', 'vorticity']
cellDatatoPointData1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cellDatatoPointData1Display.ScaleTransferFunction.Points = [-0.49000000953674316, 0.0, 0.5, 0.0, 0.49000000953674316, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cellDatatoPointData1Display.OpacityTransferFunction.Points = [-0.49000000953674316, 0.0, 0.5, 0.0, 0.49000000953674316, 1.0, 0.5, 0.0]


# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=cellDatatoPointData1)
clip1.ClipType = 'Scalar'
clip1.Value = 0.0
clip1.Invert = 0

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'level_set'
level_setLUT = GetColorTransferFunction('level_set')

# get opacity transfer function/opacity map for 'level_set'
level_setPWF = GetOpacityTransferFunction('level_set')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['POINTS', 'level_set']
clip1Display.LookupTable = level_setLUT
clip1Display.SelectTCoordArray = 'None'
clip1Display.SelectNormalArray = 'None'
clip1Display.SelectTangentArray = 'None'
clip1Display.OSPRayScaleArray = 'level_set'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'vorticity'
clip1Display.ScaleFactor = 0.5
clip1Display.SelectScaleArray = 'level_set'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'level_set'
clip1Display.GaussianRadius = 0.025
clip1Display.SetScaleArray = ['POINTS', 'level_set']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'level_set']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = level_setPWF
clip1Display.ScalarOpacityUnitDistance = 0.1605474949353342
clip1Display.OpacityArrayName = ['POINTS', 'level_set']
clip1Display.SelectInputVectors = ['POINTS', 'vorticity']
clip1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.49000000953674316, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.49000000953674316, 1.0, 0.5, 0.0]

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)


# get 2D transfer function for 'level_set'
level_setTF2D = GetTransferFunction2D('level_set')

# current camera placement for renderView1
renderView1.CameraPosition = [0.5, 12.538195644853696, 0.5]
renderView1.CameraFocalPoint = [0.5, 2.5, 0.5]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 2.598076211353316

# reset view to fit data
renderView1.ResetCamera(False)

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'V', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(level_setLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'V'
vLUT = GetColorTransferFunction('V')

# get opacity transfer function/opacity map for 'V'
vPWF = GetOpacityTransferFunction('V')

# get 2D transfer function for 'V'
vTF2D = GetTransferFunction2D('V')

# turn off scalar coloring
ColorBy(clip1Display, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(vLUT, renderView1)

# change solid color
clip1Display.AmbientColor = [0.18823529411764706, 0.6745098039215687, 1.0]
clip1Display.DiffuseColor = [0.18823529411764706, 0.6745098039215687, 1.0]

# Properties modified on clip1Display
clip1Display.Opacity = 0.5

renderView1.ResetActiveCameraToPositiveY()

# reset view to fit data
renderView1.ResetCamera(False)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(984, 477)
