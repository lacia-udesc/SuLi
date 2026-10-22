# trace generated using paraview version 5.11.1-RC1
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

# show data in view
visuxdmfDisplay = Show(visuxdmf, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
visuxdmfDisplay.Representation = 'Outline'

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [0.25, 0.4449999928474426, 2.1553751444033646]
renderView1.CameraFocalPoint = [0.25, 0.4449999928474426, 0.125]
renderView1.CameraParallelScale = 0.5254997560743715
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [0.25, 0.4449999928474426, 2.1553751444033646]
renderView1.CameraFocalPoint = [0.25, 0.4449999928474426, 0.125]
renderView1.CameraParallelScale = 0.5254997560743715
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=visuxdmf)
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# Properties modified on clip1
clip1.ClipType = 'Scalar'
clip1.Value = 0.0

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# get 2D transfer function for 'obstaculo'
obstaculoTF2D = GetTransferFunction2D('obstaculo')

# get color transfer function/color map for 'obstaculo'
obstaculoLUT = GetColorTransferFunction('obstaculo')
obstaculoLUT.TransferFunction2D = obstaculoTF2D
obstaculoLUT.RGBPoints = [-0.25, 0.231373, 0.298039, 0.752941, -0.1249999999041285, 0.865003, 0.865003, 0.865003, 1.9174301035818075e-10, 0.705882, 0.0156863, 0.14902]
obstaculoLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'obstaculo'
obstaculoPWF = GetOpacityTransferFunction('obstaculo')
obstaculoPWF.Points = [-0.25, 0.0, 0.5, 0.0, 1.9174301035818075e-10, 1.0, 0.5, 0.0]
obstaculoPWF.ScalarRangeInitialized = 1
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# Properties modified on clip1
clip1.Invert = 0

# update the view to ensure updated data information
renderView1.Update()

# Rescale transfer function
obstaculoLUT.RescaleTransferFunction(-0.25, 0.15538574755191803)

# Rescale transfer function
obstaculoPWF.RescaleTransferFunction(-0.25, 0.15538574755191803)
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# turn off scalar coloring
ColorBy(clip1Display, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(obstaculoLUT, renderView1)
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# turn off scalar coloring
ColorBy(clip1Display, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(obstaculoLUT, renderView1)
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# set active source
SetActiveSource(visuxdmf)
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=visuxdmf)

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Outline'

# hide data in view
Hide(visuxdmf, renderView1)

# update the view to ensure updated data information
renderView1.Update()
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# create a new 'Clip'
clip2 = Clip(registrationName='Clip2', Input=cellDatatoPointData1)
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# Properties modified on clip2
clip2.ClipType = 'Scalar'
clip2.Value = 0.0

# show data in view
clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'

# show color bar/color legend
clip2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# get 2D transfer function for 'levelset'
levelsetTF2D = GetTransferFunction2D('levelset')

# get color transfer function/color map for 'levelset'
levelsetLUT = GetColorTransferFunction('levelset')
levelsetLUT.TransferFunction2D = levelsetTF2D
levelsetLUT.RGBPoints = [-0.0925000011920929, 0.231373, 0.298039, 0.752941, -0.04625000048149876, 0.865003, 0.865003, 0.865003, 2.2909538155424514e-10, 0.705882, 0.0156863, 0.14902]
levelsetLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'levelset'
levelsetPWF = GetOpacityTransferFunction('levelset')
levelsetPWF.Points = [-0.0925000011920929, 0.0, 0.5, 0.0, 2.2909538155424514e-10, 1.0, 0.5, 0.0]
levelsetPWF.ScalarRangeInitialized = 1
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# Properties modified on clip2
clip2.Invert = 0

# update the view to ensure updated data information
renderView1.Update()

# Rescale transfer function
levelsetLUT.RescaleTransferFunction(-0.0925000011920929, 0.15361972153186798)

# Rescale transfer function
levelsetPWF.RescaleTransferFunction(-0.0925000011920929, 0.15361972153186798)
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

# set scalar coloring
ColorBy(clip2Display, ('POINTS', 'V', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(levelsetLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip2Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip2Display.SetScalarBarVisibility(renderView1, True)

# get 2D transfer function for 'V'
vTF2D = GetTransferFunction2D('V')

# get color transfer function/color map for 'V'
vLUT = GetColorTransferFunction('V')
vLUT.TransferFunction2D = vTF2D
vLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.0946240467429111, 0.865003, 0.865003, 0.865003, 0.1892480934858222, 0.705882, 0.0156863, 0.14902]
vLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'V'
vPWF = GetOpacityTransferFunction('V')
vPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.1892480934858222, 1.0, 0.5, 0.0]
vPWF.ScalarRangeInitialized = 1
# Adjust camera

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1550, 789)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [2.009623947391358, -0.4778974942573588, 0.5426206640882948]
renderView1.CameraFocalPoint = [0.25000000000000033, 0.44499999284744257, 0.12499999999999994]
renderView1.CameraViewUp = [-0.19365588983101847, 0.07346747905334328, 0.978314839842115]
renderView1.CameraParallelScale = 0.5254997560743715

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
