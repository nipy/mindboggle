#!/usr/bin/env python
"""
VTK Viewer
Written 2012-2013 Hal Canary <http://cs.unc.edu/~hal>
Copyright 2012-2013 University of North Carolina at Chapel Hill.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   LICENSE.md in this repository or
   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
implied.  See the License for the specific language governing
permissions and limitations under the License.

CONTROLS:

left mouse button:  rotation
right mouse button:  zooming
middle mouse button:  panning
ctrl + left mouse button:  spinning

ctrl + shift + left mouse button:  zooming
shift + left mouse button:  panning

j:  joystick
t:  trackball

c:  camera mode
a:  actor mode

3: toggle  stereo mode

e: exit the application
e: quit the application

f: fly to the picked point
r: reset the camera

p: pick

s: surface representation
w: wireframe representation

Set the STEREO_TYPE environment variable to control stereo type.
    STEREO_TYPE=CRYSTAL_EYES
    STEREO_TYPE=RED_BLUE
    STEREO_TYPE=INTERLACED
    STEREO_TYPE=LEFT
    STEREO_TYPE=RIGHT
    STEREO_TYPE=DRESDEN
    STEREO_TYPE=ANAGLYPH
    STEREO_TYPE=CHECKERBOARD
    STEREO_TYPE=SPLITVIEWPORT_HORIZONTAL

Set the COLORMAP environment variable to change the scalar color map.
It should be in the format of a ParaView-style xml colormap file.

"""
import os
import vtk
import xml.etree.ElementTree


class VTKViewer(object):
    def __init__(self):
        self.renWin = vtk.vtkRenderWindow()
        self.renWin.SetSize(800, 600)
        self.renWin.SetWindowName("VTK Viewer")
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetRenderWindow(self.renWin)
        self.renderer = vtk.vtkRenderer()
        self.renWin.AddRenderer(self.renderer)
        interactorStyle = self.iren.GetInteractorStyle()
        interactorStyle.SetCurrentStyleToTrackballCamera()
        if "STEREO_TYPE" in os.environ:
            vtkStereoType = VTKViewer.GetVTKStereoType(
                os.environ["STEREO_TYPE"])
            if vtkStereoType is not None:
                self.renWin.SetStereoType(vtkStereoType)
            else:
                print('vtkStereoType is None') #.format(stereoType)

    def Start(self):
        self.renWin.Render()
        self.iren.Start()

    @staticmethod
    def GetDefaultColorMap(dataRange):
        colorMap = vtk.vtkColorTransferFunction()
        colorMap.SetColorSpaceToLab()
        colorMap.AddRGBPoint(dataRange[0], 0.865, 0.865, 0.865)
        colorMap.AddRGBPoint(dataRange[1], 0.706, 0.016, 0.150)
        colorMap.Build()
        return colorMap

    @staticmethod
    def LoadColorMap(file_name):
        """
        ParaView has an XML colormap format:

        <ColorMap space="RGB">
        <Point x="0.0"  r="0.0"      g="0.0"      b="0.0"/>
        <Point x="0.4"  r="0.901961" g="0.0"      b="0.0"/>
        <Point x="0.8"  r="0.901961" g="0.901961" b="0.0"/>
        <Point x="1.0"  r="1.0"      g="1.0"      b="1.0"/>
        <NaN            r="0.0"      g="0.498039" b="1.0"/>
        </ColorMap>

        """

        colorMap = vtk.vtkColorTransferFunction()
        root = xml.etree.ElementTree.parse(file_name).getroot()
        if root.tag != "ColorMap":
            raise Exception('Wrong Kind of XML File')
            return None
        if "space" in root.attrib:
            space = root.attrib["space"]
            if space == "RGB":
                colorMap.SetColorSpaceToRGB()
            elif space == "Lab":
                colorMap.SetColorSpaceToLab()
            elif space == "Wrapped":
                colorMap.SetColorSpaceToHSV()
                colorMap.HSVWrapOn()
            elif space == "Diverging":
                colorMap.SetColorSpaceToDiverging()
            else:
                colorMap.SetColorSpaceToHSV()
        else:
            colorMap.SetColorSpaceToHSV()
        point_found = False
        for point in root:
            if point.tag == "Point":
                point_found = True
                a, r, g, b, h, s, v = 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                if 'x' not in point.attrib:
                    continue
                x = float(point.attrib['x'])
                if 'o' in point.attrib:
                    a = float(point.attrib['o'])
                if 'r' in point.attrib:
                    if 'g' not in point.attrib or 'b' not in point.attrib:
                        continue
                    r = float(point.attrib['r'])
                    g = float(point.attrib['g'])
                    b = float(point.attrib['b'])
                    colorMap.AddRGBPoint(x, r, g, b)
                elif 'h' in point.attrib:
                    if 's' not in point.attrib or 'v' not in point.attrib:
                        continue
                    h = float(point.attrib['h'])
                    s = float(point.attrib['s'])
                    v = float(point.attrib['v'])
                    colorMap.AddHSVPoint(x, h, s, v)
                else: ## 'r' or 'h' required
                    continue
            elif point.tag == "NaN":
                r, g, b = 0.25, 0.0, 0.0
                if 'r' in point.attrib:
                    r = float(point.attrib['r'])
                if 'g' in point.attrib:
                    g = float(point.attrib['g'])
                if 'b' in point.attrib:
                    b = float(point.attrib['b'])
                colorMap.SetNanColor(r, g, b)
            ## NAN doesn't support HSV.  Why not?
        if not point_found:
            return None
        colorMap.Build()
        return colorMap

    def AddPolyData(self, polyData, colorMap=None):
        """
        colorMap should be a vtkScalarsToColors (or derived class) object
        """
        if colorMap is None:
            colorMap = VTKViewer.GetDefaultColorMap(polyData.GetScalarRange())
        polyDataMapper = vtk.vtkPolyDataMapper()
        polyDataMapper.SetLookupTable(colorMap)

        if polyData.GetPointData().GetNormals() is None:
            polyDataNormals = vtk.vtkPolyDataNormals()

            # Migrate to VTK6:
            # http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
            # Old: polyDataNormals.SetInput(polyData)
            polyDataNormals.SetInputData(polyData)

            polyDataNormals.SetFeatureAngle(90.0)
            polyDataMapper.SetInputConnection(polyDataNormals.GetOutputPort())
        else:
            # Migrate to VTK6:
            # http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
            # Old: polyDataMapper.SetInput(polyData)
            polyDataMapper.SetInputData(polyData)

        actor = vtk.vtkActor()
        actor.GetProperty().SetPointSize(3)
        actor.SetMapper(polyDataMapper)
        self.renderer.AddActor(actor)

    def AddFile(self, file_name, colorMap=None):
        file_name_lower = file_name.lower()
        if file_name_lower.endswith('.vtk'):
            polyData = VTKViewer.ReadLegacyVTK(file_name)
        elif file_name_lower.endswith(".vtp"):
            polyData = VTKViewer.readPolyData(
                file_name, vtk.vtkXMLPolyDataReader)
        elif file_name_lower.endswith(".ply"):
            polyData = VTKViewer.readPolyData(
                file_name, vtk.vtkPLYReader)
        elif file_name_lower.endswith(".obj"):
            polyData = VTKViewer.readPolyData(
                file_name, vtk.vtkOBJReader)
        elif file_name_lower.endswith(".stl"):
            polyData = VTKViewer.readPolyData(
                file_name, vtk.vtkSTLReader)
        elif file_name_lower.endswith(".vtu"):
            polyData = VTKViewer.readDataSet(
                file_name, vtk.vtkXMLUnstructuredGridReader)
        elif file_name_lower.endswith(".pdb"):
            polyData = VTKViewer.ReadPDB(file_name)
        elif file_name_lower.endswith(".vti"):
            polyData = VTKViewer.readDataSet(
                file_name, vtk.vtkXMLImageDataReader)
        elif file_name_lower.endswith(".vts"):
            polyData = VTKViewer.readDataSet(
                file_name, vtk.vtkXMLStructuredGridReader)
        elif file_name_lower.endswith(".vtr"):
            polyData = VTKViewer.readDataSet(
                file_name, vtk.vtkXMLRectilinearGridReader)
        else:
            print('{0}: BAD FILE NAME.  Should end in VTK, VTP, PLY, OBJ, '
                  'STL, VTU, or PDB.'.format(file_name))
            raise Exception()
        self.AddPolyData(polyData, colorMap)
        return

    @staticmethod
    def ReadPDB(file_name):
        pdb = vtk.vtkPDBReader()
        pdb.SetFileName(file_name)
        pdb.SetHBScale(1.0)
        pdb.SetBScale(1.0)
        pdb.Update()

        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(0, 0, 0)
        sphere.SetRadius(1)

        glyph = vtk.vtkGlyph3D()
        glyph.SetInputConnection(pdb.GetOutputPort())
        glyph.SetSourceConnection(sphere.GetOutputPort())
        glyph.SetOrient(1)
        glyph.SetColorMode(1)
        glyph.SetScaleMode(2)
        glyph.SetScaleFactor(.25)
        glyph.Update()

        tube = vtk.vtkTubeFilter()
        tube.SetInputConnection(pdb.GetOutputPort())
        tube.SetNumberOfSides(6)
        tube.CappingOff()
        tube.SetRadius(0.2)
        tube.SetVaryRadius(0)
        tube.SetRadiusFactor(10)
        tube.Update()

        tubeMesh = vtk.vtkPolyData()
        tubeMesh.ShallowCopy(tube.GetOutput())
        N = tubeMesh.GetNumberOfPoints()

        rgb_colors = tubeMesh.GetPointData().GetArray("rgb_colors")
        if rgb_colors is not None:
            if rgb_colors.GetNumberOfComponents() == 3:
                for i in range(N):
                    rgb_colors.SetTupleValue(i, (127, 127, 127))

        appendFilter = vtk.vtkAppendPolyData()
        appendFilter.AddInputConnection(glyph.GetOutputPort())
        try:
            appendFilter.AddInputData(tubeMesh)
        except:
            appendFilter.AddInput(tubeMesh)
        appendFilter.Update()

        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(appendFilter.GetOutput())
        return polyData

    @staticmethod
    def ConvertDataSetToSurface(algorithmOutputPort):
        dataSetSurfaceFilter = vtk.vtkDataSetSurfaceFilter()
        dataSetSurfaceFilter.SetInputConnection(algorithmOutputPort)
        dataSetSurfaceFilter.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(dataSetSurfaceFilter.GetOutput())
        return polyData

    @staticmethod
    def readPolyData(file_name, readerType):
        reader = readerType()
        reader.SetFileName(file_name)
        reader.Update()
        polyData = vtk.vtkPolyData()
        polyData.ShallowCopy(reader.GetOutput())
        return polyData

    @staticmethod
    def readDataSet(file_name, readerType):
        reader = readerType()
        reader.SetFileName(file_name)
        reader.Update()
        return VTKViewer.ConvertDataSetToSurface(
            reader.GetOutputPort())

    @staticmethod
    def ReadLegacyVTK(file_name):
        reader = vtk.vtkDataSetReader()
        reader.SetFileName(file_name)
        reader.Update()
        if None != reader.GetPolyDataOutput():
            polyData = vtk.vtkPolyData()
            polyData.ShallowCopy(reader.GetPolyDataOutput())
            return polyData
        if None != reader.GetUnstructuredGridOutput():
            return VTKViewer.ConvertDataSetToSurface(reader.GetOutputPort())
        if None != reader.GetStructuredPointsOutput():
            return VTKViewer.ConvertDataSetToSurface(reader.GetOutputPort())
        if None != reader.GetStructuredGridOutput():
            return VTKViewer.ConvertDataSetToSurface(reader.GetOutputPort())
        if None != reader.GetRectilinearGridOutput():
            return VTKViewer.ConvertDataSetToSurface(reader.GetOutputPort())
        else:
            raise Exception("unsupported: ????????\n")

    @staticmethod
    def GetVTKStereoType(stereoType):
        if (stereoType == "CRYSTAL_EYES"):
            return (vtk.VTK_STEREO_CRYSTAL_EYES)
        elif (stereoType == "RED_BLUE"):
            return (vtk.VTK_STEREO_RED_BLUE)
        elif (stereoType == "INTERLACED"):
            return (vtk.VTK_STEREO_INTERLACED)
        elif (stereoType == "LEFT"):
            return (vtk.VTK_STEREO_LEFT)
        elif (stereoType == "RIGHT"):
            return (vtk.VTK_STEREO_RIGHT)
        elif (stereoType == "DRESDEN"):
            return (vtk.VTK_STEREO_DRESDEN)
        elif (stereoType == "ANAGLYPH"):
            return (vtk.VTK_STEREO_ANAGLYPH)
        elif (stereoType == "CHECKERBOARD"):
            return (vtk.VTK_STEREO_CHECKERBOARD)
        elif (stereoType == "SPLITVIEWPORT_HORIZONTAL"):
            return (vtk.VTK_STEREO_SPLITVIEWPORT_HORIZONTAL)
        else:
            return None


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("FILES",
                        help=('Example: "%(prog)s FILE1 FILE2 FILE3" -- '
                              '"FILE1",... are the names of VTK files.'),
                        nargs='+') #, metavar='')
    parser.add_argument("--use_colormap", action='store_true',
                        help="Use the Paraview-style xml colormap file "
                             "whose location is set by "
                             "the environment variable $COLORMAP.")
    args = parser.parse_args()
    fileNames = args.FILES
    use_colormap = args.use_colormap

    usage = """
    Usage: vtkviewer.py FILE [MORE FILES...] [--use_colormap]
    Supported File Formats:
      *.vtk - VTK Legacy File
      *.vtp - VTK Polygonal Data File
      *.vtu - VTK Unstructured Grid Data File
      *.ply - Stanford Polygon File
      *.obj - Wavefront Object file
      *.stl - Stereolithography File
      *.pdb - Protein Data Bank File
    Controls:
      's' - surface
      'w' - wireframe
      'r' - reset and center camera
      'q' - quit
      '3' - toggle stereo mode
    More Info:
      https://github.com/HalCanary/vtkviewer
    """

    vtkviewer = VTKViewer()

    print(usage)

    if use_colormap and "COLORMAP" in os.environ:
        colormap = VTKViewer.LoadColorMap(os.environ["COLORMAP"])
    else:
        colormap = None

    for fileName in fileNames:
        if os.path.isfile(fileName):
            vtkviewer.AddFile(fileName, colormap)
        else:
            print("Huh?: {0}".format(fileName))
            import sys
            sys.exit()

    vtkviewer.Start()
