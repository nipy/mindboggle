# Directory containing class headers.
SET(VTK_HYBRID_HEADER_DIR "${VTK_INSTALL_PREFIX}/include/vtk-5.6")

# Classes in vtkHybrid.
SET(VTK_HYBRID_CLASSES
  "vtk3DSImporter"
  "vtkArcPlotter"
  "vtkAnnotatedCubeActor"
  "vtkAxisActor"
  "vtkAxesActor"
  "vtkBarChartActor"
  "vtkCaptionActor2D"
  "vtkCornerAnnotation"
  "vtkCubeAxesActor"
  "vtkCubeAxesActor2D"
  "vtkDepthSortPolyData"
  "vtkEarthSource"
  "vtkFacetReader"
  "vtkGreedyTerrainDecimation"
  "vtkGridTransform"
  "vtkImageDataLIC2D"
  "vtkImageDataLIC2DExtentTranslator"
  "vtkImageToPolyDataFilter"
  "vtkImplicitModeller"
  "vtkIterativeClosestPointTransform"
  "vtkLandmarkTransform"
  "vtkLegendBoxActor"
  "vtkLegendScaleActor"
  "vtkLSDynaReader"
  "vtkPCAAnalysisFilter"
  "vtkPieChartActor"
  "vtkPolyDataSilhouette"
  "vtkPolyDataToImageStencil"
  "vtkProcrustesAlignmentFilter"
  "vtkProjectedTerrainPath"
  "vtkRIBExporter"
  "vtkRIBLight"
  "vtkRIBProperty"
  "vtkRenderLargeImage"
  "vtkSpiderPlotActor"
  "vtkStructuredExtent"
  "vtkStructuredGridLIC2D"
  "vtkTemporalDataSetCache"
  "vtkTemporalInterpolator"
  "vtkTemporalShiftScale"
  "vtkTemporalSnapToTimeStep"
  "vtkThinPlateSplineTransform"
  "vtkTransformToGrid"
  "vtkVRMLImporter"
  "vtkVectorText"
  "vtkVideoSource"
  "vtkWeightedTransformFilter"
  "vtkXYPlotActor"
  "vtkX3D"
  "vtkX3DExporter"
  "vtkX3DExporterWriter"
  "vtkX3DExporterXMLWriter"
  "vtkX3DExporterFIWriter"
  "vtkExodusIICache"
  "vtkExodusIIReader"
  "vtkExodusIIReaderParser"
  "vtkExodusIIReaderVariableCheck"
  "vtkExodusReader"
  "vtkDSPFilterDefinition"
  "vtkExodusModel"
  "vtkDSPFilterGroup")

# Abstract classes in vtkHybrid.
SET(VTK_HYBRID_CLASSES_ABSTRACT)

# Wrap-exclude classes in vtkHybrid.
SET(VTK_HYBRID_CLASSES_WRAP_EXCLUDE
  "vtkStructuredGridLIC2D"
  "vtkX3D"
  "vtkX3DExporterWriter"
  "vtkX3DExporterXMLWriter"
  "vtkX3DExporterFIWriter"
  "vtkExodusIIReaderParser"
  "vtkExodusIIReaderVariableCheck")

# Set convenient variables to test each class.
FOREACH(class ${VTK_HYBRID_CLASSES})
  SET(VTK_CLASS_EXISTS_${class} 1)
ENDFOREACH(class)
FOREACH(class ${VTK_HYBRID_CLASSES_ABSTRACT})
  SET(VTK_CLASS_ABSTRACT_${class} 1)
ENDFOREACH(class)
FOREACH(class ${VTK_HYBRID_CLASSES_WRAP_EXCLUDE})
  SET(VTK_CLASS_WRAP_EXCLUDE_${class} 1)
ENDFOREACH(class)
