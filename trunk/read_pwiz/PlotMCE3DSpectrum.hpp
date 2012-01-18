//
// PlotMCE3DSpectrum.hpp
// 

#ifndef PLOTMCE3DSPECTRUM_HPP
#define PLOTMCE3DSPECTRUM_HPP

#include "pwiz/data/msdata/MSDataFile.hpp"
#include "boost/shared_ptr.hpp"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkDecimatePro.h"
#include "vtkPolyDataNormals.h"
#include "vtkDelaunay2D.h"
#include "vtkLookupTable.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPolyDataMapper.h"
#include "vtkMath.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkShrinkFilter.h"
//#include "vtkDoubleArray.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkPointData.h"

using namespace pwiz::msdata;

//typedef boost::shared_ptr<MSDataList> msdListPtr;
// Rename the smart pointer version of the vtk objects.
typedef vtkSmartPointer<vtkPoints> VtkPointsPtr;
typedef vtkSmartPointer<vtkPolyData> VtkPolyDataPtr;
typedef vtkSmartPointer<vtkDecimatePro> VtkDecimateProPtr;
typedef vtkSmartPointer<vtkSmoothPolyDataFilter> VtkSmoothPolyDataFilterPtr;
typedef vtkSmartPointer<vtkPolyDataNormals> VtkPolyDataNormalsPtr;
typedef vtkSmartPointer<vtkDelaunay2D> VtkDelaunay2DPtr;
typedef vtkSmartPointer<vtkLookupTable> VtkLookupTablePtr;
typedef vtkSmartPointer<vtkUnsignedCharArray> VtkUnsignedCharArrayPtr;
typedef vtkSmartPointer<vtkPolyDataMapper> VtkPolyDataMapperPtr;
typedef vtkSmartPointer<vtkActor> VtkActorPtr;
typedef vtkSmartPointer<vtkRenderer> VtkRendererPtr;
typedef vtkSmartPointer<vtkRenderWindow> VtkRenderWindowPtr;
typedef vtkSmartPointer<vtkRenderWindowInteractor> VtkRenderWindowInteractorPtr;
typedef vtkSmartPointer<vtkPointData> VtkPointDataPtr;

class PlotMCE3DSpectrum 
{
	// The details of the private member information has been hidden!
	// Currently this is a vtkPoints object, but the interface should not care about the exact type.
	public:
		PlotMCE3DSpectrum() {}
		~PlotMCE3DSpectrum() {}
		// PlotMCE3DSpectrumPrivate contains the implementation details.
		VtkPointsPtr passMSData(MSData&);
		VtkPointsPtr passMSDataInPieces(MSData&);
		// displayVtkData() will return the status of the program.
		int displayVtkData(VtkPointsPtr);

};


#endif
