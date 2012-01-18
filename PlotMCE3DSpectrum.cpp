//
// PlotMCE3DSpectrum.cpp
//
#include "read_pwiz/PlotMCE3DSpectrum.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/Std.hpp"
// Not quite sure if this is the solution. But at least give it a try.
//#include "pwiz/analysis/passsive/MSDataAnalyzer.hpp"
//#include "pwiz/utility/misc/unit.hpp"
#include "pwiz/analysis/passive/MSDataCache.hpp"
//#include "PeakExtractor.hpp"
#include <vector>

using namespace pwiz::msdata;
using namespace pwiz::analysis;

// The MSData object will be converted into vtkPoints object, or even other structure:)
VtkPointsPtr PlotMCE3DSpectrum::passMSData(MSData& msd)
{
	SpectrumListPtr sl = msd.run.spectrumListPtr;
	
	if(!sl.get())
		throw runtime_error("Correct me!");
	else
		cout << "# of spectra: " << sl->size() << endl;
	
	VtkPointsPtr points = VtkPointsPtr::New();

	// Iterate over each spectrum.
	for (size_t i=0, size=sl->size(); i!=size; i++)
	{
		const bool getBinaryData = true;
		cout << "Successful entered: " << (int) i << endl;
		SpectrumPtr sp = sl->spectrum(i, getBinaryData);
		cout << "Successful get the spectrum: " << (int) i << endl;

		if(sp->precursors.empty())
			continue;

		vector<MZIntensityPair> pairs;
		sp->getMZIntensityPairs(pairs);

		double retention_time = sp->scanList.scans[0].cvParam(MS_scan_start_time).timeInSeconds();
		
		for (vector<MZIntensityPair>::const_iterator it=pairs.begin(), end=pairs.end(); it!=end; ++it)
		{
			// Different ways to get the subset of MZIntensityPair
			// Reduce useless points.
			if(it->intensity < .01) 
				continue;
			points->InsertNextPoint(retention_time, it->mz, it->intensity);
		}
	}
	return points;

}

// PassData in pieces.
// The input parameters include the amount of pieces to read and the msd object.
// Each time, the pieces will be stored in the memory, visualized and appeneded to previous results.
VtkPointsPtr PlotMCE3DSpectrum::passMSDataInPieces(MSData &msd)
{
	MSDataCache msdCache;
	msdCache.open(msd);
	
	const size_t spectrum_count = msdCache.size();
	
	VtkPointsPtr points = VtkPointsPtr::New();
	
	for (size_t index=0; index<spectrum_count; index++)
	{
		const SpectrumInfo& spectrum_info = msdCache.spectrumInfo(index, true);
		vector<MZIntensityPair> pairs = spectrum_info.data;

		for (vector<MZIntensityPair>::const_iterator it=pairs.begin(), end=pairs.end(); it!=end; ++it)
		{
			if (it->intensity < .01)
				continue;
			points->InsertNextPoint(spectrum_info.retentionTime, it->mz, it->intensity);
		}
	}

	return points;
}


int PlotMCE3DSpectrum::displayVtkData(VtkPointsPtr points)
{
	// Check if points is a null object
	//
	if (!points)
		throw runtime_error("vtkPoints object has not been assigned any value!\n");

	VtkPolyDataPtr inputpolydata = VtkPolyDataPtr::New();
	
	inputpolydata->SetPoints(points);
	
	VtkDecimateProPtr deci = VtkDecimateProPtr::New();
	deci->SetInput(inputpolydata);
	deci->SetTargetReduction(0.9);
	deci->PreserveTopologyOn();

	VtkSmoothPolyDataFilterPtr smoother = VtkSmoothPolyDataFilterPtr::New();
	smoother->SetInput(deci->GetOutput());
	smoother->SetNumberOfIterations(50);

	VtkPolyDataNormalsPtr normals = VtkPolyDataNormalsPtr::New();
	normals->SetInput(smoother->GetOutput());
	normals->FlipNormalsOn();

	VtkDelaunay2DPtr delaunay = VtkDelaunay2DPtr::New();
	delaunay->SetInput(inputpolydata);
	delaunay->Update();

	VtkPolyDataPtr outputpolydata = delaunay->GetOutput();

	double bounds[6];
	outputpolydata->GetBounds(bounds);

	double minz = bounds[4];
	double maxz = bounds[5];

	VtkLookupTablePtr colorlookuptable = VtkLookupTablePtr::New();
	colorlookuptable->SetTableRange(minz, maxz);
	colorlookuptable->Build();

	VtkUnsignedCharArrayPtr colors = VtkUnsignedCharArrayPtr::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("colors");
	
	cout << "# of PolyData points: " << outputpolydata->GetNumberOfPoints() << endl;

	for(int i=0; i < outputpolydata->GetNumberOfPoints(); i++) {
		double p[3];
		outputpolydata->GetPoint(i,p);
	
		double dcolor[3];
		colorlookuptable->GetColor(p[2], dcolor);
		unsigned char color[3];
		for(unsigned int j=0; j<3; j++)
		{
			color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
		}
		colors->InsertNextTupleValue(color);
	}
	
	outputpolydata->GetPointData()->SetScalars(colors);
	
	VtkPolyDataMapperPtr mapper = VtkPolyDataMapperPtr::New();
	mapper->SetInputConnection(outputpolydata->GetProducerPort());
	
	VtkActorPtr actor = VtkActorPtr::New();
	actor->SetMapper(mapper);
	
	VtkRendererPtr renderer = VtkRendererPtr::New();
	VtkRenderWindowPtr renderwindow = VtkRenderWindowPtr::New();
	
	renderwindow->AddRenderer(renderer);
	VtkRenderWindowInteractorPtr renderwindowinteractor = VtkRenderWindowInteractorPtr::New();
	renderwindowinteractor->SetRenderWindow(renderwindow);
	
	renderer->AddActor(actor);
	renderer->SetBackground(.1, .2, .3);
	
	renderwindow->Render();
	renderwindowinteractor->Start();
	
	return EXIT_SUCCESS;

}


