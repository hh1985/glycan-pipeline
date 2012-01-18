#include "read_pwiz/PlotMCE3DSpectrum.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include <vector>
using namespace pwiz::msdata;

void testPlot(MSData& msd)
{
		PlotMCE3DSpectrum msdtovtk;
		// Initialization
		//VtkPointsPtr points = msdtovtk.passMSData(msd);
		VtkPointsPtr points = msdtovtk.passMSDataInPieces(msd);
		// Plot
		msdtovtk.displayVtkData(points);
	
}

void testRead(const string& filename)
{
	FullReaderList readers;
	
	vector<MSDataPtr> msdlist; 

	readers.read(filename, msdlist);

	for (size_t i=0; i < msdlist.size(); i++)
	{
		MSData& msd = *msdlist[i];
		testPlot(msd);
	}
	
}



int main(int argc, char* argv[])
{
	try
	{
		if (argc != 2)
			throw runtime_error("Usage: PlotMCE3DSpectrumTest.cpp infile outfile\n");
		testRead(argv[1]);
		return 0;

	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
	}
	catch (...)
	{
		cerr << "Caught unknown exception. \n";
	}
	return 1;
}

