// 
// RawToMSData.cpp
//

#include "read_pwiz/RawToMSData.hpp"
#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/Std.hpp"

using namespace pwiz::msdata;
using namespace pwiz::cv;
using namespace pwiz::data;

// Might support reading of multiple files.
// Needs modification.

void RawToMSData::ReadRaw(const string& filename)
{
	//convert raw file(s) to MSData object
	
	FullReaderList readers;
	vector<MSDataPtr> msdList;

	readers.read(filename, msdList);

	for (size_t i=0; i < msdList.size(); i++)
	{
		MSData& msd = *msdList[i];

		DataManipulation(msd);

	}


}

void RawToMSData::DataManipulation(MSData& msd)
{
	// Manipulate your MSData object here, or even the rawfile object.
	// Try to get a list of precursor ions.
	
	SpectrumListPtr sl = msd.run.spectrumListPtr;

	if(!sl.get())
		throw runtime_error("[mCE.cpp] No spectra found.");
	else
		cout << "# of spectra: " << sl->size() << endl;

	// write header
	const size_t columnWidth = 14;

	cout << "#" 
	     << setw(columnWidth) << "scanNumber"
	     << setw(columnWidth) << "PreMZ"
	     << setw(columnWidth) << "PreZ"
	     << setw(columnWidth) << "PreInt"
	     << setw(columnWidth) << "CE"
	     << setw(columnWidth) << "msLevel"
	     << setw(columnWidth) << "ProMZ"
	     << setw(columnWidth) << "intensity" << endl;

	// iterate through the spectra in the SpectrumList
	for (size_t i=0, size=sl->size(); i!=size; i++)
	{
		// retrieve the spectrum, with binary data
		const bool getBinaryData = true;
		SpectrumPtr sp = sl->spectrum(i, getBinaryData);

		// Only spectrum with precursor information is retrieved.
		if(sp->precursors.empty())
			continue;
		
		// fill in MZIntensityPair vector for convenient access to binary data
		vector<MZIntensityPair> pairs;
		sp->getMZIntensityPairs(pairs);

		Precursor& precursor = sp->precursors[0];
		
		//cout << sp->precursors.activation.cvParam(MS_collision_energy).value << endl;
		
		// iterate through the m/z-intensity pairs
		for (vector<MZIntensityPair>::const_iterator it=pairs.begin(), end=pairs.end(); it!=end; ++it)
		{
			
			cout << " "
			     << setw(columnWidth) << sp->index
			     << setw(columnWidth) << precursor.selectedIons[0].cvParam(MS_selected_ion_m_z).valueAs<double>()
			     << setw(columnWidth) << precursor.selectedIons[0].cvParam(MS_charge_state).valueAs<int>()
			     << setw(columnWidth) << precursor.selectedIons[0].cvParam(MS_peak_intensity).valueAs<double>()
			     << setw(columnWidth) << precursor.activation.cvParam(MS_collision_energy).valueAs<double>()
			     << setw(columnWidth) << "ms" + sp->cvParam(MS_ms_level).value
			     << setw(columnWidth) << fixed << setprecision(4) << it->mz
			     << setw(columnWidth) << fixed << setprecision(2) << it->intensity << endl;
		}
	}	

	return; 

}


