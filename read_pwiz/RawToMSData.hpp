//
// RawToMSData.hpp
//

#ifndef _RAWTOMSDATA_HPP
#define _RAWTOMSDATA_HPP

#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/Std.hpp"

using namespace pwiz::msdata;

class RawToMSData {
//	private:
//		vector<MSDataPtr> msdList;
//		//vector<MSDataPtr> msdList;
	public:	
		void ReadRaw(const string& filename);
		void DataManipulation(MSData&);
		// The role of read raw should be to fill the MSData object

};
#endif // RAWTOMSDATA_INCLUDED
