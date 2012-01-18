#include "read_pwiz/RawToMSData.hpp"

using namespace std;


int main(int argc, char* argv[])
{
	try
	{
		if (argc != 2)
			throw runtime_error("Usage: RawToMSDataTest.cpp infile outfile");
		RawToMSData D;
		D.ReadRaw(argv[1]);
		return 0;
	}
	catch (exception& e)
	{
		cerr << e.what() << endl;
	}
	catch(...)
	{
		cerr << "Caught unknown exception.\n";
	}
	return 1;
}
