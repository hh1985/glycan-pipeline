/********************************************************************
	created:	2014/10/19
	created:	19:10:2014   19:52
	filename: 	LibraryMatchMain.cpp
	file path:	GAG\src\GAGPL\GAGLIBRARY
	file base:	LibraryMatchMain
	file ext:	cpp
	author:		Han Hu
	
	purpose:	Generation of fragment library and peak matching.
*********************************************************************/
#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <GAGPL/MATH/MassConversion.h>
#include <GAGPL/IO/SequenceReader.h>
#include "GAGPL/MISC/Param.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <time.h>

using namespace std;
using namespace gag;
using namespace msmath;
using namespace gag::io;

namespace po = boost::program_options;

void exportAssignmentList(const std::multimap<MonoPeakPtr, NodeItem>& node_map, const std::string& filename)
{
  std::ofstream outfile(filename.c_str());
  outfile.precision(5);

  if(outfile.is_open()) {
    // print title.
    outfile << "MZ\tCharge\tEXP_Mass\tTheo_Mass\tError\tCleavage_Type\tCleavage_Num\tComposition_Shift\tComposition\tAc\tSO3\n";
    for(std::multimap<MonoPeakPtr, NodeItem>::const_iterator iter = node_map.begin(); iter != node_map.end(); iter++)
    {
      double mass = msmath::calculateMass(iter->first->mz, -1 * iter->first->z);
      double ppm = (iter->second.getMass() - mass)/mass * 1e6;

      outfile << std::fixed << iter->first->mz << "\t" << iter->first->z << "\t" << mass << "\t" << iter->second.getMass() << "\t" << ppm << "\t" << iter->second.getCleavageType() << "\t" << iter->second.getCleavageNum() << "\t" << iter->second.getCompositionShift() << "\t" << iter->second.getCompositionString() << "\t" << iter->second.getModificationNum("Ac") << "\t" << iter->second.getModificationNum("SO3")  << "\n";
    }
    outfile.close();
  } else {
    throw runtime_error("Unable to open file: " + filename);
  }
}

int main(int argc, char *argv[]) 
{
  try 
  {
    po::options_description desc("Allowed options");
    
    string compo_file, mono_file, out_file, lib_file;
    double error_window;
    desc.add_options()
      ("help,h", "Use -h or --help to list all arguments")
      ("composition,c", po::value<string>(&compo_file), "HS sequence composition and derivation")
      ("mono-peak-list,m", po::value<string>(&mono_file), "Monoisotopic peak list: m/z, intensity, z")
      ("output,o", po::value<string>(&out_file), "File for peak assignments")
      ("library,l", po::value<string>(&lib_file), "File for Fragment library")
      ("error,e", po::value<double>(&error_window)->default_value(2.0), "Matching error (ppm)");

    po::variables_map vm;

    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if(argc == 1 || vm.count("help")) {
      cout << desc << endl;
      return 0;
    }

    if(!vm.count("composition") || !vm.count("mono-peak-list")) {
      cerr << "Option \"composition\"(-c) and \"mono-peak-list\"(-m) are required!" << endl;
      return 1;
    }

    // Parse the mono_file string and generate library file.
    boost::filesystem::path p(mono_file);
    

    if(!vm.count("output")) {
      boost::filesystem::path assign_path = p.parent_path() / (p.stem().string() + "_assignments.txt");
      out_file = assign_path.string();
    }

    if(!vm.count("library")) {
      boost::filesystem::path lib_path = p.parent_path() / (p.stem().string() + "_library.txt");
      lib_file = lib_path.string();
    }

    cout << "Fragment library will be stored in " << lib_file << endl;
    cout << "Assignment list will be stored in " << out_file << endl;
    cout << "The match error is " << error_window << " ppm" << endl;
 
    // Modify the parameter in xml.
    Param& param = Param::Instance();
    param.setParameter<double>("matching_error", error_window);

    clock_t t = clock();

    set<MonoPeakPtr> peak_list;
    loadMonoPeakList(mono_file, peak_list);

    SequenceReader seq_reader;
    GlycanSequencePtr glyco_seq = boost::make_shared<GlycanSequence>();
    seq_reader.readSequenceStructure(compo_file, glyco_seq);

    FragmentTree glyco_tree(glyco_seq);
    cout << "Fragment library size: " << glyco_tree.getNodeSize() << "\n";

    glyco_tree.exportLibrary(lib_file);

    std::multimap<MonoPeakPtr, NodeItem> node_map = glyco_tree.searchLibrary(peak_list);
    cout << "Assignment data size: " << node_map.size() << endl; 

    exportAssignmentList(node_map, out_file);

    t = clock() - t;
    cout << "Elapsed time: " << ((double)t)/CLOCKS_PER_SEC << endl;

    return 0;
  }
  catch(std::exception& e) {
    cout << "Exception: " << e.what() << endl;
  }
  
}
