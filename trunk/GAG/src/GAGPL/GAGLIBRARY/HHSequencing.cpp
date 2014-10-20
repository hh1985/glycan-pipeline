/********************************************************************
	created:	2014/01/15
	created:	15:1:2014   13:24
	filename: 	HHSequencing.cpp
	file path:	GAGPL\GAGLIBRARY
	file base:	HHSequencing
	file ext:	cpp
	author:		Han Hu
	
	purpose:	The main function for HS-SEQ algorithm
*********************************************************************/
//#include <iostream>
#include <GAGPL/GAGLIBRARY/SequencePrediction.h>
#include <GAGPL/IO/SequenceReader.h>
#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <GAGPL/MISC/Param.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <time.h>

using namespace std;
using namespace gag;
using namespace gag::io;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  try
  {
    po::options_description desc("Allowed options");

    string compo_file, mono_file, out_file;
    double error_window;
    desc.add_options()
      ("help,h", "Use -h or --help to list all arguments")
      ("composition,c", po::value<string>(&compo_file), "HS sequence composition and derivation")
      ("mono-peak-list,m", po::value<string>(&mono_file), "Monoisotopic peak list: m/z, intensity, z")
      ("output,o", po::value<string>(&out_file), "File for peak assignments")
      //("library,l", po::value<string>(&lib_file), "File for Fragment library")
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
      boost::filesystem::path assign_path = p.parent_path() / (p.stem().string() + "_dist.txt");
      out_file = assign_path.string();
    }

    // Modify the parameter in xml.
    Param& param = Param::Instance();
    param.setParameter<double>("matching_error", error_window);

    clock_t t = clock();

    SequenceReader seq_reader;

    static GlycanSequencePtr gs = boost::make_shared<GlycanSequence>();
    seq_reader.readSequenceStructure(compo_file, gs);

    set<MonoPeakPtr> mono_set;
    seq_reader.readMonoPeakList(mono_file, mono_set);

    FragmentTree glyco_tree(gs);
    cout << "Fragment library size: " << glyco_tree.getNodeSize() << "\n";

    std::multimap<MonoPeakPtr, NodeItem> node_map = glyco_tree.searchLibrary(mono_set);

    cout << "Assignment data size: " << node_map.size() << "\n";

    // Predict the structure from the backbone and mass list.
    SequencePrediction seq_predictor(node_map, gs);

    set<string> mod_types = gs->getModificationTypes();

    std::ofstream outfile(out_file.c_str());

    if(outfile.is_open()) {
      outfile << "Residue\tPosition(ResidueID - SiteID)\tModification\tIntensity\n";
      for(auto iter = mod_types.begin(); iter != mod_types.end(); iter++)
      {
        ModificationDistribution mod_dist = seq_predictor.getModificationDistribution(*iter);
        for(ModificationDistribution::const_iterator mod_iter = mod_dist.begin(); mod_iter != mod_dist.end(); mod_iter++)
        {
          outfile << gs->getBranchByID(0).getUnitByID(mod_iter->first.getMonosaccharideID()).getSymbol() << "\t";
          outfile << mod_iter->first.printString() << "\t" << *iter << "\t" << mod_iter->second << "\n";
        }

      }
      outfile.close();

      cout << "Results were stored in " << out_file << endl;
      cout << "The match error is " << error_window << " ppm" << endl;

      t = clock() - t;

      cout << "Elapsed time: " << (double)t/CLOCKS_PER_SEC << "s\n";
    } else {
      throw std::runtime_error("Unable to open the file: " + out_file);
    }

    
    //cin.get();
    return 0;
  }
  catch(std::exception& e)
  {
    cout << "Exception: " << e.what() << endl;
  }
}