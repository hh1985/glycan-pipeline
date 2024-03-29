/********************************************************************
	created:	2012/11/16
	created:	16:11:2012   16:59
	filename: 	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\test\EnvelopFinderTest.cpp
	file path:	E:\Sync\Dropbox\Codes\glycan_pipeline\GAG\test
	file base:	SimpleFinderTest
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "GAGPL/SPECTRUM/SimpleFinder.h"
#include "GAGPL/SPECTRUM/PeakList.h"
#include <GAGPL/IO/SequenceReader.h>
#include <GAGPL/MISC/Param.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <time.h>

using namespace gag;
using namespace std;
using namespace param;
using namespace gag::io;

namespace po = boost::program_options;

RichList readSpectrum(const std::string& filename)
{
	std::cout << "Load spectrum..." << std::endl;

	std::ifstream infile(filename.c_str());

	std::string line;

	RichList spec;
	if(infile.is_open())
	{
		while(std::getline(infile, line))
		{
			std::istringstream is;
			is.str(line);
			double k1, k2, k3, k4; 
			// The four columns are: m/z, intensity, resolution, s/n
			is >> k1 >> k2 >> k3 >> k4;
			//gag::Peak pk(k1, k2);
			RichPeakPtr pk = createRichPeak(k1, k2, k3, k4);
			spec.addPeak(pk);
		}
	}
	infile.close();
	infile.clear();
	return spec;

}

// argv[1] -- filename without suffix.
// argv[2] -- precursor m/z
// argv[3] -- precursor charge

int main(int argc, char *argv[])
{
	try {
    po::options_description desc("Allowed options");

    string spec_file, out_file;
    double pre_mz;
    int pre_z;
    double error_window;

    desc.add_options()
      ("help,h", "Use -h or --help to list all arguments")
      ("spectrum,s", po::value<string>(&spec_file), "Spectrum peak list")
      ("output,o", po::value<string>(&out_file), "File for monoisotopic peak list")
      //("library,l", po::value<string>(&lib_file), "File for Fragment library")
      ("error,e", po::value<double>(&error_window)->default_value(2.0), "Matching error (ppm)")
      ("precursor,p", po::value<double>(&pre_mz), "Precursor m/z.")
      ("charge,z", po::value<int>(&pre_z), "Charge state.");

    po::variables_map vm;

    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if(argc == 1 || vm.count("help")) {
      cout << desc << endl;
      return 0;
    }

    if(!vm.count("spectrum")) {
      cerr << "Option \"spectrum\"(-s) is required!" << endl;
      return 1;
    }

    if(!vm.count("precursor")) {
      cerr << "Precursor m/z is requried!"<< endl;
      return 1;
    }

    if(!vm.count("charge")) {
      cerr << "Charge state is required!" << endl;
      return 1;
    }

    // Parse the mono_file string and generate library file.
    boost::filesystem::path p(spec_file);

    if(!vm.count("output")) {
      boost::filesystem::path assign_path = p.parent_path() / (p.stem().string() + "_mono.txt");
      out_file = assign_path.string();
    }

    boost::filesystem::path assign_path(out_file);
    auto alter_path = assign_path.parent_path() / (assign_path.stem().string() + "_undetermined.txt");
    string alter_file = alter_path.string();

    clock_t t = clock();

    // Modify the parameter in xml.
    Param& param = Param::Instance();
    param.setParameter<double>("matching_error", error_window);
    param.setParameter<double>("precursor_mz", pre_mz);
    param.setParameter<int>("precursor_charge", pre_z);
    param.setParameter<std::string>("instrument_type", "ft");

    // threshold for observed peak.
    param.setParameter<double>("lower_bound_sn",5.0);
    param.setParameter<double>("learning_rate", 0.01);
    param.setParameter<double>("max_intensity_shift", 0.3);
    param.setParameter<int>("max_sulfur_num",6);
    param.setParameter<int>("num_missing_peaks",1);


    RichList spectrum = readSpectrum(spec_file);

    SimpleFinder env_finder2(spectrum);

    std::multimap<RichPeakPtr, EnvelopPtr> env_pool = env_finder2.getEnvelops();

    std::ofstream outfile1(out_file.c_str());

    if(outfile1.is_open()) {
      outfile1 << "MZ\tIntensity\tCharge\n"; 
      outfile1.precision(5);

      for(std::multimap<RichPeakPtr, EnvelopPtr>::iterator iter = env_pool.begin(); iter != env_pool.end(); iter++)
      {
        env_finder2.printEnvelop(iter->second);

        outfile1 << fixed << iter->first->mz << "\t" << iter->first->intensity << "\t" << iter->second->charge_state << "\n"; 		
      }

      outfile1.close();
    } else {
      throw std::runtime_error("Unable to open the file: " + out_file);
    }

    std::ofstream outfile2(alter_file.c_str());

    if(outfile2.is_open()) {
      outfile2 << "MZ\tIntensity\tCharge\n";
      outfile2.precision(5);

      std::multimap<RichPeakPtr, int>& pks = env_finder2.getUndeterminedPeaks();
      for(std::multimap<RichPeakPtr, int>::iterator iter = pks.begin(); iter!= pks.end(); iter++) 
      {
        if(iter->first->pk_type == "NOISE")
          continue;
        outfile2 << fixed << iter->first->mz << "\t" << iter->first->intensity << "\t" << iter->second << "\n";
      }
      outfile2.close();

    } else {
      throw std::runtime_error("Unable to open the file: " + alter_file);
    }

    cout << "Monoisotopic peak list was stored in " << out_file << endl;
    cout << "Undetermined monoisotopic peak list was stored in " << alter_file << endl;
    cout << "The match error is " << error_window << " ppm" << endl;

    t = clock() - t;

		/*std::string filename = "./data/spectrum/";
		filename.append(argv[1]);
		filename.append(".txt");

		RichList spectrum = readSpectrum(filename);

		Param& param = Param::Instance();
		param.setParameter<double>("precursor_mz", atof(argv[2]));
		param.setParameter<int>("precursor_charge", atoi(argv[3]));

		 Threshold for most abundant peak in an envelop (mono peak in current case.
		param.setParameter<double>("signal_noise", 15.0);
		param.setParameter<std::string>("instrument_type", "ft");

		 threshold for observed peak.
		param.setParameter<double>("lower_bound_sn",5.0);
		param.setParameter<double>("learning_rate", 0.01);
		param.setParameter<double>("max_intensity_shift", 0.3);
		param.setParameter<int>("max_sulfur_num",6);
		param.setParameter<int>("num_missing_peaks",1);

		SimpleFinder env_finder2(spectrum);

		 2. print the identified envelops.
		std::multimap<RichPeakPtr, EnvelopPtr> env_pool = env_finder2.getEnvelops();



		cout << "Output the results into file." << endl;

		std::string path = "./data/output/";
		path.append(argv[1]);
		path.append("_mono_list.txt");

		ofstream outfile(path.c_str());

		if(outfile.is_open()) {
			outfile << "MZ\tIntensity\tCharge\n"; 
			outfile.precision(5);
			cout.precision(5);

			for(std::multimap<RichPeakPtr, EnvelopPtr>::iterator iter = env_pool.begin(); iter != env_pool.end(); iter++)
			{
				env_finder2.printEnvelop(iter->second);

				outfile << fixed << iter->first->mz << "\t" << iter->first->intensity << "\t" << iter->second->charge_state << "\n"; 		
			}

			outfile.close();
		} else {
			std::cout << "Unable to open file!\n" << std::endl;
		}

		std::string und_pks = "./data/output/";
		und_pks.append(argv[1]);
		und_pks.append("_undetermined_list.txt");
		ofstream outfile2(und_pks.c_str());

		if(outfile2.is_open()) {
			outfile2 << "MZ\tIntensity\tCharge\n";
			outfile2.precision(5);

			std::multimap<RichPeakPtr, int>& pks = env_finder2.getUndeterminedPeaks();
			for(std::multimap<RichPeakPtr, int>::iterator iter = pks.begin(); iter!= pks.end(); iter++) 
			{
				if(iter->first->pk_type == "NOISE")
					continue;
				outfile2 << fixed << iter->first->mz << "\t" << iter->first->intensity << "\t" << iter->second << "\n";
			}
			outfile2.close();

		} else {
			std::cout << "Unable to open file " << und_pks << std::endl;
		}*/
	} catch(std::exception& e) {
		cout << "Exception: " << e.what() << endl;
	}

	return 0;

}