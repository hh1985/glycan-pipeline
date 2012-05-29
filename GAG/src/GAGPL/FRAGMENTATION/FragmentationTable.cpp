/*
 * =====================================================================================
 *
 *       Filename:  FragmentationTable.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/29/2012 10:16:24 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/FRAGMENTATION/FragmentationTable.h>

namespace gag
{
	FragmentationTable& FragmentationTable::Instance()
	{
		static FragmentationTable ftt;
		return ftt;
	}

	void FragmentationTable::load(const std::string& filename)
	{
		//FragmentationTable& ftable = FragmentationTable::Instance();
		//ftable.load();

		using boost::property_tree::ptree;
		ptree pt;

		read_xml(filename, pt);

		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.CleavageTypes"))
		{
			FragmentationParams fp;
			fp.type = v.second.get<std::string>("Name");
			BOOST_FOREACH(ptree::value_type &s, v.second.get_child("Shift"))
			{
				if(s.first == "Value")
					fp.cleavage_shift = s.second.data();
				else if(s.first == "Dissociation") {
					std::string name = s.second.get<std::string>("Name");
					std::vector<std::string> shift_vals;
					std::pair<ptree::assoc_iterator, ptree::assoc_iterator> ret = s.second.equal_range("Value");
					for(ptree::assoc_iterator iter = ret.first; iter != ret.second; iter++)
					{
						shift_vals.push_back(iter->second.data());
					}
					fp.dis_shift.insert(std::make_pair(name, shift_vals));
				}
			}


			//if(v.first == "Name")
			//	fp.type = v.second.data();
			//else if(v.first == "Shift")	{
			//	fp.cleavage_shift = v.second.get<std::string>("Value");
			//	std::vector<std::string> shift_vals;
			//	std::string name;
			//	BOOST_FOREACH(ptree::value_type &s, v.second.get_child("Dissociation"))
			//	{
			//		if(s.first == "Value")
			//			shift_vals.push_back(s.second.data());
			//		if(s.first == "Name")
			//			name = s.second.data();
			//	}
			//	fp.dis_shift.insert(std::make_pair(name, shift_vals));
			//}
			fragmentation_params.insert(std::make_pair(fp.type, fp));	
			
		}

		pt.erase("parameters.CleavageTypes");

	}

	FragmentationParams FragmentationTable::getFragmentationParams(const std::string& cleavage_type) const
	{
		std::map<std::string, FragmentationParams>::const_iterator i = fragmentation_params.find(cleavage_type);
		return i != fragmentation_params.end()? i->second : FragmentationParams();
	}

	CompositionShift FragmentationTable::getCompositionShift(const std::string& str) const
	{
		if(str.empty())
			return CompositionShift();
		CompositionShift compo_shift;
		std::string::const_iterator it = str.begin();
		if((*it) == '-') { // Negative
			compo_shift.first = -1;
			compo_shift.second = str.substr(1);
		}	else {
			compo_shift.first = 1;
			compo_shift.second= str;
		}
		return compo_shift;
	}
	
	CompositionShift FragmentationTable::getCleavageShift(const std::string& cleavage_type) const
	{
		std::map<std::string, FragmentationParams>::const_iterator i = fragmentation_params.find(cleavage_type);

		return i!=fragmentation_params.end() ? (*this).getCompositionShift(i->second.cleavage_shift) : CompositionShift();
	}

	std::vector<CompositionShift> FragmentationTable::getCleavageShift(const std::string& cleavage_type, const std::string& dis_type) const
	{
		std::map<std::string, FragmentationParams>::const_iterator i = fragmentation_params.find(cleavage_type);
		
		std::vector<CompositionShift> compo_shift_vec;

		if(i != fragmentation_params.end()) {
			std::map<std::string, std::vector<std::string> >::const_iterator j = i->second.dis_shift.find(dis_type);
			
			if(j != i->second.dis_shift.end())
				
				for(std::vector<std::string>::const_iterator iter = j->second.begin(); iter != j->second.end(); iter++)
				{
					compo_shift_vec.push_back((*this).getCompositionShift(*iter));
				}
		} 

		return compo_shift_vec;
	}
}