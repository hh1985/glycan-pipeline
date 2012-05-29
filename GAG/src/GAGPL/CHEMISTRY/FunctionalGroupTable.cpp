/*
 * =====================================================================================
 *
 *       Filename:  FunctionalGroupTable.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/29/2012 10:16:24 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>

namespace gag
{
	FunctionalGroupTable& FunctionalGroupTable::Instance()
	{
		static FunctionalGroupTable cgt;
		return cgt;
	}

	void FunctionalGroupTable::load(const std::string& filename)
	{
		PeriodicTable& ptable = PeriodicTable::Instance();
		ptable.load();

		using boost::property_tree::ptree;
		ptree pt;

		read_xml(filename, pt);

		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.FunctionalGroupSets"))
		{
			if(v.first == "FunctionalGroup")
			{
				FunctionalGroup fg;
				fg.symbol = v.second.get<std::string>("Symbol");
				//std::cout << fg.symbol << std::endl;
				fg.name = v.second.get<std::string>("Name");
				const std::string compo_string = v.second.get<std::string>("Composition");
				fg.getComposition().update(compo_string );
				
				if(v.second.count("Sites")>0)
				{
					BOOST_FOREACH(ptree::value_type &s, v.second.get_child("Sites"))
					{
						// Composition nodes.
						if(s.first == "Composition")
						{
							std::string cp_str = s.second.data();
							Composition temp_cp(cp_str);
							fg.site_gps.push_back(temp_cp);
						}
					}
				}
				functionalgroups.insert(std::make_pair(fg.symbol, fg));
				
			}
		}

		pt.erase("parameters.FunctionalGroupSets");

	}

	FunctionalGroup FunctionalGroupTable::getFunctionalGroupBySymbol(const std::string& symbol) const
	{
		std::map<std::string, FunctionalGroup>::const_iterator i = functionalgroups.find(symbol);

		return i != functionalgroups.end() ? i->second : FunctionalGroup();
	}
}
