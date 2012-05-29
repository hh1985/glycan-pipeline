/*
 * =====================================================================================
 *
 *       Filename:  PeriodicTable.cpp
 *
 *    Description:  The specific implementation for loading periodic table data.  
 *
 *        Version:  1.0
 *        Created:  4/24/2012 3:06:40 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com 
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/CHEMISTRY/PeriodicTable.h>
#include <iostream>

namespace gag
{

	PeriodicTable& PeriodicTable::Instance()
	{
		static PeriodicTable pt;
		return pt;
	}

	void PeriodicTable::load(const std::string& filename)
	{
		using boost::property_tree::ptree;
		ptree pt;
		
		read_xml(filename, pt);

		//ptree es = pt.get_child("parameters.ElementIsotopes");

		//std::pair<ptree::assoc_iterator, ptree::assoc_iterator> pcld = 
		//	es.equal_range("Element");
		BOOST_FOREACH(ptree::value_type &v, pt.get_child("parameters.ElementIsotopes"))
		{
			//std::cout << v.first << std::endl;
			if(v.first == "Element")
			{
				Element elem;
				elem.symbol = v.second.get<std::string>("Symbol");
				elem.name = v.second.get<std::string>("Name");
				elem.atomicity = v.second.get<int>("Atomicity");
				//std::cout << "Before Isotope" << std::endl;
				std::pair<ptree::assoc_iterator, ptree::assoc_iterator> ei = v.second.equal_range("Isotope");
				for(ptree::assoc_iterator iter = ei.first; iter != ei.second; iter++)
				{
					//std::cout << "After Isotope" << std::endl;
					//std::cout << iter->first << std::endl;
					//std::cout << iter->second.get<double>("Mass") << " " << iter->second.get<float>("Probability") << std::endl;
					Isotope iso = {iter->second.get<double>("Mass"), iter->second.get<float>("Probability")};
					//std::cout << iso.mass << " " << iso.abundance << std::endl;
					elem.isotopes.push_back(iso);					
				}

				elements.insert(std::make_pair(elem.symbol, elem));
			}
		}	
		// Delete params and keep all information only in elements.
		pt.erase("parameters.ElementIsotopes");
	
	}


	const Element PeriodicTable::getElementBySymbol(const std::string& symbol) const
	{
		std::map<std::string, Element>::const_iterator i = elements.find(symbol);
		return i != elements.end() ? i->second : Element();
	}

}

