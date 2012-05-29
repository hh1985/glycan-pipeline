/*
 * =====================================================================================
 *
 *       Filename:  Branch.cpp
 *
 *    Description:  Implementation file of class Branch.
 *
 *        Version:  1.0
 *        Created:  05/ 4/2012  9:07:08 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/GLYCAN/Branch.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>

namespace gag
{

	void Branch::addUnit(Monosaccharide& mono_unit)
	{	
		
		// Remove OH at link.end.
		if(mono_chain.size() != 0)
		{
			Linkage& lk = links.back();
			//Composition& cp = is.getFunctionalGroups().at(0).getComposition();
			FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
			//:) TBD The functional group should be specified by name.
			FunctionalGroup& fg = fgt.getFunctionalGroupBySymbol("OH");
			mono_unit.remove(lk.end, fg);
		}
		mono_chain.push_back(mono_unit);
		compo.add(mono_unit.getComposition());

	}

	void Branch::addLinkage(const Linkage& link)
	{
		// Modify the unit.
		//std::vector<Monosaccaride>::iterator iter = mono_chain.back();
		//Monosaccharide& mono = mono_chain.back();
		Monosaccharide& mono = mono_chain.at(link.nre_id);
		
		//InternalSite& is = mono.getInternalSites().at(link.start); 
		FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
		//:) TBD The functional group should be specified by name.
		FunctionalGroup& fg = fgt.getFunctionalGroupBySymbol("H");
		mono.remove(link.start, fg);
		// Retrieve the composition on specified position.
		//Composition& cp = is.getFunctionalGroups().at(0).site_gps.at(1);
		//FunctionalGroup fg("OH");
		//is.remove(fg, "H");

		links.push_back(link);
		compo.deduct(fg.getComposition());
	}

	std::vector<Linkage> Branch::getNeighborLinks(const size_t mono_id)
	{
		// The structure is guaranteed to be linear.
		std::vector<Linkage> lk;
		if(links.size() == 1){
		}
		else if(mono_id == links.size()-1) {
			lk.push_back(links.at(mono_id-1));
		} else if(mono_id == 0) {
			lk.push_back(links.at(mono_id));
		} else {
			lk.push_back(links.at(mono_id - 1));
			lk.push_back(links.at(mono_id));
		}
		return lk;
	}

	Composition Branch::getSubComposition(const size_t start, const size_t end)
	{
		Composition sub_compo;

		for(std::vector<Monosaccharide>::iterator iter = mono_chain.begin() + start; 
				iter != mono_chain.begin()+end+1; iter++)
		{
			sub_compo.add(iter->getComposition());
		}

		return sub_compo;
	}
	
	// Recalculate the overall mass.
	void Branch::update()
	{
		compo.clear();

		for(std::vector<Monosaccharide>::iterator iter = mono_chain.begin(); 
				iter!=mono_chain.end(); iter++)
		{
			compo.add(iter->getComposition());
		}

	}

	void Branch::addModification(const size_t mono_id, const size_t site_id, FunctionalGroup& ori, const Composition& plus, const Composition& minus)
	{
		Monosaccharide& mono = mono_chain.at(mono_id);
		mono.add(site_id, ori, plus, minus);
		compo.add(plus);
		compo.deduct(minus);
	}

	void Branch::addModification(const size_t mono_id, const size_t site_id, const std::string& ori, const std::string& plus, const std::string& minus)
	{
		Monosaccharide& mono = mono_chain.at(mono_id);
		mono.add(site_id, ori, plus, minus);
		compo.add(plus);
		compo.deduct(minus);
	}

	void Branch::removeModification(const size_t mono_id, const size_t site_id, FunctionalGroup& ori, const Composition& minus)
	{
		Monosaccharide& mono = mono_chain.at(mono_id);
		mono.remove(site_id, ori, minus);
		compo.deduct(minus);
	}
	void Branch::removeModification(const size_t mono_id, const size_t site_id, const std::string& ori, const std::string& minus)
	{
		Monosaccharide& mono = mono_chain.at(mono_id);
		mono.remove(site_id, ori, minus);
		compo.deduct(minus);
	}
}





