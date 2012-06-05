/*
 * =====================================================================================
 *
 *       Filename:  Fragmentation.cpp
 *
 *    Description:  The implement file for class Fragment.
 *
 *        Version:  1.0
 *        Created:  05/ 8/2012  2:24:18 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <set>

namespace gag
{
 // // if true, simply remove the corresponding part instead of recalculation.
	// Only works for A type.
	bool Fragment::isInternalCleavage(const FragmentPosition& fp)
	{
		Branch& bc = glyco_seq.getBranchByID(fp.branch_id);
		if(fp.mono_id != 0) { // Within branch.
			if(bc.getLinkages().size() < fp.mono_id+1) // Reducing end;
				return true;
			if((bc.getLinkages().at(fp.mono_id-1).end > bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_second)&& bc.getLinkages().at(fp.mono_id).start <= bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_first)) ||
				(bc.getLinkages().at(fp.mono_id-1).end <= bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_second)&& bc.getLinkages().at(fp.mono_id).start > bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_first)))
				return true;
			else 
				return false;
		} else {
			BranchMap::right_const_iterator right_iter = glyco_seq.getBranchLinks().right.find(fp.branch_id);
			if(right_iter != glyco_seq.getBranchLinks().right.end()){ // No leaf
				if((glyco_seq.getBranchByID(fp.branch_id-1).getLinkages().back().end > bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_second) && bc.getLinkages().at(fp.mono_id).start <= bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_first)) ||
					(glyco_seq.getBranchByID(fp.branch_id-1).getLinkages().back().end <= bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_second) && bc.getLinkages().at(fp.mono_id).start > bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_first)))
					return true;
				else 
					return false;
			} else { // Leaf.
				if(bc.getLinkages().at(fp.mono_id).start > bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_first) && 
					bc.getLinkages().at(fp.mono_id).start < bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_second))
					return true;
				else 
					return false;
			}
		}

	}
	// Record cleavage type and fragmentation position.
	// 1. Each time when the function is called, the validity of the cleavage should be checked.
	// 2. The composition should be updated immediately.
	void Fragment::setFragmentation(const std::string& type, const FragmentPosition& site)
	{
		CleavageCollection::iterator iter1 = cleavage_sites.find(type);
		if(iter1 != cleavage_sites.end()) // The type of cleavage has been created.
		{
			(*iter1).second.push_back(site);
		} else {
			std::vector<FragmentPosition> site_vec;
			site_vec.push_back(site);
			cleavage_sites.insert(std::make_pair(type, site_vec));
		}
		/* Update the composition. */
		
		//Composition tmp_compo = (*this).updateFragmentComposition(site);



		Branch& bc = glyco_seq.getBranchByID(site.branch_id);
		Composition ori_compo = glyco_seq.getComposition();
		Composition tmp_compo1 = this->updateFragmentComposition(site);
		Composition tmp_compo2 = this->updateFragmentComposition(site, false);
		//size_t pos1 = bc.getUnitByID(site.mono_id).getCarbonID(site.xring_first);
		//size_t pos2 = bc.getLinkages().at(site.mono_id).start;
		if(type == "Y" || type == "Z") { // The RE part will be kept.
			compo.deduct(tmp_compo1);
		} else if(type == "B" || type == "C") {				// The NRE part will be kept.
			ori_compo.deduct(tmp_compo1);
			compo.deduct(ori_compo);
		} else if(type == "X") {
			if(bc.getUnitByID(site.mono_id).getCarbonID(site.xring_first) < bc.getUnitByID(site.mono_id).getRingStart()){
				compo.deduct(tmp_compo2);
			} else {
				compo.deduct(tmp_compo1);
			}
		} else if(type == "A") {
			if(bc.getUnitByID(site.mono_id).getCarbonID(site.xring_first) < bc.getUnitByID(site.mono_id).getRingStart()){
				ori_compo.deduct(tmp_compo2);
				compo.deduct(ori_compo);
			} else {
				//compo.add(this->updateFragmentComposition(site, false));
				ori_compo.deduct(tmp_compo1);
				compo.deduct(ori_compo);
			}
		} 
		
		FragmentationParams& fp = ft.getFragmentationParams(type);
		CompositionShift& cs = ft.getCleavageShift(type); 
		for(CompositionShift::iterator iter = cs.begin(); iter != cs.end(); iter++) {
			if(!iter->second.empty()) {
				if(iter->first == 1)
					compo.add(iter->second);
				else if(iter->first == -1)
					compo.deduct(iter->second);
			}
		}
	}

	// Generate generic segment from NRE for each fragmentation site.
	// The start id can be and endid can be the 
	Composition Fragment::updateFragmentComposition(const FragmentPosition& fp, bool cw)
	{
		Branch& bc = glyco_seq.getBranchByID(fp.branch_id);
		Monosaccharide& ms = glyco_seq.getBranchByID(fp.branch_id).getUnitByID(fp.mono_id);

		Composition tmp_compo;
		Composition& cp1 = bc.getUnitByID(fp.mono_id).getSubCompositionByRingID(fp.xring_first, fp.xring_second);
		if(cw == true){
			tmp_compo.add(cp1);
			if(fp.mono_id != 0) {
				if(bc.getLinkages().at(fp.mono_id-1).end <= bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_second) || fp.xring_second == 0) {
					tmp_compo.add(glyco_seq.getSubTreeComposition(fp.branch_id));
					tmp_compo.add(bc.getSubComposition(0,fp.mono_id-1));
				}
				
			} else {
				BranchMap::right_const_iterator right_iter = glyco_seq.getBranchLinks().right.find(fp.branch_id);
				BOOST_FOREACH(BranchMap::right_reference right_ref, glyco_seq.getBranchLinks().right.equal_range(fp.branch_id))
				{
					size_t hinge_id = glyco_seq.getBranchByID(right_ref.second).getLinkages().back().end;

					if(hinge_id > ms.getCarbonID(fp.xring_first) && hinge_id <= ms.getCarbonID(fp.xring_second))
						tmp_compo.add(glyco_seq.getTreeComposition(right_ref.second));
				}
			}
		} else { // Anti-clockwise
			Composition mono_compo = bc.getUnitByID(fp.mono_id).getComposition();
			tmp_compo.add(mono_compo);
			tmp_compo.deduct(cp1);
			if(fp.mono_id != 0){
				if(bc.getLinkages().at(fp.mono_id-1).end > bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_second)) {
					tmp_compo.add(glyco_seq.getSubTreeComposition(fp.branch_id));
					tmp_compo.add(bc.getSubComposition(0,fp.mono_id-1));
				}
			} else {
				BranchMap::right_const_iterator right_iter = glyco_seq.getBranchLinks().right.find(fp.branch_id);
				BOOST_FOREACH(BranchMap::right_reference right_ref, glyco_seq.getBranchLinks().right.equal_range(fp.branch_id))
				{
					size_t hinge_id = glyco_seq.getBranchByID(right_ref.second).getLinkages().back().end;

					if(hinge_id > ms.getCarbonID(fp.xring_second))
						tmp_compo.add(glyco_seq.getTreeComposition(right_ref.second));
				}
			}
		}

		return tmp_compo;
	}


	
	Fragment& Fragment::operator=(const Fragment& rhs)
	{
		if(this != &rhs)
		{
			cleavage_sites = rhs.cleavage_sites;
		}
		return *this;
	}
}
