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
	// Record cleavage type and fragmentation position.
	// 1. Each time when the function is called, the validity of the cleavage should be checked.
	// 2. The composition should be updated immediately.
	void Fragment::setFragmentation(const std::string& type, const FragmentPosition& site)
	{
		CleavageCollection::iterator iter1 = cleavage_sites.find(type);
		if(iter1 != cleavage_sites.end()) // The type of cleavage has been created.
		{
			//// Validity check.			
			//if(type == 'A' || type == 'B' || type == 'C') {
			//	// 1. For A/B/C cleavage, there should be only in the whole set.
			//	if(cleavage_sites.count('A') > 0 || cleavage_sites.count('B') > 0 || cleavage_sites.count('C') > 0)
			//		throw std::runtime_error("The number of A/B/C cleavage should be at most 1.");

			//} else if(type == 'X' || type == 'Y' || type == 'Z'){
			//	// 2. For X/Y/Z cleavage, every branch can has at most one
			//	for(char temp_type = 'X'; temp_type <= 'Z'; temp_type++)
			//	{
			//		CleavageCollection::iterator type_iter = cleavage_sites.find(temp_type);
			//		if(type_iter != cleavage_sites.end()) {
			//			// For each FragmentationPosition, check if there is any overlapping between branch_ids.
			//			for(std::vector<FragmentPosition>::iterator pos_iter = type_iter->second.begin(); pos_iter != type_iter->second.end(); pos_iter++)
			//			{
			//				if(site.branch_id == pos_iter->branch_id)
			//					throw std::runtime_error("No X/Y/Z cleavages allowed to occur on the same branch simultaneously.");
			//			}
			//		}
			//	}
			//	// 3. X/Y/Z can only appear in the downstream of the NRE cleavage.
			//	if(cleavage_sites.count('A') > 0 || cleavage_sites.count('B') > 0 || cleavage_sites.count('C') > 0)
			//	
			//}
			(*iter1).second.push_back(site);
		} else {
			std::vector<FragmentPosition> site_vec;
			site_vec.push_back(site);
			cleavage_sites.insert(std::make_pair(type, site_vec));
		}
		/* Update the composition. */
		Composition ori_compo = glyco_seq.getComposition();
		Composition tmp_compo = (*this).updateFragmentComposition(site);

		FragmentationParams& fp = ft.getFragmentationParams(type);
		CompositionShift& cs = ft.getCleavageShift(type); 
		for(CompositionShift::iterator iter = cs.begin(); iter != cs.end(); iter++) {
			if(!iter->second.empty()) {
				if(iter->first == 1)
					tmp_compo.add(iter->second);
				else if(iter->first == -1)
					tmp_compo.deduct(iter->second);
			}
		}

		if(type == "X" || type == "Y" || type == "Z") { // The RE part will be kept.
			compo.deduct(tmp_compo);
		} else {				// The NRE part will be kept.
			ori_compo.deduct(tmp_compo);
			compo.deduct(ori_compo);
		}
	}

	// Generate generic segment from NRE for each fragmentation site.
	// The start id can be and endid can be the 
	Composition Fragment::updateFragmentComposition(const FragmentPosition& fp)
	{
		Branch& bc = glyco_seq.getBranchByID(fp.branch_id);

		Composition tmp_compo;
		Composition& cp1 = bc.getUnitByID(fp.mono_id).getSubCompositionByRingID(fp.xring_first, fp.xring_second);
		tmp_compo.add(cp1);

		if(fp.mono_id != 0) { 
			// It is necessory to judge if there is any internal cleavage here.
			if(bc.getLinkages().at(fp.mono_id-1).end > bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_second)
				&& bc.getLinkages().at(fp.mono_id).start <= bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_first))
			{
				// Do nothing here.
			} else {
				// Remove related monosaccharide residues.
				tmp_compo.add(bc.getSubComposition(0,fp.mono_id-1));
				// Remove descendant tree compositions.
				tmp_compo.add(glyco_seq.getSubTreeComposition(fp.branch_id));
			}
		} else {
			// Selectively clean descendant branches. xring_first+1 to xring_end.
			//for(size_t r = fp.xring_first+1; r <= fp.xring_second; r++)
			//{
				// If the ring ID has branches attached to it.
				//std::pair<BranchMap::right_const_iterator, BranchMap::right_const_iterator> right_range = glyco_seq.getBranchLinks().right.equal_range(fp.branch_id);
				//for(BranchMap::right_const_iterator right_iter = right_range.first; right_iter != right_range.second; right_iter++)
			// Consider about the possibility of internal cleavage.
			BranchMap::right_const_iterator right_iter = glyco_seq.getBranchLinks().right.find(fp.branch_id);
			// Not the leaf node.
			if(right_iter != glyco_seq.getBranchLinks().right.end())
			{
				if(glyco_seq.getBranchByID(fp.branch_id-1).getLinkages().back().end > bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_second)
					&& bc.getLinkages().at(fp.mono_id).start <= bc.getUnitByID(fp.mono_id).getCarbonID(fp.xring_first))
				{
					// Do nothing here
				} else {
					BOOST_FOREACH(BranchMap::right_reference right_ref, glyco_seq.getBranchLinks().right.equal_range(fp.branch_id))
					{
						size_t hinge_id = glyco_seq.getBranchByID(right_ref.second).getLinkages().back().end;
						Monosaccharide& ms = glyco_seq.getBranchByID(fp.branch_id).getUnitByID(fp.mono_id);
						if(ms.getRingID(hinge_id) > fp.xring_first && ms.getRingID(hinge_id) <= fp.xring_second)
						{
							tmp_compo.add(glyco_seq.getTreeComposition(right_ref.second));
						}
					}				
				}
			}
				
			//}
		}

		return tmp_compo;
	}
	
	// 1. Base on the cleavage type, decide which part to keep and which part to discard.
	// 2. Base on the cleavage type, append the composition shift.
	void Fragment::updateFragmentByType(const std::string& type)
	{
		std::map<std::string, std::vector<FragmentPosition> >::iterator iter1 = cleavage_sites.find(type);
		// For A/B/C cleavage, there should be only one cleavage.
		// For X/Y/Z cleavage, there might be more than one cleavage.
		/*if(iter1->second.size() > 1 && (type == "A" || type == "B" || type == "C")) {
			throw std::runtime_error("There should be only one cleavage on the reducing end!");
		}*/
		
		if(iter1 != cleavage_sites.end()) {
			Composition ori_compo = glyco_seq.getComposition();
			Composition tmp_compo;
			// Add composition shift. 
			// Only cleavage type is considered. The dissociation type will be dealed when mathching against library.
			FragmentationParams& fp = ft.getFragmentationParams(type);
			CompositionShift& cs = ft.getCleavageShift(type); 

			for(std::vector<FragmentPosition>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
			{
				tmp_compo.add(this->updateFragmentComposition(*iter2));
				for(CompositionShift::iterator iter = cs.begin(); iter != cs.end(); iter++) {
					if(!iter->second.empty()) {
						if(iter->first == 1)
							tmp_compo.add(iter->second);
						else if(iter->first == -1)
							tmp_compo.deduct(iter->second);
					}
				}
			}
			
			if(type == "X" || type == "Y" || type == "Z") { // The RE part will be kept.
				compo.deduct(tmp_compo);
			} else {				// The NRE part will be kept.
				ori_compo.deduct(tmp_compo);
				compo.deduct(ori_compo);
			}

		} else {
			// Do nothing.
		}

	}

	//void Fragment::update()
	//{
	//	///* Validity check. */
	//	//// 1. At most one NRE cleavage.
	//	//size_t re_count = 0;
	//	//char flag_type;
	//	//for(char temp_type = 'A'; temp_type <= 'C'; temp_type++)
	//	//{
	//	//	CleavageCollection::iterator iter = cleavage_sites.find(temp_type);
	//	//	if(iter != cleavage_sites.end()) {
	//	//		re_count += iter->second.size();
	//	//		flag_type = temp_type;
	//	//	}
	//	//}
	//	//if(re_count > 1)
	//	//	throw std::runtime_error("The number of A/B/C cleavage should be at most 1.");

	//	//// 2. At most one RE cleavage on each branch.
	//	//std::set<size_t> ids;
	//	//// The downid should include the id where NRE happens.
	//	//std::set<size_t> down_ids;
	//	//for(char temp_type = 'X'; temp_type <= 'Z'; temp_type++)
	//	//{
	//	//	CleavageCollection::iterator iter = cleavage_sites.find(temp_type);
	//	//	if(iter != cleavage_sites.end()) {
	//	//		std::vector<FragmentPosition>::iterator pos_iter = iter->second.begin();
	//	//		for(; pos_iter != iter->second.end(); pos_iter++)
	//	//		{
	//	//			std::set<size_t>::iterator set_iter = ids.find(pos_iter->branch_id);
	//	//			if(set_iter != ids.end())
	//	//				throw std::runtime_error("No X/Y/Z cleavages allowed to occur on the same branch simultaneously.");
	//	//			else if(re_count == 1 && !down_ids.empty()) {
	//	//				std::set<size_t>::iterator down_iter = down_ids.find(pos_iter->branch_id);
	//	//				if(down_iter == down_ids.end())
	//	//					throw std::runtime_error("X/Y/Z cleavage has to be on the subtree of A/B/C cleavage");
	//	//				// If the NRE cleavage has the same branch id as the RE cleavage.
	//	//				if(pos_iter->branch_id == )
	//	//			}
	//	//			ids.insert(pos_iter->branch_id);
	//	//
	//	//		}
	//	//	}
	//	//}
	//	//for(char temp_type = 'A'; temp_type < 'C'; temp_type++)
	//	//{
	//	//	CleavageCollection::iterator iter = cleavage_sites.find(temp_type);
	//	//	if(iter != cleavage_sites.end()) {
	//	//		updateFragmentByType(temp_type);
	//	//		// Assume the validity check has been passed.
	//	//		break;
	//	//	}
	//	//}
	//	CleavageCollection::iterator iter = cleavage_sites.begin();

	//	for(; iter != cleavage_sites.end(); iter++)
	//	{
	//		(*this).updateFragmentByType(iter->first);
	//	}
	//}

}
