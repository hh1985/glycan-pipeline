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

#include<GAGPL/FRAGMENTATION/Fragmentation.h>

namespace gag
{
	// Record cleavage type and fragmentation position.
	// 1. Each time when the function is called, the validity of the cleavage should be checked.
	// 2. The composition should be updated immediately.
	void Fragment::setFragmentation(const std::string& type, const FragmentPosition& site)
	{
		CleavageCollection::iterator iter = cleavage_sites.find(type);
		if(iter != cleavage_sites.end()) // The type of cleavage has been created.
		{
			// TBD:) Validity check.
			(*iter).second.push_back(site);
		} else {
			// TBD:) Validity check.
			std::vector<FragmentPosition> site_vec;
			site_vec.push_back(site);
			cleavage_sites.insert(std::make_pair(type, site_vec));
		}
		// TBD:) Update the composition of the fragment.
	}


	// Here the Unit is a monosaccharide unit. If other types of Units are attached to it, they will
	// be treated separately.
	// No H loss during Xring cleavage.
	void Fragment::generateXType()
	{
		std::map<std::string, std::vector<FragmentPosition> >::iterator iter1 = cleavage_sites.find("X");
		if(iter1 != cleavage_sites.end()) {
			// Traversing each fragment position and produce the corresponding type of fragment.
			for(std::vector<FragmentPosition>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
			{
				//Branch& bc = glyco_seq.getBranchByID(iter2->branch_id);

				//// Remove the related internal sites.
				//// The default composition of the fragment is compo, which comes from the composition of the glycosequence.
				//Composition& cp1 = bc.getUnitByID(iter2->mono_id).getSubCompositionByRingID(iter2->xring_first, iter2->xring_second);
				//compo.deduct(cp1);

				//if(iter2->mono_id != 0) {
				//	// Remove related monosaccharide residues.
				//	compo.deduct(bc.getSubComposition(0,iter2->mono_id-1));
				//	// Remove descendant tree compositions.
				//	compo.deduct(glyco_seq.getSubTreeComposition(iter2->branch_id));
				//} else {
				//	// Selectively clean descendant branches. xring_first+1 to xring_end.
				//	for(size_t r = iter2->xring_first+1; r <= iter2->xring_second; r++)
				//	{
				//		// If the ring ID has branches attached to it.
				//		std::pair<BranchMap::right_const_iterator, BranchMap::right_const_iterator> right_range = glyco_seq.getBranchLinks().right.equal_range	(iter2->branch_id);
				//		for(BranchMap::right_const_iterator right_iter = right_range.first; right_iter != right_range.second; right_iter++)
				//		{
				//			if(right_iter->second > iter2->xring_first && right_iter->second <= iter2->xring_second)
				//			{
				//				Composition& compo_tmp = glyco_seq.getTreeComposition(right_iter->second);
				//				compo.deduct(compo_tmp);
				//			}
				//		}
				//	}
				//}
				Composition compo_tmp = this->getFragmentComposition(*iter2);
				compo.deduct(compo_tmp);

			}
		}
	}

	void Fragment::generateYType()
	{
		std::map<std::string, std::vector<FragmentPosition> >::iterator iter1 = cleavage_sites.find("Y");
		if(iter1 != cleavage_sites.end()) {
			for(std::vector<FragmentPosition>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
			{
				Composition compo_tmp = this->getFragmentComposition(*iter2);
				compo.deduct(compo_tmp);
				// Calculate the "O" acquirement.
				compo.add("O");

			}
		}

	}
		
	void Fragment::generateZType()
	{
		std::map<std::string, std::vector<FragmentPosition> >::iterator iter1 = cleavage_sites.find("Z");
		if(iter1 != cleavage_sites.end()) {
			for(std::vector<FragmentPosition>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
			{
				Composition compo_tmp = this->getFragmentComposition(*iter2);
				compo.deduct(compo_tmp);
			}
		}
	}
		
	void Fragment::generateAType()
	{
		// In fact, there will only be one A type cleavage.
		std::map<std::string, std::vector<FragmentPosition> >::iterator iter1 = cleavage_sites.find("A");
		if(iter1->second.size() > 1)
			throw std::runtime_error("There should be only one A type cleavage.");
		// In this case, we only need to consider the sum up of the composition.
		
		if(iter1 != cleavage_sites.end()) {
			for(std::vector<FragmentPosition>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
			{
				compo = this->getFragmentComposition(*iter2);
		//		Branch& bc = glyco_seq.getBranchByID(iter2->branch_id);

		//		// Remove the related internal sites.
		//		Monosaccharide& mono = bc.getUnitByID(iter2->mono_id);
		//		compo.add(mono.getSubCompositionByRingID(iter2->xring_first+1, iter2->xring_second));

		//		size_t last_id  = mono.getRingPosByRingID().second;
		//		if(iter2->xring_second < last_id) {
		//			compo.deduct(mono.getSubCompositionByRingID(iter2->xring_second+1, last_id));
		//		}

		//		if(iter2->mono_id != 0) {
		//			// Remove related monosaccharide residues.
		//			compo.add(bc.getSubComposition(0,iter2->mono_id-1));
		//			// Remove descendant tree compositions.
		//			compo.add(glyco_seq.getSubTreeComposition(iter2->branch_id));
		//		} else {
		//			// Selectively clean descendant branches. xring_first+1 to xring_end.
		//			for(size_t r = iter2->xring_first+1; r <= iter2->xring_second; r++)
		//			{
		//				// If the ring ID has branches attached to it.
		//				std::pair<BranchMap::right_const_iterator, BranchMap::right_const_iterator> right_range = glyco_seq.getBranchLinks().right.equal_range(iter2->branch_id);
		//				for(BranchMap::right_const_iterator right_iter = right_range.first; right_iter != right_range.second; right_iter++)
		//				{
		//					if(right_iter->second > iter2->xring_first && right_iter->second <= iter2->xring_second)
		//					{
		//						Composition& compo_tmp = glyco_seq.getTreeComposition(right_iter->second);
		//						compo.add(compo_tmp);
		//					}
		//				}
		//			}
		//		}

			}
		}
	}

	void Fragment::generateBType()
	{
		std::map<std::string, std::vector<FragmentPosition> >::iterator iter1 = cleavage_sites.find("Y");
		if(iter1 != cleavage_sites.end()) {
			for(std::vector<FragmentPosition>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
			{
				compo = this->getFragmentComposition(*iter2);
				// Calculate the "O" acquirement.
				compo.deduct("O");
			}
		}

	}
		
	void Fragment::generateCType()
	{
		std::map<std::string, std::vector<FragmentPosition> >::iterator iter1 = cleavage_sites.find("Z");
		if(iter1 != cleavage_sites.end()) {
			for(std::vector<FragmentPosition>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
			{
				compo = this->getFragmentComposition(*iter2);
			}
		}
	}

	// Generate generic segment from NRE.
	// The start id can be and endid can be the 
	Composition Fragment::getFragmentComposition(FragmentPosition& fp)
	{
		Branch& bc = glyco_seq.getBranchByID(fp.branch_id);

		Composition tmp_compo;
		Composition& cp1 = bc.getUnitByID(fp.mono_id).getSubCompositionByRingID(fp.xring_first, fp.xring_second);
		tmp_compo.add(cp1);

		if(fp.mono_id != 0) {
			// Remove related monosaccharide residues.
			tmp_compo.add(bc.getSubComposition(0,fp.mono_id-1));
			// Remove descendant tree compositions.
			tmp_compo.add(glyco_seq.getSubTreeComposition(fp.branch_id));
		} else {
			// Selectively clean descendant branches. xring_first+1 to xring_end.
			for(size_t r = fp.xring_first+1; r <= fp.xring_second; r++)
			{
				// If the ring ID has branches attached to it.
				std::pair<BranchMap::right_const_iterator, BranchMap::right_const_iterator> right_range = glyco_seq.getBranchLinks().right.equal_range(fp.branch_id);
				for(BranchMap::right_const_iterator right_iter = right_range.first; right_iter != right_range.second; right_iter++)
				{
					if(right_iter->second > fp.xring_first && right_iter->second <= fp.xring_second)
					{
						//Composition& compo_tmp = glyco_seq.getTreeComposition(right_iter->second);
						tmp_compo.add(glyco_seq.getTreeComposition(right_iter->second));
					}
				}
			}
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
		if(iter1->second.size() > 1 && (type == "A" || type == "B" || type == "C")) {
			throw std::runtime_error("There should be only one cleavage on the reducing end!");
		}
		
		if(iter1 != cleavage_sites.end()) {
			Composition tmp_compo;
			for(std::vector<FragmentPosition>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++)
			{
				tmp_compo.add(this->getFragmentComposition(*iter2));
			}
			
			if(type == "X" || type == "Y" || type == "Z") { // The RE part will be kept.
				compo.deduct(tmp_compo);
			} else {				// The NRE part will be kept.
				compo = tmp_compo;
			}

			// Add composition shift. 
			// Only cleavage type is considered. The dissociation type will be dealed when mathching against library.
			FragmentationParams& fp = ft.getFragmentationParams(type);
			CompositionShift& cs = ft.getCleavageShift(type); 
			
			if(!cs.second.empty()) {
				if(cs.first == 1)
					compo.add(cs.second);
				else
					compo.deduct(cs.second);
			}

		} else {
			// Do nothing.
		}

	}

	// There will be a H loss for B,Y,C,Z type.
	void correctMassForCAD(const std::string& dis)
	{
		
	}

}
