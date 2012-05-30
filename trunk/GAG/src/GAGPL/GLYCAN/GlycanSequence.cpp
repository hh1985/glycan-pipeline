/*
 * =====================================================================================
 *
 *       Filename:  GlycanSequence.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/ 5/2012 10:15:09 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#include <GAGPL/GLYCAN/GlycanSequence.h>

namespace gag
{
	void GlycanSequence::addBranch(Branch& bc)
	{
		glycan_sequence.push_back(bc);
		compo.add(bc.getComposition());
	}
	void GlycanSequence::addBranchLink(const size_t re_id, const size_t nre_id)
	{
		branch_links.insert(BranchMap::value_type(re_id, nre_id));
		// Update composition.
		Branch& bc = (*this).getBranchByID(nre_id);
		Linkage& lk = (*this).getBranchByID(re_id).getLinkages().back();
		bc.removeModification(0,lk.end,"OH","OH");
		compo.deduct("OH");
	}

	void GlycanSequence::updateChildrenIDs(const size_t branch_id)
	{
		BranchMap::right_const_iterator right_lower_iter = branch_links.right.lower_bound(branch_id);

		size_t min = branch_id;
		BranchDescendants::const_iterator d_iter = children.find(right_lower_iter->second);
		
		if(d_iter == children.end())
			min = right_lower_iter->second;
		else
			min = d_iter->second;

		children.insert(std::make_pair(branch_id, min));
	}

	void GlycanSequence::update()
	{
		compo.clear();
		
		for(std::vector<Branch>::iterator iter = glycan_sequence.begin(); 
				iter!=glycan_sequence.end(); iter++)
		{
			compo.add(iter->getComposition());
		}
	}

	void GlycanSequence::addModification(const size_t bc_id, const size_t mono_id, const size_t site_id, FunctionalGroup& fg, const Composition& modi, const Composition& loss)
	{
		Branch& bc = glycan_sequence.at(bc_id);
		bc.addModification(mono_id, site_id, fg, modi, loss);
		compo.add(modi);
		compo.deduct(loss);
	}
	void GlycanSequence::addModification(const size_t bc_id, const size_t mono_id, const size_t site_id, const std::string& ori, const std::string& plus, const std::string& loss)
	{
		Branch& bc = glycan_sequence.at(bc_id);
		bc.addModification(mono_id, site_id, ori, plus, loss);
		compo.add(plus);
		compo.deduct(loss);
	}

	std::pair<size_t, size_t> GlycanSequence::getDescendantBranchIDs(const size_t& branch_id)
	{	
		BranchDescendants::const_iterator d_iter = children.find(branch_id);
		if(d_iter != children.end())
			return std::pair<size_t, size_t>(d_iter->second, branch_id-1);
		else // Leaf branch.
			return std::pair<size_t, size_t>(branch_id, branch_id-1);
	}

	Composition GlycanSequence::getSubComposition(const size_t& id1, const size_t& id2)
	{
		Composition sub_compo;
		for(size_t i = id1; i < id2+1; i++)
		{
			sub_compo.add((*this).getBranchByID(i).getComposition());
		}
		return sub_compo;
	}

	Composition GlycanSequence::getSubTreeComposition(const size_t& branch_id)
	{
		std::pair<size_t, size_t>& range = (*this).getDescendantBranchIDs(branch_id);

			return (*this).getSubComposition(range.first, range.second);
	}
	
	Composition GlycanSequence::getTreeComposition(const size_t& branch_id)
	{
		std::pair<size_t, size_t>& range = (*this).getDescendantBranchIDs(branch_id);

			return (*this).getSubComposition(range.first, range.second+1);
	}
}
