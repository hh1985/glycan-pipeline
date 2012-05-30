/*
 * =====================================================================================
 *
 *       Filename:  GlycanSequence.h
 *
 *    Description:  The class for building glycan chain.
 *
 *        Version:  1.0
 *        Created:  05/ 3/2012  3:49:50 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#ifndef  GAG_GLYCANSEQUENCE_INC
#define  GAG_GLYCANSEQUENCE_INC

#include <algorithm>
#include <GAGPL/GLYCAN/Branch.h>
//#include <GAGPL/CHEMISTRY/FunctionalGroup.h>
#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/foreach.hpp>

namespace gag
{
	typedef boost::bimap< boost::bimaps::set_of<size_t, std::less<size_t> >, boost::bimaps::multiset_of<size_t> > BranchMap;
	typedef BranchMap::left_map::const_iterator left_const_iterator;
	typedef BranchMap::right_map::const_iterator right_const_iterator;

	typedef std::map<size_t, size_t> BranchDescendants;

	class GlycanSequence: public Unit
	{	

		private:
			std::vector<Branch> glycan_sequence;
			// The left value should be the NRE branch and the right value should be the RE branch.
			BranchMap branch_links;
			// The key is the target branch id, the value is lowest children id, the range of the children ids is [lowest, branch_id -1].
			BranchDescendants children;

		public:
			// When a new branch is added, the masses of the corresponding 
			// monosaccharide unit will be modified.
			void addBranch(Branch&);

			void addBranchLink(const size_t re_id, const size_t nre_id);

			void updateChildrenIDs(const size_t branch_id);

			inline std::vector<Branch>& getBranches()
			{
				return glycan_sequence;
			}

			inline Branch& getBranchByID(const size_t id)
			{
				return glycan_sequence.at(id);
			}

			inline BranchMap& getBranchLinks()
			{
				return branch_links;
			}

			void addModification(const size_t bc_id, const size_t mono_id, const size_t site_id, FunctionalGroup& fg, const Composition& modi, const Composition& loss);
			void addModification(const size_t bc_id, const size_t mono_id, const size_t site_id, const std::string& ori, const std::string& plus, const std::string& loss);

			void update();

			// The second value is used to evaluate if this is a leaf branch.
			std::pair<size_t, size_t> getDescendantBranchIDs(const size_t& branch_id);

			Composition getSubComposition(const size_t& id1, const size_t& id2);

			Composition getSubTreeComposition(const size_t& branch_id);
			Composition getTreeComposition(const size_t& branch_id);

			
	};
}

#endif   /* ----- #ifndef GAG_GLYCANSEQUENCE_INC  ----- */
