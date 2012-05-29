#ifndef  GAG_BRANCH_INC
#define  GAG_BRANCH_INC

#include <algorithm>
#include <GAGPL/GLYCAN/Linkage.h>
#include <GAGPL/GLYCAN/Monosaccharide.h>

namespace gag
{
	class Branch : public Unit
	{
		// using boost::any_cast;
		// The connecting branch ids and the linkage.
		// typedef std::map<size_t, Linkage> LinkEnv; 

		private:
			size_t branch_id;
			
			std::vector<Monosaccharide> mono_chain;

			// For glycan, there will be only one reducing end.
			//LinkEnv nre_branch_ids;
			//LinkEnv re_branch_ids; 
			
			// The mono_id and the linkage.
			std::vector<Linkage> links;

		public:
			Branch(const size_t id)
				: branch_id(id), mono_chain(), links()
			{}
			Branch(){}

			void addUnit(Monosaccharide& mono_unit);
			
			void addLinkage(const Linkage& link);

			std::vector<Linkage> getNeighborLinks(const size_t mono_id);

			inline size_t getBranchID() const
			{
				return branch_id;
			}

			inline std::vector<Monosaccharide>& getGlycanChainUnits()
			{
				return mono_chain;
			}

			inline std::vector<Linkage>& getLinkages()
			{
				return links;
			}

			inline Monosaccharide& getUnitByID(const size_t id)
			{
				return mono_chain.at(id);
			}
			
			// Calculate the mass and composition. The correction of the terminal mass will be considered.
			void update();

			// Useful for calculating fragment mass.
			Composition getSubComposition(const size_t start, const size_t end);
			
			void addModification(const size_t mono_id, const size_t site_id, FunctionalGroup& ori, const Composition& plus, const Composition& minus);
			void addModification(const size_t mono_id, const size_t site_id, const std::string& ori, const std::string& plus, const std::string& minus);

			void removeModification(const size_t mono_id, const size_t site_id, FunctionalGroup& ori, const Composition& minus);
			void removeModification(const size_t mono_id, const size_t site_id, const std::string& ori, const std::string& minus);
	};
}

#endif   /* ----- #ifndef GAG_BRANCH_INC  ----- */
