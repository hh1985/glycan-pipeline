/*
 * =====================================================================================
 *
 *       Filename:  Fragmentation.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  4/23/2012 4:39:53 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef  GAG_FRAGMENTATION_H_INC
#define  GAG_FRAGMENTATION_H_INC

#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/FRAGMENTATION/FragmentationTable.h>
//#include <boost/tuple/tuple.hpp>

namespace gag
{
  //typedef std::vector<std::vector<Branch> > Fragments;
  
	struct FragmentPosition
	{
		// The branch ID should follow the leaf-to-root style.
		size_t branch_id;
		size_t mono_id;
		// Ring ID. Please be careful of the conversion from Carbon ID to
		// Ring ID.
		size_t xring_first;
		size_t xring_second;
	};


	//using boost::tuples;
	// Cleavage type, fragmentation position.
	typedef std::map<std::string, std::vector<FragmentPosition> > CleavageCollection;

	class Fragment: public Unit
  {
    private:
			// The sequence should not be modified.
			// All calculation is based on position.
			GlycanSequence& glyco_seq;
			FragmentationTable& ft;

			//std::vector<FragmentPosition> cleavage_sites;
			// Cleavage type and site.
			//std::map<std::string, std::vector<FragmentPosition> > cleavage_sites;
			CleavageCollection cleavage_sites;
			//CleavageShift cleavage_correction;
		
		protected:
			// A generic function which can be called by generate[X]Type() ([X] = X,Y,Z,A,B,C here)
			// using information from cleavage_sites.
			Composition updateFragmentComposition(const FragmentPosition& fp);
		
		public:
			// Constructor.
			Fragment(GlycanSequence& seq)
				: glyco_seq(seq), ft(FragmentationTable::Instance()), cleavage_sites(), Unit(seq.getComposition())
			{
			}

			// Check if the cleavage is a type of internal cleavage.
			bool InternalCleavage(const FragmentPosition& site);

			// Repeatedly store all the fragmentation information.
			void setFragmentation(const std::string& type, const FragmentPosition& site);

			// A generic function for fragmentation.
			void updateFragmentByType(const std::string& type);
	
			// The update process can be divided into several steps:
			// 1. Locate the cleavage position using information from Fragmentation.
			// 2. Update the mass and composition of that Branch.
			// 3. Clear the composition and mass of all related branch.
			// 4. Recalculate the overall composition and mass.
			// void update();


  };


}





#endif   /* ----- #ifndef GAG_FRAGMENTATION_H_INC  ----- */
