/*
 * =====================================================================================
 *
 *       Filename:  InternalSite.h
 *
 *    Description:  The internal site inside a monosaccharide unit.
 *
 *        Version:  1.0
 *        Created:  04/30/2012  3:50:55 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#ifndef  GAG_INTERNALSITE_H
#define  GAG_INTERNALSITE_H

#include <GAGPL/CHEMISTRY/FunctionalGroup.h>
#include <set>

namespace gag
{

	class InternalSite: public Unit 
	{
		private:
			// This information should come from the xml file.
			// But the CompositionGroup should be a copy instead of singleton.
			// The only reason using FunctionalGroup instead of Composition is that
			// The name can be used for verification purposes.
			std::multiset<FunctionalGroup> _groups;
		
		public:
	
		InternalSite() {}
		InternalSite(Composition& compo)
			: Unit(compo)
		{}
		InternalSite(std::string str)
			: Unit(str)
		{}

		inline void addFunctionalGroup(const FunctionalGroup& fg)
		{
			// TBD: Exception of repeating keys.
			//_groups.insert(std::make_pair(fg, num));
			_groups.insert(fg);
		}

		// Specify the original fg, the fg that will add to it, and the lost composition of all.
		// High level.
		void modify(FunctionalGroup& ori, FunctionalGroup& sub, Composition& loss);
		// Specify the original fg, and the overall composition that will add to it.
		// Low level.
		//void add(FunctionalGroup& ori, const Composition& compo);

		// The problem is: how to specify which functionalgroup it is going to remove.
		// If there is repeating fgs.
		//void remove(FunctionalGroup& ori, const Composition& compo);
		//void remove(FunctionalGroup& ori);

		inline std::multiset<FunctionalGroup>& getFunctionalGroups()
		{
			return _groups;
		}
		
	};
}

#endif   /* ----- #ifndef GAG_INTERNALSITE_H ---- */
