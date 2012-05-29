/*
 * =====================================================================================
 *
 *       Filename:  FunctionGroup.h
 *
 *    Description:  FunctionGroup
 *
 *        Version:  1.0
 *        Created:  4/23/2012 3:11:45 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Han Hu 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef  GAG_FUNCTIONALGROUP_H
#define  GAG_FUNCTIONALGROUP_H

#include <GAGPL/CHEMISTRY/Composition.h>
#include <GAGPL/CHEMISTRY/Unit.h>

namespace gag
{

	class FunctionalGroup: public Unit
	{
		// This should be ordered from NRE to RE.
		typedef std::vector<Composition> Sites;
		
		public:
			std::string symbol;
			std::string name;
			// Sites can be easily converted to cleavage units.
			Sites site_gps;
			// The functionalgroup on the substrate and composition loss.
			//std::map<FunctionalGroup, Composition> rules;
		
			// The function is useful when the functional group is treated as a combination of units for fragmentation.
			inline std::vector<Unit> convertToUnits()
			{
				std::vector<Unit> units;

				for(std::vector<Composition>::iterator iter = site_gps.begin(); iter!=site_gps.end(); iter++)
				{
					Unit ut(*iter);
					units.push_back(ut);
				}
				return units;
			}
			
			friend bool operator<(const FunctionalGroup& left, const FunctionalGroup& right)
			{
				return left.symbol < right.symbol;
			}
	};
}

#endif   /* ----- #ifndef GAG_FUNCTIONALGROUP_H_INC  ----- */
