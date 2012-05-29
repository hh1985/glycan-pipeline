/*
 * =====================================================================================
 *
 *       Filename:  FunctionalGroupTable.h
 *
 *    Description:  Functional group table loaded from xml drive.
 *
 *        Version:  1.0
 *        Created:  04/29/2012  9:51:00 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */


#ifndef	GAG_FUNCTIONALGROUPTABLE_H
#define GAG_FUNCTIONALGROUPTABLE_H

#include <GAGPL/MISC/ConfigLoader.h>
#include <GAGPL/CHEMISTRY/FunctionalGroup.h>
#include <boost/noncopyable.hpp>

namespace gag
{
	class FunctionalGroupTable: public ConfigLoader, private boost::noncopyable
	{
		private:
			//boost::property_tree::ptree params;
			std::map<std::string, FunctionalGroup> functionalgroups;

		protected:
			FunctionalGroupTable() {}

		public:
			
			static FunctionalGroupTable& Instance();

			void load(const std::string& filename = "../config/commongroup.xml");

			FunctionalGroup getFunctionalGroupBySymbol(const std::string&) const;

	};
}



#endif   /* ----- #ifndef GAG_FUNCTIONALGROUPTABLE_H ----- */
