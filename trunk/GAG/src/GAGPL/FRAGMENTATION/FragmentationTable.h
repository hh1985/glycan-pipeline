/*
 * =====================================================================================
 *
 *       Filename:  FragmentationTable.h
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


#ifndef	GAG_FRAGMENTATIONTABLE_H
#define GAG_FRAGMENTATIONTABLE_H

#include <GAGPL/MISC/ConfigLoader.h>
#include <GAGPL/CHEMISTRY/Composition.h>
#include <GAGPL/FRAGMENTATION/FragmentationParams.h>
#include <boost/noncopyable.hpp>

namespace gag
{
	typedef std::pair<int, Composition> CompositionShift;

	class FragmentationTable: public ConfigLoader, private boost::noncopyable
	{
		private:
			std::map<std::string, FragmentationParams> fragmentation_params;

		protected:
			FragmentationTable() {}

		public:
			
			static FragmentationTable& Instance();

			void load(const std::string& filename = "../config/fragmentation.xml");

			FragmentationParams getFragmentationParams(const std::string& cleavage_type) const;
			CompositionShift getCompositionShift(const std::string& str) const;
			CompositionShift getCleavageShift(const std::string& cleavage_type) const;
			std::vector<CompositionShift> getCleavageShift(const std::string& cleavage_type, const std::string& dis_type) const;			 

	};
}



#endif   /* ----- #ifndef GAG_FUNCTIONALGROUPTABLE_H ----- */
