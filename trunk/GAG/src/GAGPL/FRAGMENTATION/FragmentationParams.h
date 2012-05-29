/*
 * =====================================================================================
 *
 *       Filename:  FragmentationParams.h
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

#ifndef GAG_FRAGMENTATIONPARAMS_H
#define GAG_FRAGMENTATIONPARAMS_H

#include <string>
#include <map>
#include <vector>

namespace gag
{
	//typedef std::pair<int, Composition> CompositionShift;
	struct FragmentationParams
	{
			std::string type;
			std::string cleavage_shift;
			std::map<std::string, std::vector<std::string> > dis_shift;			
	};
}

#endif /* GAG_FRAGMENTATIONPARAMS_H */