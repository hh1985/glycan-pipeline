/*
 * =====================================================================================
 *
 *       Filename:  FunctionalGroupTableTest.cpp
 *
 *    Description:  Test file for class FunctionalGroupTable.
 *
 *        Version:  1.0
 *        Created:  04/30/2012  1:39:20 AM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */
#include <iostream>

#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>

int main ( int argc, char *argv[] )
{

	using namespace std;
	using namespace gag;

	//PeriodicTable& pt = PeriodicTable::Instance();
	//pt.load();

	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();
	fgt.load();

	FunctionalGroup& fg = fgt.getFunctionalGroupBySymbol("OH");

	std::cout << "Symbol: " << fg.symbol << std::endl;
	std::cout << "Name: " << fg.name << std::endl;
	std::cout << "Composition: " << fg.getComposition().getCompositionString() << std::endl;

	for(std::vector<Composition>::const_iterator iter = fg.site_gps.begin(); iter != fg.site_gps.end(); iter++)
	{
		std::cout << "Site: " << (*iter).getMass() << std::endl;
	}

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
