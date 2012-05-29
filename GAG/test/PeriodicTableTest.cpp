/*
 * =====================================================================================
 *
 *       Filename:  PeriodicTableTest.cpp
 *
 *    Description:  Test file for class Periodic
 *
 *        Version:  1.0
 *        Created:  04/26/2012  4:17:25 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */
#include <iostream>

#include <GAGPL/CHEMISTRY/PeriodicTable.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
int main ( int argc, char *argv[] )
{

	using namespace std;
	using namespace gag;

	PeriodicTable& pt = PeriodicTable::Instance();
	pt.load();
	//cout << pt.get("parameters.PeakParameters.PeakBackgroundRatio", 0.0) << '\n';
	//cout << pt.get<std::string>("parameters.PeakParameters.WritePeaksToTextFile", "True")
	cout << "We are going to add Ca" << endl;

	Element e = pt.getElementBySymbol("Ca");
	
	cout << e.name << ' ' << e.isotopes.size() << endl;

	for (size_t i = 0; i < e.isotopes.size(); i++)
	{
		cout << e.isotopes[i].mass << ": " << e.isotopes[i].abundance << endl;
	}

	return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
