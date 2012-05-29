/*
 * =====================================================================================
 *
 *       Filename:  PeriodicTable.h
 *
 *    Description:  Periodic Table loaded from xml file.
 *
 *        Version:  1.0
 *        Created:  4/24/2012 2:53:40 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com 
 *   Organization:  Boston University
 *
 * =====================================================================================
 */

#ifndef  GAG_PERIODICTABLE_H_INC
#define  GAG_PERIODICTABLE_H_INC

#include <GAGPL/CHEMISTRY/Element.h>
#include <GAGPL/MISC/ConfigLoader.h>
#include <boost/noncopyable.hpp>
#include <map>

namespace gag
{
  class PeriodicTable: public ConfigLoader, private boost::noncopyable
  {
    private:
    
			//boost::property_tree::ptree params;
			std::map<std::string, Element> elements;

		protected:

			PeriodicTable() {}
			
		public:

			static PeriodicTable& Instance();
			
			void load(const std::string& filename = "../config/SampleParameterFile.xml");
			
			// Currently we only consider about the situation of retrieving element by symbol.
			const Element getElementBySymbol(const std::string& symbol) const;

      inline Isotope getIsotope(const std::string& symbol, int isotope_order=0) const
      {
        return (*this).getElementBySymbol(symbol).isotopes.at(isotope_order);
      }
      size_t getNumIsotope(const std::string symbol) const
      {
        return (*this).getElementBySymbol(symbol).isotopes.size();
      }
  };
}




#endif   /* ----- #ifndef GAG_PERIODICTABLE_H_INC  ----- */
