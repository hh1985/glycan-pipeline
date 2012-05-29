/*
 * =====================================================================================
 *
 *       Filename:  Elements.h
 *
 *    Description:  This class stores the basic information of elements including
 *    							Symbol, mass, Alias name and so on. This information can come from 
 *    							config files storing in a separate element.conf file.
 *
 *        Version:  1.0
 *        Created:  4/22/2012 3:06:04 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu 
 *   Organization:  
 *
 * =====================================================================================
 */


#ifndef  GAG_ELEMENT_H_INC
#define  GAG_ELEMENT_H_INC

#include <string>
#include <vector>

namespace gag
{
  struct Isotope
  {
    double mass;
    float abundance;
  };
  struct Element 
  {
    std::string symbol;
    std::string name;
		int atomicity;
    std::vector<Isotope> isotopes;
  };
}

#endif   /* ----- #ifndef GAG_ELEMENT_H_INC  ----- */
