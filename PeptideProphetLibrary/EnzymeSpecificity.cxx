#include "EnzymeSpecificity.h"
#include <string>
#include <stdlib.h>
#include <ctype.h>
#using <mscorlib.dll>

using namespace std ;
using namespace System ;
using namespace System::Collections ;


/*

Program       : EnzymeSpecificity for PeptideProphet                                                       
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 

Copyright (C) 2003 Andrew Keller

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Andrew Keller
Insitute for Systems Biology
1441 North 34th St. 
Seattle, WA  98103  USA
akeller@systemsbiology.org

Institute for Systems Biology, hereby disclaims all copyright interest 
in PeptideProphet written by Andrew Keller

*/

EnzymeSpecificity::EnzymeSpecificity() { }


// factory for producing enzyme digestions
// register all new enzyme digestions here
EnzymeDigestion* EnzymeSpecificity::getEnzymeDigestion(char* enz) {
	//Xiuxia, 07/20/2006, convert enz to lowercase so that the original enz can be any case or a mixture of lower and upper case
	char *p ;
	char enz_lowercase[32] ;
	int i = 0 ;

	for (p = enz; p < enz + strlen(enz); p++)
	{
		if (isupper(*p))
		{
			enz_lowercase[i] = _tolower(*p) ;
		}
		else
		{
			enz_lowercase[i] = enz[i] ;
		}
		i++ ;
	}
	enz_lowercase[strlen(enz_lowercase)] = '\0' ;

  if(enz_lowercase == NULL) // default
    return new TrypticEnzymeDigestion();
  if(strcmp(enz_lowercase, "tryptic") == 0)
    return new TrypticEnzymeDigestion();
  if(strcmp(enz_lowercase, "gluC") == 0)
    return new GluCEnzymeDigestion();
  if(strcmp(enz_lowercase, "gluC_bicarb") == 0)
    return new GluC_bicarbEnzymeDigestion();
  if(strcmp(enz_lowercase, "chymotryptic") == 0)
    return new ChymotrypticEnzymeDigestion();
  if(strcmp(enz_lowercase, "elastase") == 0)
    return new Elastase();
  if(strcmp(enz_lowercase, "nonspecific") == 0) 
    return new NonspecificEnzymeDigestion();
  if(strcmp(enz_lowercase, "tca") == 0)
    return new TrypChymAspnEnzymeDigestion();

  // enzyme not registered
  return NULL;
}
