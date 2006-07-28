
#include "GlycMixtureDistr.h"
  

/*

Program       : GlycMixtureDistr for PeptideProphet                                                       
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

                
GlycMixtureDistr::GlycMixtureDistr(int charge, char* name, char* tag) : DiscreteMixtureDistr(charge, 2, name, tag) {
  char* bindefs[] = {"without glyc motif", "with glyc motif"};
  maxdiff_ = 0.005;
  negOnly_ = True;
  DiscreteMixtureDistr::init(bindefs);
}

int GlycMixtureDistr::inttranslate(char* val) {
  if(hasGlycMotif(val)) {
    return 1;
  }
  return 0;
}

Boolean2 GlycMixtureDistr::hasGlycMotif(char* pep) {
  for(int k = 0; k < strlen(pep)-2; k++) {
    if(pep[k] == 'N') {

      if(pep[k+1] != 'P' && (pep[k+2] == 'T' || pep[k+2] == 'S')) 
	return True;
      else if(k < strlen(pep) - 3 && isModification(pep[k+1]) && pep[k+2] != 'P' && (pep[k+3] == 'T' || pep[k+3] == 'S'))
	return True;
      else if(k < strlen(pep) - 3 && pep[k+1] != 'P' && isModification(pep[k+2]) && (pep[k+3] == 'T' || pep[k+3] == 'S')) 
	return True;
      else if(k < strlen(pep) - 4 && isModification(pep[k+1]) && pep[k+2] != 'P' && isModification(pep[k+3]) &&
	    (pep[k+4] == 'T' || pep[k+4] == 'S')) 
	return True;
    }
  } // next position in pep
  return False; // still here
}

Boolean2 GlycMixtureDistr::isModification(char c) {
  return (c == '*' || c == '#' || c == '@');
}
