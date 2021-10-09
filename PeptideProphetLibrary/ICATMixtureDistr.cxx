
#include "ICATMixtureDistr.h"


/*

Program       : ICATMixtureDistr for PeptideProphet                                                       
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


ICATMixtureDistr::ICATMixtureDistr(int charge, char* name, char* tag) : DiscreteMixtureDistr(charge, 2, name, tag) {
  char* bindefs[] = {"icat incompatible", "icat compatible"};
  maxdiff_ = 0.005;
  negOnly_ = True;
  DiscreteMixtureDistr::init(bindefs);
}

int ICATMixtureDistr::inttranslate(char* val) {
  if(icatCompatible(val)) {
    return 1;
  }
  return 0;
}

Boolean2 ICATMixtureDistr::icatCompatible(const char* pep) {
  Boolean2 light = False;
  Boolean2 heavy = False;
  for(int k = 0; k < strlen(pep)-1; k++) {
    if(pep[k] == 'C') {
      if(pep[k+1] == '*' || pep[k+1] == '#' || pep[k+1] == '@') {
        heavy = True;
      }
      else {
        light = True;
      }
    }
  }
  if(pep[strlen(pep)-1] == 'C') {
    light = True;
  }
  if(light && heavy) {
    return False;
  }
  if(! light && ! heavy) {
    return False;
  }
  return True;
}
