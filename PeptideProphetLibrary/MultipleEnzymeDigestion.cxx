#include "MultipleEnzymeDigestion.h"

/*

Program       : MultipleEnzymeDigestion for PeptideProphet
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

MultipleEnzymeDigestion::MultipleEnzymeDigestion(): EnzymeDigestion("", "", "", 1, 0) {
  enzs_ = new Array<EnzymeDigestion*>;
}

void MultipleEnzymeDigestion::addEnzyme(EnzymeDigestion* enz) {
  enzs_->insertAtEnd(enz);
}

int MultipleEnzymeDigestion::numMissedCleavages(std::string pep) {
  // minimum for all enz's
  int min = 10;
  int next;
  for(int k = 0; k < enzs_->length(); k++) {
    next = (*enzs_)[k]->numMissedCleavages(pep);
    if(next < min)
      min = next;
  }
  return min;
}

int MultipleEnzymeDigestion::numCompatibleTermini(std::string peptide)
{
    const char *pep = peptide.c_str() ;
    // max for all enzs
    int max = 0;
    int next;
    for(int k = 0; k < enzs_->length(); k++)
    {
        next = (*enzs_)[k]->numCompatibleTermini(pep);
        if(next > max)
            max = next;
    }
    return max;
}


