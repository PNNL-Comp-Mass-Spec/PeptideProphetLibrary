#ifndef MULT_ENZ_DIG_H
#define MULT_ENZ_DIG_H

#include "EnzymeDigestion.h"
#include "Array.h"

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


class MultipleEnzymeDigestion : public EnzymeDigestion {

 public:

  MultipleEnzymeDigestion();
  int numMissedCleavages(std::string peptide);
  int numCompatibleTermini(std::string peptide);
  void addEnzyme(EnzymeDigestion* enz);

 protected:
  Array<EnzymeDigestion*>* enzs_;

};


#endif