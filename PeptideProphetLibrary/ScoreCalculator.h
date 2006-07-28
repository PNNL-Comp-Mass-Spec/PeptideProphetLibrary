#ifndef SCORE_CALC_H
#define SCORE_CALC_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "EnzymeDigestion.h"
#include "EnzymeSpecificity.h"
#include "SequestResult.h"
#include "sysdepend.h"


/*

Program       : ScoreCalculator for discr_calc of PeptideProphet 
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


class ScoreCalculator {

 public:
  
  ScoreCalculator();
  ScoreCalculator(char* enz_spec);
  virtual void processData(std::vector<SequestResult> &results) = 0;
  int numMissedCleavages(std::string peptide);
  int numCompatibleTermini(std::string peptide);
  std::string strip(std::string pep, Boolean2 remove_mods);
  char* parseSpectrum(char* spec);
  Boolean2 isNumber(char c);


 protected:

  Boolean2 maldi_set_;
  Boolean2 maldi_;
  char* enzyme_;
  EnzymeDigestion* enzyme_digestion_;

};


#endif
