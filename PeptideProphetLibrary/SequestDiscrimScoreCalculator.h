#ifndef SEQ_DISCR_CALC_H
#define SEQ_DISCR_CALC_H

#include <stdio.h>
#include <iostream>
#include <fstream>

#include "ScoreCalculator.h"
#include "SequestDiscrimFunction.h"
#include "AbbrevSequestDiscrimFunction.h"
#include "ToftofAbbrevSequestDiscrimFunction.h"

/*

Program       : SequestScoreDiscrimCalculator for discr_calc of PeptideProphet 
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


class SequestDiscrimScoreCalculator : public ScoreCalculator {

 public:

  SequestDiscrimScoreCalculator();
  SequestDiscrimScoreCalculator(std::vector<SequestResult> &results, Boolean2 exclude_deltastars, Boolean2 windows, Boolean2 massd, Boolean2 mod_delt, Boolean2 maldi);  
  SequestDiscrimScoreCalculator(std::vector<SequestResult> &results, Boolean2 exclude_deltastars, Boolean2 windows, Boolean2 massd, Boolean2 mod_delt, Boolean2 maldi, char* enz_spec);
  void processData(std::vector<SequestResult> &results) ;

  void Abort() { abort = true; }
 private:

  void init(std::vector<SequestResult> &results, Boolean2 exclude_deltastars, Boolean2 windows, Boolean2 massd);
  SequestDiscrimFunction** discr_funcs_; // one for each charge
  Boolean2 exclude_deltastars_;
  Boolean2 set_deltastars_zero_;
  Boolean2 windows_;
  Boolean2 modify_deltas_;
  bool abort;
};


#endif
