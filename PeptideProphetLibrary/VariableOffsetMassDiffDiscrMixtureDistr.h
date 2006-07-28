#ifndef VAR_OR_MASSD_H
#define VAR_OR_MASSD_H

#include <assert.h>

#include "MassDifferenceDiscrMixtureDistr.h"




/*

Program       : VariableOffsetMassDiffDiscrMixtureDistr for PeptideProphet                                                       
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 

/*

Program       : ChymotrypticEnzymeDigestion for PeptideProphet                                                       
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

class VariableOffsetMassDiffDiscrMixtureDistr : public MassDifferenceDiscrMixtureDistr {

 public:

  VariableOffsetMassDiffDiscrMixtureDistr(int charge, char* name, char* tag, float range, float window, float orig);
  int getIntegralValue(float val);
  float getMode(float window, Array<float>* probs);
  void enter(int index, char* val);
  Boolean2 update(Array<float>* probs);
  void writeDistr(FILE* fout);

 protected:
  float offset_;
  float offset_init_;
  Array<float>* vals_;
  int update_ctr_;
  int min_ctr_;
  Boolean2 offset_set_;


};


#endif
