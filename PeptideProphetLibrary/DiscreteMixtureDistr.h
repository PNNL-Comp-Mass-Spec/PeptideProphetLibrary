#ifndef DISCR_MIX_DISTR_H
#define DISCR_MIX_DISTR_H

#include "MixtureDistr.h"

/*

Program       : DiscreteMixtureDistr for PeptideProphet                                                       
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


class DiscreteMixtureDistr : public MixtureDistr {
 public:

  DiscreteMixtureDistr(int charge, int numbins, char* name, char* tag);
 
  virtual void init(char** bindefs);
  virtual void printDistr();
  virtual void printPosDistribution();
  virtual void printNegDistribution();
  Boolean2 update(Array<float>* probs);
  void writeDistr(FILE* fout);

 protected:

  void initializeBinDefs(char** bindefs);
  float maxdiff_;
  Array<char*>* bindefs_;
  int numbins_;
  float* priors_;
  float numpos_priors_;
  float numneg_priors_;



};

#endif
