#ifndef DECAY_MULT_H
#define DECAY_MULT_H

#include "ContinuousMultimixtureDistr.h"
#include "NTTMixtureDistr.h"

/*

Program       : DecayContinuousMultimixtureDistr for PeptideProphet                                                       
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 

/*

Program       : DecayContinuousMultimixtureDistr for PeptideProphet                                                       
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


class DecayContinuousMultimixtureDistr : public ContinuousMultimixtureDistr {

 public:

  DecayContinuousMultimixtureDistr(float maxdiff, Boolean2 qtof);
  ~DecayContinuousMultimixtureDistr();

  void initWithNTTs(NTTMixtureDistr* ntt);
  int presetZvalue(int index);
  float getMixtureProb(int index, float val);
  float getMixtureProbWithNTT(int ntt, int index);
  void addVal(float wt, float val);
  Boolean2 update();
  void printDistr();
  void writeDistr(FILE* fout);
  Boolean2 updateDistr(int k);
  void commence();
  Boolean2 oneProb(float val);
  void removeViolatingDistrs();
  Boolean2 violation(int leftdistr, int rightdistr);
  Boolean2 unmixed(int leftdistr, int rightdistr);
  void reset();
  float sliceWithNTT(float left_val, float right_val, int ntt);



 protected:

  NTTMixtureDistr* ntt_;
  Array<float*>* nttpriors_;
  float* newtotnttpriors_;
  float* newnttpriors_;
  Boolean2 qtof_;
  float min_total_wts_;
  Boolean2 equal_stdevs_;
  float* nttdistrtots_;
};

#endif
