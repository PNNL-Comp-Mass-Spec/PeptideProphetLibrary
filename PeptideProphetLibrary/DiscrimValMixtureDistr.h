#ifndef DISCRVAL_MIX_H
#define DISCRVAL_MIX_H

#include "MixtureDistr.h"
#include "NTTMixtureDistr.h"
#include "ContinuousDistribution.h"
#include "DecayContinuousMultimixtureDistr.h"

/*

Program       : DiscrimValMixtureDistr for PeptideProphet                                                       
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



class DiscrimValMixtureDistr : public MixtureDistr {

 public:

  DiscrimValMixtureDistr(int charge, char* name, char* tag, Boolean2 gamma, Boolean2 maldi, Boolean2 qtof);
  DiscrimValMixtureDistr() { }

  virtual Boolean2 initializeNegDistribution(NTTMixtureDistr* nttdistr);

  float getPosProb(int index);
  float getNegProb(int index);

  void writePosDistribution();
  void enter(int index, float val);
  void printDistr();
  int slice(float left_val, float right_val);
  float posSlice(float left_val, float right_val);
  float negSlice(float left_val, float right_val);
  float getRightCumulativeNegProb(int index, float right_val);
  float getRightCumulativeNegProb(float total, int index, float right_val);
  void writeDistr(FILE* fout);
  void reset();
  virtual Boolean2 noDistr();
  Boolean2 update(Array<float>* probs);
  void resetTot();
  Boolean2 decayMultimixture();
  float posSliceWithNTT(float left_val, float right_val, int ntt);

  float getnegmean() ;//Xiuxia added

 protected:

  float* copy(float* init, int num);
  float maxdiff_;
  Boolean2 gamma_;
  float min_dataval_;
  float minval_; // for gamma distribution, Xiuxia, {-2.0, -5.0, -5.0} for charge = 1, 2, 3
  void setNegativeDistr(float mean, float stdev, float zero);
  void setPositiveDistr(float mean, float stdev);
  float negmean_;
  int MIN_NUM_PSEUDOS_;
  int ZERO_SET_;
  int NUM_DEVS_;
  Boolean2 USE_TR_NEG_DISTR_;
  float* posinit_;
  float* neginit_;
  Boolean2 maldi_;
  Boolean2 all_negs_;
  Boolean2 qtof_;
}; 



#endif
