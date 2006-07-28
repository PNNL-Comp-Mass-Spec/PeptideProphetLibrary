#ifndef CONT_MULT_MIX_H
#define CONT_MULT_MIX_H

#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "ContinuousDistribution.h"
#include "GaussianDistribution.h"
#include "GammaDistribution.h"
#include "Array.h"

/*

Program       : ContinuousMultimixtureDistr for PeptideProphet                                                       
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

class ContinuousMultimixtureDistr : public ContinuousDistribution {


 public:

  ContinuousMultimixtureDistr(float maxdiff);
  ~ContinuousMultimixtureDistr();
  void initialize();
  void addComponent(float* settings, char* distr, char* def);
  float getProb(float val);
  float getProb(int val);
  void init(float* priors); // set index_ to 0
  void addVal(float wt, float val);
  void addVal(float wt, int val);
  int getNumDistributions();
  Boolean2 update();
  void initUpdate(float* prior);
  void printDistr();
  void writeDistr(FILE* fout);
  void setMinVal(float min); 
  float slice(float num, float left_val, float right_val);
  float slice(float left_val, float right_val);
  virtual int presetZvalue(int index);
  float getMean();
  float getStdev();
  float getMixtureDistrProb(int k, float val);
  virtual void removeViolatingDistrs();
  virtual Boolean2 violation(int leftdistr, int rightdistr);

 protected:

  Array<ContinuousDistribution*>* distrs_;
  //Array<GaussianDistribution*>* distrs_;
  float* priors_;
  Array<float*>* zvals_;
  int index_;
  Array<float>* vals_; // keep copy of data
  Array<char*>* defs_; // describe each distr
  Array<float>* wts_;
  float totalwts_;
  Array<char*>* distrtypes_;
};

#endif
