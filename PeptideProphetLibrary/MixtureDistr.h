#ifndef MIX_DISTR
#define MIX_DISTR


#include <string.h>
#include <iostream>
#undef abs
#include <string>

#include "Array.h"
#include "GammaDistribution.h"
#include "GaussianDistribution.h"
#include "DiscreteDistribution.h"
#include "Distribution.h"

/*

Program       : MixtureDistr for PeptideProphet                                                       
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


class MixtureDistr {

 public:

  MixtureDistr();
  MixtureDistr(int charge, char* name, char* tag);
  ~MixtureDistr();
  virtual float getPosProb(int index);
  virtual float getNegProb(int index);
  virtual void enter(int index, char* val);
  virtual Boolean2 update(Array<float>* probs);
  virtual char* getName();
  virtual char* getTag();
  virtual void printVal(int index);
  virtual void enter(int index, float val);
  virtual void enter(int index, int val);
  virtual void printDistr();
  int getNumVals();
  Boolean2 isValue(int index, int val);
  virtual void setPosDistr(MixtureDistr* distr);
  virtual Distribution* getPosDistr();
  virtual Boolean2 negOnly();
  virtual void writeDistr(FILE* fout);


  protected:

  virtual int inttranslate(char* val);
  virtual float floattranslate(char* val);
  void initializeDistr(int charge, char* name, char* tag);
  Distribution* posdistr_;
  Distribution* negdistr_;
  char* name_;
  char* tag_; // for data entry
  Array<int>* intvals_;
  Array<float>* floatvals_;
  int charge_;
  float* pospriors_;
  float* negpriors_;
  float minval_;
  Boolean2 negOnly_;
};

#endif
