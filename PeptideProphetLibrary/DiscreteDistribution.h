#ifndef DISCR_DISTR_H
#define DISCR_DISTR_H

#include "Distribution.h"

/*

Program       : DiscreteDistribution for PeptideProphet                                                       
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


class DiscreteDistribution : public Distribution {
  public:
  DiscreteDistribution();
  DiscreteDistribution(float maxdiff);
  DiscreteDistribution(int numbins, float maxdiff);
  ~DiscreteDistribution();

  float getProb(float val);
  float getProb(int val);
  void init(float* priors);
  void addVal(float wt, int val);
  void addVal(float wt, float val);
  Boolean2 update();
  void initUpdate(float* prior);
  void printDistr();
  virtual void setPriors(float* priors, float numpriors);
  void writeDistr(FILE* fout);
  int getNumBins();

  protected:

  int num_bins_;
  float* priors_;
  float num_priors_;

};




#endif 
