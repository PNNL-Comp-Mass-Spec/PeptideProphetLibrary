#include "DiscreteDistribution.h"


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

DiscreteDistribution::DiscreteDistribution() { }

DiscreteDistribution::~DiscreteDistribution() { 
  if(priors_ != NULL) {
    delete [] priors_;
  }
}


DiscreteDistribution::DiscreteDistribution(float maxdiff) : Distribution(maxdiff) { 
}

DiscreteDistribution::DiscreteDistribution(int numbins, float maxdiff) : Distribution(maxdiff) { 
  num_bins_ = numbins;

  tot_ = new float[num_bins_];
  for(int k = 0; k < num_bins_; k++) {
    tot_[k] = 1.0 / num_bins_;
  }
  totwt_ = 1.0;
  priors_ = NULL;
  num_priors_ = 0.0;
}

int DiscreteDistribution::getNumBins() { return num_bins_; }

float DiscreteDistribution::getProb(float val) {
  return -1;
}

void DiscreteDistribution::setPriors(float* priors, float numpriors) 
{
	if (priors_ != NULL)
	{
		delete [] priors_ ; 
	}
	priors_ = priors;
	num_priors_ = numpriors;
}

float DiscreteDistribution::getProb(int val) {
  return tot_[val] / totwt_;
}
void DiscreteDistribution::initUpdate(float* prior) {

  if(newtot_ != NULL) 
  {
	  delete [] newtot_ ; 
  }
  newtot_ = new float[num_bins_];
  
  for(int k = 0; k < num_bins_; k++) {
    if(priors_ == NULL) {
      newtot_[k] = 0.0;
    }
    else {
      newtot_[k] = priors_[k] * num_priors_;
    }
  }
  if(priors_ == NULL) {
    newtotwt_ = 0.0;
  }
  else {
    newtotwt_ = num_priors_;
  }
}

void DiscreteDistribution::init(float* priors) {
    if(newtot_ != NULL) 
	{
		delete [] newtot_ ; 
    }
    newtot_ = new float[num_bins_];
    newtot_[0] = 0.0;
    newtotwt_ = 0.0;
    for(int k = 0; k < num_bins_; k++) {
      newtot_[k] = 0.0;
    }
    set_ = False;
}

void DiscreteDistribution::addVal(float wt, float val) {
}

void DiscreteDistribution::addVal(float wt, int val) {
  newtot_[val] += wt;
  newtotwt_ += wt;
}

Boolean2 DiscreteDistribution::update() {

  Boolean2 output = False;
  float delta;
  for(int k = 0; k < num_bins_; k++) {
    if(! set_ || (totwt_ == 0 && newtotwt_ > 0) || (totwt_ > 0 && newtotwt_ > 0)) {
      delta = (tot_[k]/totwt_) - (newtot_[k]/newtotwt_);
      if(abs(delta) > maxdiff_) {
	output = True;
      }
    }

  }
  if(output) {
    for(int k = 0; k < num_bins_; k++) {
      tot_[k] = newtot_[k];
    }
    totwt_ = newtotwt_;
    set_ = True;
  }
  if(! output) {
    //std::cout << "*** no change in discrete distribution..." << std::endl;
  }

  return output;
}
void DiscreteDistribution::writeDistr(FILE* fout) { 
}

void DiscreteDistribution::printDistr() {
}
