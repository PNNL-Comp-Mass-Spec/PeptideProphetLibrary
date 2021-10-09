#include "ExtremeValueDistribution.h"

/*

Program       : ExtremeValueDistribution for PeptideProphet                                                       
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

ExtremeValueDistribution::ExtremeValueDistribution(float maxdiff) : ContinuousDistribution(maxdiff) {

}

void ExtremeValueDistribution::init(float* prior) {
  ContinuousDistribution::init(NULL);
  mean_ = prior[0];
  stdev_ = prior[1];
  computeMoments(mean_, stdev_);
}

void ExtremeValueDistribution::initUpdate(float* prior) {
  if(newtot_ != NULL) {
      delete  [] newtot_ ; 
  }
  newtot_ = new float[1];
  newtot_[0] = 0.0;
  newtotsq_ = 0.0;
  newtotwt_ = 0.0;
}

void ExtremeValueDistribution::addVal(float wt, float val) {
  newtot_[0] += wt * val;
  newtotsq_ += wt * val * val;
  newtotwt_ += wt;
}

void ExtremeValueDistribution::addVal(float wt, int val) { }

Boolean2 ExtremeValueDistribution::update() {
  float newmean = newtot_[0]/newtotwt_;
  float newstdev = (newtotsq_ / newtotwt_) - newmean * newmean;
  newstdev = sqrt(newstdev);
  Boolean2 output = False;
  float delta = newmean - mean_;

  if(abs(delta) >= maxdiff_) {
    output = True;
  }

  delta = newstdev - stdev_;

  if(abs(delta) >= maxdiff_) {
    output = True;
  }
  if(output) {
    mean_ = newmean;
    stdev_ = newstdev;
    computeMoments(mean_, stdev_);
  }

  return output;
}

float ExtremeValueDistribution::getProb(int val) { return 0.0; }

float ExtremeValueDistribution::getExtremeValueProb(float val, float mu, float beta) { 
  if(beta == 0.0) {
    if(val == mu) {
      return 1.0;
    }
    else {
      return 0.0;
    }
  }
  float negexponent = exp(-(val - mu)/beta);
  return negexponent * exp(-negexponent) / beta;
}

float ExtremeValueDistribution::getProb(float val) {
  return getExtremeValueProb(val, mu_, beta_);
}

void ExtremeValueDistribution::computeMoments(float mean, float std) {
  beta_ = 0.779697  * std;
  mu_ = mean - 0.5772 * beta_;
}


float ExtremeValueDistribution::cumulative(float val, float mu, float beta) {
  if(beta == 0.0)
    return 1.0;
  float negexponent = exp(-(val - mu)/beta);
  return exp(-negexponent);

}

float ExtremeValueDistribution::extremeValueSlice(float num, float left_val, float right_val, float mu, float beta) {
  float diff = cumulative(right_val, mu, beta) - cumulative(left_val, mu, beta);
  return num * diff;
}

float ExtremeValueDistribution::slice(float num, float left_val, float right_val) {
  return extremeValueSlice(num, left_val, right_val, mu_, beta_);
}

void ExtremeValueDistribution::printDistr() {
  printf("(evd mean: %0.2f, stdev: %0.2f, mu: %0.2f, beta: %0.2f)\n", mean_, stdev_, mu_, beta_);

}

void ExtremeValueDistribution::writeDistr(FILE* fout) {
  fprintf(fout, "(evd mean: %0.2f, stdev: %0.2f, mu: %0.2f, beta: %0.2f)\n", mean_, stdev_, mu_, beta_);

}
