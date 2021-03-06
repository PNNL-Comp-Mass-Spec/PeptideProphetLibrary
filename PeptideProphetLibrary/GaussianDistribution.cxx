#include "GaussianDistribution.h"

/*

Program       : GaussianDistribution for PeptideProphet                                                       
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

GaussianDistribution::GaussianDistribution(float maxdiff) : ContinuousDistribution(maxdiff) {
}

Boolean2 GaussianDistribution::oneProb(float val) {
  return val > mean_ + stdev_;
}

void GaussianDistribution::init(float* prior) {
  ContinuousDistribution::init(NULL);
  mean_ = prior[0];
  stdev_ = prior[1];
}

void GaussianDistribution::initUpdate(float* prior) {
  if(newtot_ == NULL) {
    newtot_ = new float[1];
  }
  newtot_[0] = 0.0;
  newtotsq_ = 0.0;
  newtotwt_ = 0.0;
}


void GaussianDistribution::addVal(float wt, float val) {
  newtot_[0] += wt * val;
  newtotsq_ += wt * val * val;
  newtotwt_ += wt;
}


void GaussianDistribution::addVal(float wt, int val) { }

Boolean2 GaussianDistribution::update() {
    float newmean = newtot_[0]/newtotwt_;
    float newstdev = (newtotsq_ / newtotwt_) - newmean * newmean;
    newstdev = sqrt(newstdev);

    Boolean2 output = False;

    float delta = newmean - mean_;
    if(abs(delta) >= maxdiff_) {
        output = True;
    }
    if(use_stdev_)
        delta = newstdev - stdev_;
    if(abs(delta) >= maxdiff_) {
        output = True;
    }    
    if(output) {
        mean_ = newmean;
        stdev_ = newstdev;
    }
    return output;
}


float GaussianDistribution::getProb(int val) { return 0.0; }

// add real thing here...
float GaussianDistribution::getGaussianProb(float val, float mean, float stdev) { 
  if(stdev == 0) {
    if(val == mean) {
      return 1.0;
    }
    else {
      return 0.0;
    }
  }
  float exponent = (-0.5) * (val - mean) * (val - mean) / (stdev * stdev);
  return ((exp (exponent)) / (stdev * sqrt(6.28318)));
}

// add real thing here...
float GaussianDistribution::getProb(float val) {
  return getGaussianProb(val, mean_, stdev_);
}

float GaussianDistribution::slice(float num, float left_val, float right_val) {
  return gaussianSlice(num, left_val, right_val, mean_, stdev_);
}
  
float GaussianDistribution::gaussianSlice(float num, float left_val, float right_val, float mean, float stdev) {
  if(stdev == 0) {
    if(mean >= left_val && mean <= right_val) {
      return num;
    }
    else {
      return 0.0;
    }
  }
  double left = N((left_val - mean_)/stdev_);
  double right = N((right_val - mean_)/stdev_);
  return num * (right - left);
}


// cumulative univariate normal distribution.
// This is a numerical approximation to the normal distribution.  
// See Abramowitz and Stegun: Handbook of Mathemathical functions
// for description.  The arguments to the functions are assumed 
// normalized to a (0,1 ) distribution. 

double GaussianDistribution::N(double z) {
    double b1 =  0.31938153; 
    double b2 = -0.356563782; 
    double b3 =  1.781477937;
    double b4 = -1.821255978;
    double b5 =  1.330274429; 
    double p  =  0.2316419; 
    double c2 =  0.3989423; 

    if (z >  6.0) { return 1.0; }; // this guards against overflow 
    if (z < -6.0) { return 0.0; };
    double a=fabs(z); 
    double t = 1.0/(1.0+a*p); 
    double b = c2*exp((-z)*(z/2.0)); 
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t; 
    n = 1.0-b*n; 
    if ( z < 0.0 ) n = 1.0 - n; 
    return n; 
}



void GaussianDistribution::printDistr() {
  printf("(gaussian mean: %0.2f, stdev: %0.2f)\n", mean_, stdev_);

}
void GaussianDistribution::writeDistr(FILE* fout) {
  fprintf(fout, "(gaussian mean: %0.2f, stdev: %0.2f)\n", mean_, stdev_);

}
