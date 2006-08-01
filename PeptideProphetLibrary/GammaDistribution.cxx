#include "GammaDistribution.h"

/*

Program       : GammaDistribution for PeptideProphet                                                       
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


GammaDistribution::GammaDistribution(float maxdiff) : ContinuousDistribution(maxdiff) {
  minval_ = NULL;
  mean_ = 0.0;
  stdev_ = 0.0;
  zero_ = 0.0;
  alpha_ = 0.0;
  beta_ = 0.0;
}

GammaDistribution::~GammaDistribution() {
  if(minval_ != NULL) {
    delete [] minval_;
  }
}

void GammaDistribution::init(float* prior) {
  ContinuousDistribution::init(NULL);
  if(prior != NULL) {
    mean_ = prior[0];
    stdev_ = prior[1];
    zero_ = prior[2];
    computeAlphaBeta();
  }
}

float GammaDistribution::getZero() { return zero_; }

Boolean2 GammaDistribution::zeroProb(float val) {
  return (getProb(val) < 0.01 && val < alpha_ * (beta_ - 1.0) + 
zero_);
}


void GammaDistribution::computeAlphaBeta() {
  float maxbeta = 125.0;

  if(mean_ == 0) {
    alpha_ = 999.0;
  }
  else {
    alpha_ = (stdev_ - mean_ * mean_) / mean_;
  }
  if(alpha_ == 0.0) {
    beta_ = 999.0;
  }
  else {
    beta_ = mean_ / alpha_;
  }
  if(beta_ > maxbeta) {
    beta_ = maxbeta;
    alpha_ = mean_ / beta_;
  }
}


void GammaDistribution::initUpdate(float* prior) { //Xiuxia, prior is not used here
  if(newtot_ == NULL) {
    newtot_ = new float[1];
  }
  newtot_[0] = 0.0;
  newtotsq_ = 0.0;
  newtotwt_ = 0.0;
}

void GammaDistribution::addVal(float wt, float val) {
  if(val < zero_ || ! aboveMin(val)) {
    return;
  }
  float adjval = val - zero_;
  newtot_[0] += wt * adjval;
  newtotsq_ += wt * adjval * adjval;
  newtotwt_ += wt;
}

void GammaDistribution::addVal(float wt, int val) { }

Boolean2 GammaDistribution::update() {
	//Xiuxia, update the mean_, stdev_, alpha, beta

	if(newtotwt_ == 0.0)
		return False;
	float newmean = newtot_[0]/newtotwt_;
	float newstdev = newtotsq_/newtotwt_;

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
		computeAlphaBeta();
	}
	return output;
}

float GammaDistribution::getProb(int val) { return 0.0; }

float GammaDistribution::getGammaProb(float val, float alpha, float beta, float zero) {
  double prob = 0;
  double gamm;
  //double helper;
  float value = val - zero;
  double loggamma;
  double first;
  if(beta <= 0) {
    return 0.0;
  }
  if(value > 0 && alpha > 0 && aboveMin(val)) {
    loggamma = gammln(beta);
	
	//Xiuxia, commented out
    //helper = ((beta - 1) * log(value)) - (beta * log(alpha)) - loggamma;
    //helper = exp(helper) * exp(- value / alpha);
    prob = pow((double)value, (double)(beta - 1)) * exp(- (double)value / (double)alpha);

    gamm = exp(loggamma);//Xiuxia, how can gamm be zero from this exponential function?
    if(gamm == 0.0) {
      return 1.0;
    }
    first = pow((double)alpha, (double) beta) * gamm;
    if(first == 0) {
      return 0.0;
    }
    prob /= first;
  }
  return (float)prob;
}


// add real thing here...
float GammaDistribution::getProb(float val) { 
  return getGammaProb(val, alpha_, beta_, zero_);
}



float GammaDistribution::gammln(float xx) {
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
        int j;

        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
}

void GammaDistribution::setDistrMinval(float val) {
	if (minval_ == NULL)
	  minval_ = new float[1];
  minval_[0] = val;	//Xiuxia, probability will be zero below the minval_
}

Boolean2 GammaDistribution::aboveMin(float val) {
  return (minval_ == NULL || val >= minval_[0]);
}

float GammaDistribution::slice(float num, float left_val, float right_val) {
  return gammaSlice(num, left_val, right_val, mean_, stdev_);
}

float GammaDistribution::gammaSlice(float num, float left_val, float right_val, float m1, float m2) {
  if(m2 == 0) {
    if(m2 >= left_val && m2 <= right_val) {
      return num;
    }
    else {
      return 0.0;
    }
  }
  int num_windows = 20;
  float window_width = (right_val - left_val) / num_windows;
  float next;
  float tot = 0.0;
  for(int k = 0; k <= num_windows; k++) {
    next = getProb((k * window_width) + left_val);
    if(k == 0 || k == num_windows) {
      tot += next / 2;
    }
    else {
      tot += next;
    }
  }
  return (tot * window_width * num);

}


void GammaDistribution::printDistr() {
  printf("(gamma m1: %0.2f, m2: %0.2f, alpha: %0.2f, beta: %0.2f, zero: %0.2f)\n", mean_, stdev_, alpha_, beta_, zero_);
}
void GammaDistribution::writeDistr(FILE* fout) {
  fprintf(fout, "(gamma m1: %0.2f, m2: %0.2f, alpha: %0.2f, beta: %0.2f, zero: %0.2f)\n", mean_, stdev_, alpha_, beta_, zero_);
}
