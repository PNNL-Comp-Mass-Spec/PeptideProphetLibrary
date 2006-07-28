#include "MascotDiscrimValMixtureDistr.h"

/*

Program       : MascotDiscrimValMixtureDistr for PeptideProphet                                                       
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

MascotDiscrimValMixtureDistr::MascotDiscrimValMixtureDistr(int charge, char* name, char* tag, Boolean2 maldi, Boolean2 qtof) {  
  initializeDistr(charge, name, tag);
  all_negs_ = False;
  maldi_ = maldi;
  maxdiff_ = 0.002;
  gamma_ = False;
  qtof_ = qtof;

  floatvals_ = new Array<float>();
  posdistr_ = new GaussianDistribution(maxdiff_);
  //posdistr_ = new ExtremeValueDistribution(maxdiff_);
  negdistr_ = new ExtremeValueDistribution(maxdiff_);
  
  // evd with discriminant score (Mascot + delta)
  //float singlyposprior[] = {2.0, 0.4};
  //float doublyposprior[] = {3.754, 2.3};
  //float triplyposprior[] = {4.117, 2.9};
  //float singlynegprior[] = {-2.0, 0.9};
  //float doublynegprior[] = {-0.377, 0.64};
  //float triplynegprior[] = {-0.224, 0.67};

  float singlyposmaldiprior[] = {5.0, 0.6};
  float singlynegmaldiprior[] = {1.0, 1.5};

  // single score
  float singlyposprior[] = {-0.191, 0.833};
  float doublyposprior[] = {1.109, 1.767};
  float triplyposprior[] = {0.411, 1.545};
  float singlynegprior[] = {-1.748, 0.539};
  float doublynegprior[] = {-2.468, 0.374};
  float triplynegprior[] = {-2.179, 0.363};
 

  negmean_ = 0.0;
  MIN_NUM_PSEUDOS_ = 50;
  ZERO_SET_ = 100; // ?
  NUM_DEVS_ = 6;
  USE_TR_NEG_DISTR_ = False;
  posinit_ = NULL;
  neginit_ = NULL;

  if(charge == 0) {
    if(maldi_) {
      posinit_ = copy(singlyposmaldiprior, 2);
    }
    else {
      posinit_ = copy(singlyposprior, 2);
    }

    if(maldi_) {
      neginit_ = copy(singlynegmaldiprior, 2);
    }
    else {
      neginit_ = copy(singlynegprior, 2);
    }

  }
  else if(charge == 1) {
    posinit_ = copy(doublyposprior, 2);
    neginit_ = copy(doublynegprior, 2);
  }
  else if(charge == 2) {
    posinit_ = copy(triplyposprior, 2);
    neginit_ = copy(triplynegprior, 2);
  }
  reset();

}


float MascotDiscrimValMixtureDistr::getPosProb(int index) { 
  if(all_negs_ || (*floatvals_)[index] < negmean_) {
    return 0.0;
  }
  return MixtureDistr::getPosProb(index);
}

float MascotDiscrimValMixtureDistr::getNegProb(int index) {
  if(all_negs_ || (*floatvals_)[index] < negmean_) {
    return 1.0;
  }
  return MixtureDistr::getNegProb(index);
}

Boolean2 MascotDiscrimValMixtureDistr::initializeNegDistribution(NTTMixtureDistr* nttdistr) {
  assert(nttdistr->getNumVals() == getNumVals());
  float mean = 0.0;
  float stdev = 0.0;
  int tot = 0;
  float totsq = 0.0;
  float zero;

  // evd for discriminatn score
  //float posmean[] = { 2.0, 4.102, 4.563 };
  //float posstdev[] = { 0.4, 1.64, 1.84 };

  // for single score guassian
  float posmean[] = { -.191, 1.109, 0.411};
  float posstdev[] = { 0.833, 1.77, 1.55 };


  float MAX_SINGLY_NEGMEAN = 1.0;

  //float negmean_num_stds[] = {-1.0, 0.5, 0.5}; // by charge
  float negmean_num_stds[] = {-0.1, 1.0, 1.0}; // by charge

  float min_singly_fval = -2;

  for(int k = 0; k < getNumVals(); k++) {
    if(nttdistr->isValue(k, 0)) {
      mean += (*floatvals_)[k];
      totsq += (*floatvals_)[k] * (*floatvals_)[k];
      tot++;
    }
  }

  if(tot < MIN_NUM_PSEUDOS_) {
 
    tot = 0; // restart
    mean = 0.0;
    totsq = 0.0;
    for(int k = 0; k < getNumVals(); k++) {
      if((*floatvals_)[k] < posmean[charge_] - posstdev[charge_]) {
	mean += (*floatvals_)[k];
	totsq += ((*floatvals_)[k]) * ((*floatvals_)[k]);
	tot++;
      }
    } // next

    if(tot > 0) {
      mean /= tot;
      stdev = totsq / tot - mean * mean;
    }

    float* newsettings = new float[2];
    newsettings[0] = mean;
    newsettings[1] = stdev;
    newsettings[1] = sqrt(newsettings[1]);
    negdistr_->init(newsettings);
    delete [] newsettings;


    negmean_ = mean - negmean_num_stds[charge_] * stdev;
    if(negmean_ > MAX_SINGLY_NEGMEAN) {
      negmean_ = MAX_SINGLY_NEGMEAN;
    }

    USE_TR_NEG_DISTR_ = True;
    return False; // done
  } // if not enough pseudos
    
  mean /= tot;
  stdev = (totsq / tot) - mean * mean;
  stdev = sqrt(stdev);  // delete this later ?

  negmean_ = mean - negmean_num_stds[charge_] * stdev;

  // now recompute real mean and stdev
  mean = 0.0;
  totsq = 0.0;
  tot = 0;
    
  for(int k = 0; k < getNumVals(); k++) {
    if(nttdistr->isValue(k, 0)) {
      mean += (*floatvals_)[k];
      totsq += ((*floatvals_)[k]) * ((*floatvals_)[k]);
      tot++;
    }
  } // next
  
  if(tot > 0) {
    mean /= tot;
    stdev = (totsq / tot) - mean * mean;
    stdev = sqrt(stdev);
  }
  Boolean2 reset = False;
  float newposmean = posmean[charge_];
  float newposstdev = posstdev[charge_];;
  while(noDistr()) {
    newposmean += 0.5;
    newposstdev += 0.95;
    reset = True;
  }
  if(reset) {
    setPositiveDistr(newposmean, newposstdev);
  }

  setNegativeDistr(mean, stdev);
  return True;
}

Boolean2 MascotDiscrimValMixtureDistr::noDistr() {
  return ((ContinuousDistribution*)(posdistr_))->getMean() + ((ContinuousDistribution*)(posdistr_))->getStdev() < ((ContinuousDistribution*)(negdistr_))->getMean() + ((ContinuousDistribution*)(negdistr_))->getStdev(); 

}


void MascotDiscrimValMixtureDistr::setNegativeDistr(float mean, float stdev) {
  float* next = new float[2];
  next[0] = mean;
  next[1] = stdev;

  negdistr_->init(next);
  delete [] next;
}
