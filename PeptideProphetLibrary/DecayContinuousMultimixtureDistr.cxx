#include "DecayContinuousMultimixtureDistr.h"
#include <sstream>
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


DecayContinuousMultimixtureDistr::DecayContinuousMultimixtureDistr(float maxdiff, Boolean2 qtof) : ContinuousMultimixtureDistr(maxdiff) {
  nttpriors_ = new Array<float*>;
  qtof_ = qtof;
  min_total_wts_ = 5.0; 
  equal_stdevs_ = True;
  priors_ = NULL;
  nttdistrtots_ = NULL;
}

DecayContinuousMultimixtureDistr::~DecayContinuousMultimixtureDistr() {
  if(newnttpriors_ != NULL)
    delete[] newnttpriors_;
  if(newtotnttpriors_ != NULL)
    delete[] newtotnttpriors_;
}

float DecayContinuousMultimixtureDistr::getMixtureProb(int index, float val = -999.0) {

  if(index >= vals_->length())
    return getProb(val);
  return getMixtureProbWithNTT(ntt_->getNTTValue(index), index);
}

float DecayContinuousMultimixtureDistr::getMixtureProbWithNTT(int ntt, int index) {
  float val = (*vals_)[index];
  float prob = 0.0;
  for(int k = 0; k < distrs_->length(); k++)
    prob += (*distrs_)[k]->getProb(val) * ((*nttpriors_)[k])[ntt];
  return prob;
}

void DecayContinuousMultimixtureDistr::reset() {
  initWithNTTs(ntt_);
}

void DecayContinuousMultimixtureDistr::initWithNTTs(NTTMixtureDistr* ntt) {

  // should customize these for whether gaussian or gamma
  float nondecay_settings[] = {1.592, 1.199};
  float decay_settings[] = {0.743, 1.111};

  float qtof_nondecay_settings[] = {3.74, 1.4};
  float qtof_decay_settings[] = {2.2, 1.1};

  float* nondecay = new float[2];
  for(int k = 0; k < 2; k++)
    if(qtof_)
      nondecay[k] = qtof_nondecay_settings[k];
    else
      nondecay[k] = nondecay_settings[k];
  float* decay = new float[2];
  for(int k = 0; k < 2; k++)
    if(qtof_)
      decay[k] = qtof_decay_settings[k];
    else
      decay[k] = decay_settings[k];

  ntt_ = ntt;
  initialize();

  addComponent(nondecay, "gaussian", "nondecay");
  addComponent(decay, "gaussian", "decay");


  newtotnttpriors_ = new float[3];
  newnttpriors_ = new float[3];

  priors_ = new float[distrs_->length()];
  for(int k = 0; k < distrs_->length(); k++) {
    priors_[k] = 1.0 / distrs_->length();
    float* next = new float[3]; // one for each ntt value
    
    for(int ntt = 0; ntt < 2; ntt++)
      next[ntt] = 1.0 / distrs_->length(); // equal
    // special case
    if(k == 0)
      next[2] = 1.0;
    else
      next[2] = 0.0;
    
    nttpriors_->insertAtEnd(next);
  } // next distr

  if(equal_stdevs_)
    for(int k = 0; k < distrs_->length(); k++)
      (*distrs_)[k]->ignoreStdev();

}

Boolean2 DecayContinuousMultimixtureDistr::oneProb(float val) {
  return val > mean_ + stdev_;
}

int DecayContinuousMultimixtureDistr::presetZvalue(int index) {
  // if ntt 2, must be mixture 0
  if(ntt_->getNTTValue(index) == 2)
    return 0;
  return -1;
}

void DecayContinuousMultimixtureDistr::addVal(float wt, float val) {
  if(vals_->length() <= index_) {
    vals_->insertAtEnd(val);
    // set initial zvalues based upon distr's
    float totprob = 0.0;

    float* nextz = new float[distrs_->length()];
    for(int k = 0; k < distrs_->length(); k++) {
      nextz[k] = (*distrs_)[k]->getProb(val) * ((*nttpriors_)[k])[ntt_->getNTTValue(index_)];
      totprob += nextz[k];
    }
    if(totprob > 0)
      for(int k = 0; k < distrs_->length(); k++) 
	nextz[k] /= totprob;

    // now add it
    zvals_->insertAtEnd(nextz);
    ////std::cerr << "initial wt: " << wt << std::endl;
    wts_->insertAtEnd(wt);
    assert(index_ == vals_->length() - 1);
    
  } // if need to add
  assert((*vals_)[index_] == val);

  wts_->replace(index_, wt); // update

  // new: exclude ntt 0's from setting distributions.
  //float adj_wt = wt;
  //if(ntt_ != NULL && ntt_->getNTTValue(index_) == 0)
  //adj_wt = 0.0;
  for(int k = 0; k < distrs_->length(); k++)
    (*distrs_)[k]->addVal(wt * ((*zvals_)[index_])[k], val);

  assert(vals_->length() == zvals_->length());
  index_++;
}


Boolean2 DecayContinuousMultimixtureDistr::updateDistr(int k) {

  assert(vals_->length() == index_);

  Boolean2 output = False;
  if((*distrs_)[k]->update())
    output = True;

  return output;

}


// ready to allow other distributions into the mixture
void DecayContinuousMultimixtureDistr::commence() {
  totalwts_ = 0.0;

  int ntt;
  float nextprob;
  float totalprob = 0.0;
  for(int i = 0; i < vals_->length(); i++) {
    totalwts_ += (*wts_)[i];
    nextprob = 0.0;
    ntt = ntt_->getNTTValue(i);
    for(int k = 0; k < distrs_->length(); k++) {
      ((*zvals_)[i])[k] = (*distrs_)[k]->getProb((*vals_)[i]) * ((*nttpriors_)[k])[ntt];
      nextprob += ((*zvals_)[i])[k];
    }
  
    // now normalize
    if(nextprob > 0) 
      for(int k = 0; k < distrs_->length(); k++) 
	  ((*zvals_)[i])[k] /= nextprob;

  } // next data

  // now recompute priors
  for(int k = 0; k < 3; k++)
    newtotnttpriors_[k] = 0.0;

  for(int k = 0; k < distrs_->length(); k++) {
    priors_[k] = 0.0;
    for(int i = 0; i < vals_->length(); i++) {
      priors_[k] += ((*zvals_)[i])[k] * (*wts_)[i];
      totalprob += ((*zvals_)[i])[k] * (*wts_)[i];
    }
  }
  // now normalize priors
  if(totalprob > 0) 
    for(int k = 0; k < distrs_->length(); k++) 
      priors_[k] /= totalprob;

  for(int k = 0; k < distrs_->length(); k++) {
    for(int j = 0; j < 3; j++)
      newnttpriors_[j] = 0.0;

    for(int i = 0; i < vals_->length(); i++) {
      ntt = ntt_->getNTTValue(i);
      newnttpriors_[ntt] += ((*zvals_)[i])[k] * (*wts_)[i];
      newtotnttpriors_[ntt] += ((*zvals_)[i])[k] * (*wts_)[i]; // regardless of distribution
    }
    for(int j = 0; j < 3; j++)
      ((*nttpriors_)[k])[j] = newnttpriors_[j];

  } // next distribution
  // now normalize
  for(int j = 0; j < 3; j++)
    if(newtotnttpriors_[j] > 0)
      for(int k = 0; k < distrs_->length(); k++)
	((*nttpriors_)[k])[j] /= newtotnttpriors_[j];
  nttdistrtots_ = new float[distrs_->length()];

  set_ = True;
}

void DecayContinuousMultimixtureDistr::removeViolatingDistrs() {
  for(int k = 1; k < distrs_->length(); k++) {
    if((*distrs_)[k] != (*distrs_)[k-1] && priors_[k] > 0.0 && (violation(k, k-1) || unmixed(k, k-1))) {
      // set all nttpriors to 0.0
      for(int n = 0; n < 3; n++)
	((*nttpriors_)[k])[n] = 0.0;
    } // if violation
  } // next pair of adjacent distributions
}


Boolean2 DecayContinuousMultimixtureDistr::violation(int leftdistr, int rightdistr) {
  float num_stdevs = 0.0; //0.1;
  return (*distrs_)[leftdistr]->getMean() + num_stdevs * (*distrs_)[leftdistr]->getStdev() >= (*distrs_)[rightdistr]->getMean() - num_stdevs * (*distrs_)[rightdistr]->getStdev();
}

Boolean2 DecayContinuousMultimixtureDistr::unmixed(int leftdistr, int rightdistr) {
  float num_stdevs = 2.0; //0.1;
  return (*distrs_)[leftdistr]->getMean() + num_stdevs * (*distrs_)[leftdistr]->getStdev() <= (*distrs_)[rightdistr]->getMean() - num_stdevs * (*distrs_)[rightdistr]->getStdev();
}

Boolean2 DecayContinuousMultimixtureDistr::update() {
  if(! set_)
    commence();

  if(! set_)
    return updateDistr(0);


  if(vals_ == NULL || zvals_ == NULL || vals_->length() != zvals_->length()) {
    //std::cerr << "error in ContinuousMultimixtureDistr update...." << std::endl;
    //if(vals_ == NULL)
    //  //std::cerr << "null vals_" << std::endl;
    //else if(zvals_ == NULL)
    //  //std::cerr << "null zvals_" << std::endl;
    //else
      //std::cerr << vals_->length() << " vals_ <> " << zvals_->length() << " zvals_" << std::endl;
		std::stringstream str;
		str << "error in ContinuousMultimixtureDistr update...." << std::endl;
		if(vals_ == NULL)
			str << "null vals_" << std::endl;
		else if(zvals_ == NULL)
			str << "null zvals_" << std::endl;
		else
			str << vals_->length() << " vals_ <> " << zvals_->length() << " zvals_" << std::endl;

		throw new System::Exception(str.str().c_str());
    //exit(1);
  }

  assert(vals_->length() == index_);
  float totalprob = 0.0;
  Boolean2 output = False;


  for(int k = 0; k < distrs_->length(); k++)
    if((*distrs_)[k]->update())
      output = True;


  // here set old stdevs (just so cannot be the cause of update = true);
  float newstdev = 0.0;
  // ensure equal stdevs
  if(equal_stdevs_) {
    for(int k = 0; k < distrs_->length(); k++)
      newstdev += priors_[k] * (*distrs_)[k]->getStdev();


    // now set them
    for(int k = 0; k < distrs_->length(); k++)
      (*distrs_)[k]->setStdev(newstdev);


  } // if equal

  // now remove all violate distributions....
  ////std::cerr << "now removing distrs..." << std::endl;
  removeViolatingDistrs();

  if(totalwts_ < min_total_wts_) {
    for(int k = 0; k < distrs_->length(); k++) {
      if(k == 0)
	priors_[k] = 1.0;
      else
	priors_[k] = 0.0;
      for(int j = 0; j < 3; j++)
	if(k == 0)
	  ((*nttpriors_)[k])[j] = 1.0;
	else
	  ((*nttpriors_)[k])[j] = 0.0;
    } // next distr
  } // only use first distribution of mixture



  if(! output)
    return output;
  // now update zvals and priors


  totalwts_ = 0.0;

  float NUM_STDS = 2.5;
  int ntt;
  float nextprob;
  for(int i = 0; i < vals_->length(); i++) {
    totalwts_ += (*wts_)[i];
    nextprob = 0.0;
    ntt = ntt_->getNTTValue(i);
    for(int k = 0; k < distrs_->length(); k++) {
      /*
      if(ntt < 2 && (*vals_)[i] > (*distrs_)[0]->getMean() + NUM_STDS * (*distrs_)[0]->getStdev()) {
	if(k == 0)
	  ((*zvals_)[i])[k] = 1.0;
	else
	  ((*zvals_)[i])[k] = 0.0;
      }
      else if(ntt < 2 && (*vals_)[i] < (*distrs_)[distrs_->length()-1]->getMean() - NUM_STDS * (*distrs_)[distrs_->length()-1]->getStdev()) {
	if(k == distrs_->length()-1)
	  ((*zvals_)[i])[k] = 1.0;
	else
	  ((*zvals_)[i])[k] = 0.0;
      }
      else 
      */
	((*zvals_)[i])[k] = (*distrs_)[k]->getProb((*vals_)[i]) * ((*nttpriors_)[k])[ntt];
      nextprob += ((*zvals_)[i])[k];
    }
  
    // now normalize
    if(nextprob > 0) 
      for(int k = 0; k < distrs_->length(); k++) 
	  ((*zvals_)[i])[k] /= nextprob;

  } // next data

  // now recompute priors
  for(int k = 0; k < 3; k++)
    newtotnttpriors_[k] = 0.0;

  for(int k = 0; k < distrs_->length(); k++) {
    priors_[k] = 0.0;
    for(int i = 0; i < vals_->length(); i++) {
      priors_[k] += ((*zvals_)[i])[k] * (*wts_)[i];
      totalprob += ((*zvals_)[i])[k] * (*wts_)[i];
    }
  }
  // now normalize priors
  if(totalprob > 0) 
    for(int k = 0; k < distrs_->length(); k++) 
      priors_[k] /= totalprob;

  for(int k = 0; k < distrs_->length(); k++) {
    for(int j = 0; j < 3; j++)
      newnttpriors_[j] = 0.0;

    for(int i = 0; i < vals_->length(); i++) {
      ntt = ntt_->getNTTValue(i);
      newnttpriors_[ntt] += ((*zvals_)[i])[k] * (*wts_)[i];
      newtotnttpriors_[ntt] += ((*zvals_)[i])[k] * (*wts_)[i]; // regardless of distribution
    }
    for(int j = 0; j < 3; j++)
      ((*nttpriors_)[k])[j] = newnttpriors_[j];

  } // next distribution
  // now normalize
  for(int j = 0; j < 3; j++)
    if(newtotnttpriors_[j] > 0)
      for(int k = 0; k < distrs_->length(); k++)
	((*nttpriors_)[k])[j] /= newtotnttpriors_[j];

  ////std::cerr << " end priors: " << priors_[0] << ":" << priors_[1] << std::endl;
  ////std::cerr << "total wt: " << totalwts_ << std::endl;

  
  //if(! output)
  //return output;


  // now make adjustments based upon ntt
  if(ntt_ != NULL) {
    for(int k = 0; k < distrs_->length(); k++) {
      nttdistrtots_[k] = 0.0;
      for(int n = 0; n < 3; n++) {
	nttdistrtots_[k] += ((*nttpriors_)[k])[n] * ntt_->getPosDistr()->getProb(n);
      }
    }
  }

  if(ntt_ != NULL) {
    // now readjust nttpriors 
    if(nttdistrtots_[0] > 0.0 && ntt_->getPosDistr()->getProb(1) > 0.0 && ntt_->getPosDistr()->getProb(0) > 0.0 && nttdistrtots_[0] > 0.0) {
      float desired = (((*nttpriors_)[0])[2] * ntt_->getPosDistr()->getProb(2) + 0.5 * ((*nttpriors_)[0])[1] * ntt_->getPosDistr()->getProb(1)) / nttdistrtots_[0];
      ////std::cerr << "desired: " << desired << std::endl;
      // only move distributions to the right...
      if(((*nttpriors_)[1])[1] * ntt_->getPosDistr()->getProb(1) < desired) {
	((*nttpriors_)[1])[1] = desired * nttdistrtots_[1] / ntt_->getPosDistr()->getProb(1);
	if(((*nttpriors_)[1])[1] > 1.0)
	  ((*nttpriors_)[1])[1] = 1.0;
	if(((*nttpriors_)[1])[1] < 0.0)
	  ((*nttpriors_)[1])[1] = 0.0;
	((*nttpriors_)[0])[1] = 1.0 - ((*nttpriors_)[1])[1];
	((*nttpriors_)[1])[0] = (1 - desired) * nttdistrtots_[1] / ntt_->getPosDistr()->getProb(0);
	if(((*nttpriors_)[1])[0] > 1.0)
	  ((*nttpriors_)[1])[0] = 1.0;
	if(((*nttpriors_)[1])[0] < 0.0)
	  ((*nttpriors_)[1])[0] = 0.0;
	((*nttpriors_)[0])[0] = 1.0 - ((*nttpriors_)[1])[0];
      ////std::cerr << "(" << ((*nttpriors_)[1])[1] << "," << ((*nttpriors_)[0])[1] << "," << ((*nttpriors_)[1])[0] << ")" << std::endl;
      } // if must move to the right
    }
  } // if ntt not NULL
  	  
  // here can remove all distributions with prior below minimum
  float MIN_DISTR_PRIOR = 0.025;
  for(int k = 1; k < distrs_->length(); k++) 
    if(priors_[k] < MIN_DISTR_PRIOR) {
      ////std::cerr << "removing " << k << " distr" << std::endl;
      for(int n = 0; n < 3; n++)
	((*nttpriors_)[k])[n] = 0.0; // get rid of distribution!
    }
  ////std::cerr << "returning from update! " << std::endl;
  
  return output;

}

float DecayContinuousMultimixtureDistr::sliceWithNTT(float left_val, float right_val, int ntt) {
  float prob = 0.0;
  for(int k = 0; k < distrs_->length(); k++)
    if(strcmp((*distrtypes_)[k], "gaussian") == 0)
      prob += ((GaussianDistribution*)(*distrs_)[k])->slice(totalwts_, left_val, right_val) * ((*nttpriors_)[k])[ntt] * ntt_->getPosDistr()->getProb(ntt);
    else if(strcmp((*distrtypes_)[k], "gamma") == 0)
      prob += ((GammaDistribution*)(*distrs_)[k])->slice(totalwts_, left_val, right_val) * ((*nttpriors_)[k])[ntt] * ntt_->getPosDistr()->getProb(ntt);
    else
      //std::cerr << "error in slice" << std::endl;

  return prob;
}

void DecayContinuousMultimixtureDistr::printDistr() {
  for(int k = 0; k < distrs_->length(); k++) {
    printf("%s %d: prior: %0.2f ", (*defs_)[k], k+1, priors_[k]);
    (*distrs_)[k]->printDistr();
  }
}

void DecayContinuousMultimixtureDistr::writeDistr(FILE* fout) {

  fprintf(fout, "(mixture distr priors:");
  for(int ntt = 0; ntt < 3; ntt++) {
    fprintf(fout, " ntt%d [", ntt);
    for(int k = 0; k < distrs_->length(); k++) {
      if(nttpriors_ != NULL)
	fprintf(fout, "%0.2f", ((*nttpriors_)[k])[ntt]);
      else 
	fprintf(fout, "0");
      if(k < distrs_->length() - 1)
	fprintf(fout, ",");
    } // next distr
    fprintf(fout, "]");
  } // next ntt
  fprintf(fout, ")\n");
  float nttpos[] = {0.0, 0.0, 0.0};
  if(ntt_ != NULL) {
    for(int k = 0;  k < 3; k++) {
      nttpos[k] = ntt_->getPosDistr()->getProb(k);
    }
  }
  
  for(int k = 0; k < distrs_->length(); k++) {
    fprintf(fout, "            %s prior: %0.2f ", (*defs_)[k], priors_[k]);
    (*distrs_)[k]->writeDistr(fout);
    
    // here add ntt distributions for each distr[k]
    if(nttpriors_ != NULL && priors_ != NULL && priors_[k] > 0.0 && nttdistrtots_ != NULL && nttdistrtots_[k] > 0.0) 
      ////std::cerr << "K: " << k << " " << ((*nttpriors_)[k])[0] << std::endl; //" " << nttpos[0] << std::endl;
      fprintf(fout, "\t\t(ntt=0 %0.2f ntt=1 %0.2f ntt=2 %0.2f)\n", ((*nttpriors_)[k])[0] * nttpos[0]/ nttdistrtots_[k],((*nttpriors_)[k])[1] * nttpos[1]/ nttdistrtots_[k],((*nttpriors_)[k])[2] * nttpos[2]/ nttdistrtots_[k]); 
    else fprintf(fout, "\t\t(no nott info)\n");
    
  }

}
