#include "ContinuousMultimixtureDistr.h"
#include<sstream>

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

ContinuousMultimixtureDistr::ContinuousMultimixtureDistr(float maxdiff) : ContinuousDistribution(maxdiff) {

  distrs_ = NULL;
  priors_ = NULL;
  zvals_ = NULL;
  defs_ = NULL;

}

int ContinuousMultimixtureDistr::getNumDistributions() { 
  if(distrs_ == NULL)
    return 0;
  return distrs_->length();
}

void ContinuousMultimixtureDistr::initialize() 
{
  distrs_ = new Array<ContinuousDistribution*>;
  priors_ = NULL;
  zvals_ = new Array<float*>;
  vals_ = new Array<float>;
  defs_ = new Array<char*>;
  wts_ = new Array<float>;
  distrtypes_ = new Array<char*>;
  mean_ = 0.0;
  stdev_ = 0.0; // will never change
}

ContinuousMultimixtureDistr::~ContinuousMultimixtureDistr() {
  if(distrs_ != NULL) {
    for(int k = 0; k < distrs_->length(); k++)
      if((*distrs_)[k] != NULL)
	delete (*distrs_)[k];
    delete distrs_;
  }
  if(priors_ != NULL)
    delete [] priors_;
  if(zvals_ != NULL) {
    for(int k = 0; k < zvals_->length(); k++)
      if((*zvals_)[k] != NULL)
	delete[] (*zvals_)[k];
    delete zvals_;
  }
  if(vals_ != NULL)
    delete vals_;
  if(defs_ != NULL) {
    for(int k = 0; k < defs_->length(); k++)
      if((*defs_)[k] != NULL)
	delete (*defs_)[k];
    delete defs_;
  }
    
  
}

float ContinuousMultimixtureDistr::getProb(float val) {

  float prob = 0.0;
  for(int k = 0; k < distrs_->length(); k++)
    prob += getMixtureDistrProb(k, val) * priors_[k];
  ////std::cerr << "returning " << prob << " for " << val << std::endl;
  return prob;
}

Boolean2 ContinuousMultimixtureDistr::violation(int leftdistr, int rightdistr) {
  return (*distrs_)[leftdistr]->getMean() >= (*distrs_)[rightdistr]->getMean();
}

void ContinuousMultimixtureDistr::removeViolatingDistrs() {
  for(int k = 1; k < distrs_->length(); k++) {
    if((*distrs_)[k] != (*distrs_)[k-1] && violation(k, k-1)) {
      if((*distrs_)[k] != NULL) 
	delete (*distrs_)[k];
      distrs_->replace(k, (*distrs_)[k-1]);
    } // if violation
  } // next pair of adjacent distributions
}

float ContinuousMultimixtureDistr::getProb(int val) { return 0.0; }

void ContinuousMultimixtureDistr::addVal(float wt, int val) { }

float ContinuousMultimixtureDistr::getMean() {
  float tot = 0.0;
  for(int k = 0; k < distrs_->length(); k++)
    tot += priors_[k] * (*distrs_)[k]->getMean();
  return tot;
}
float ContinuousMultimixtureDistr::getStdev() {
  float tot = 0.0;
  for(int k = 0; k < distrs_->length(); k++)
    tot += priors_[k] * (*distrs_)[k]->getStdev();
  return tot;
}

float ContinuousMultimixtureDistr::getMixtureDistrProb(int k, float val) {
  ////std::cerr << (*defs_)[k] << std::endl;
  if(strcmp((*distrtypes_)[k], "gaussian") == 0)
      return ((GaussianDistribution*)(*distrs_)[k])->getProb(val);
  else if(strcmp((*distrtypes_)[k], "gamma") == 0)
      return ((GammaDistribution*)(*distrs_)[k])->getProb(val);
  // error

  //std::cerr << "error with getMixtureProb" << std::endl;
  return 0.0;
}

void ContinuousMultimixtureDistr::addVal(float wt, float val) {
  if(vals_->length() <= index_) {
    vals_->insertAtEnd(val);
    // set initial zvalues based upon distr's
    float totprob = 0.0;
    float* nextz = new float[distrs_->length()];
    int preset_mixture = presetZvalue(index_);
    for(int k = 0; k < distrs_->length(); k++) {
      if(preset_mixture < 0)
	nextz[k] = getMixtureDistrProb(k, val) * priors_[k];
      else if(k == preset_mixture)
	nextz[k] = 1.0;
      else
	nextz[k] = 0.0;
      totprob += nextz[k];
    }
    if(preset_mixture < 0 && totprob > 0) 
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

  for(int k = 0; k < distrs_->length(); k++)
    (*distrs_)[k]->addVal(wt * ((*zvals_)[index_])[k], val);

  assert(vals_->length() == zvals_->length());
  index_++;
}

int ContinuousMultimixtureDistr::presetZvalue(int index) { return -1; }

void ContinuousMultimixtureDistr::addComponent(float* settings, char* distr, char* def){

  if(strcmp(distr, "gaussian") == 0) {
      distrs_->insertAtEnd(new GaussianDistribution(maxdiff_));
  }
  else if(strcmp(distr, "gamma") == 0)
    distrs_->insertAtEnd(new GammaDistribution(maxdiff_));
  else {
    //std::cerr << "cannot accommodate " << distr << " distribution" << std::endl;
    //exit(1);
	  std::stringstream str;
	  str << "cannot accommodate " << distr << " distribution" << std::endl;
	  throw new System::Exception(str.str().c_str());
  }
  char* newdistr = new char[strlen(distr)+1];
  strcpy(newdistr, distr);
  newdistr[strlen(distr)] = 0;
  distrtypes_->insertAtEnd(newdistr);
  (*distrs_)[distrs_->length()-1]->init(settings);

  char* next = new char[strlen(def)+1];
  strcpy(next, def);
  next[strlen(def)] = 0;
  defs_->insertAtEnd(next);
}

// must initialize each specified distribution via initComponent
void ContinuousMultimixtureDistr::init(float* prior) {
  ContinuousDistribution::init(NULL);
  //initUpdate(prior);
}

void ContinuousMultimixtureDistr::setMinVal(float min) { }

void ContinuousMultimixtureDistr::initUpdate(float* prior) {
  index_ = 0; // reset
  for(int k = 0; k < distrs_->length(); k++) {
    (*distrs_)[k]->initUpdate(prior);
  }
}

Boolean2 ContinuousMultimixtureDistr::update() {
  if(vals_ == NULL || zvals_ == NULL || vals_->length() != zvals_->length()) {
	  std::stringstream str; 
	  str << "error in ContinuousMultimixtureDistr update...." << std::endl;
    //std::cerr << "error in ContinuousMultimixtureDistr update...." << std::endl;
    if(vals_ == NULL)
		str << "null vals_" << std::endl;
      //std::cerr << "null vals_" << std::endl;
    else if(zvals_ == NULL)
		str << "null zvals_" << std::endl;
      //std::cerr << "null zvals_" << std::endl;
    else
		str << vals_->length() << " vals_ <> " << zvals_->length() << " zvals_" << std::endl;
      //std::cerr << vals_->length() << " vals_ <> " << zvals_->length() << " zvals_" << std::endl;

	throw new System::Exception(str.str().c_str());

//    exit(1);
  }
  //set_ = True;
  ////std::cerr << " updating...." << std::endl;
  //std::cerr << " start priors: " << priors_[0] << ":" << priors_[1] << std::endl;

  assert(vals_->length() == index_);
  float totalprob = 0.0;
  Boolean2 output = False;
  for(int k = 0; k < distrs_->length(); k++)
    if((*distrs_)[k]->update())
      output = True;

  if(! output)
    return output;

  mean_ = getMean(); // update
  stdev_ = getStdev();

  // now update zvals and priors
  totalwts_ = 0.0;
  float nextprob;
  for(int i = 0; i < vals_->length(); i++) {
    totalwts_ += (*wts_)[i];
    int preset_mixture = presetZvalue(i);
    nextprob = 0.0;
    if(preset_mixture < 0) {
      for(int k = 0; k < distrs_->length(); k++) {
	((*zvals_)[i])[k] = getMixtureDistrProb(k, (*vals_)[i]);
	nextprob += ((*zvals_)[i])[k];
      }

    }
    // now normalize
    if(preset_mixture >= 0 || nextprob > 0) 
      for(int k = 0; k < distrs_->length(); k++) {
	if(preset_mixture < 0)
	  ((*zvals_)[i])[k] /= nextprob;
	else if(k == preset_mixture)
	  ((*zvals_)[i])[k] = 1.0;
	else 
	  ((*zvals_)[i])[k] = 0.0;
	//priors_[k] += ((*zvals_)[i])[k];
	//(*priors_)[k] += ((*zvals_)[i])[k];
	//totalprob += ((*zvals_)[i])[k];
      }

    ////std::cerr << "index " << i << ": [" << (*distrs_)[0]->getProb((*vals_)[i]) << " <=> " << (*distrs_)[1]->getProb((*vals_)[i]) << "]  " << (*vals_)[i] << " (" << ((*zvals_)[i])[0] << " " << ((*zvals_)[i])[1] << ")" << " wt: " << (*wts_)[i] << std::endl;

  } // next data
  // now recompute priors
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

  ////std::cerr << " end priors: " << priors_[0] << ":" << priors_[1] << std::endl;
  ////std::cerr << "total wt: " << totalwts_ << std::endl;
  return output;

}

void ContinuousMultimixtureDistr::printDistr() {
  for(int k = 0; k < distrs_->length(); k++) {
    printf("%s %d: prior: %0.2f ", (*defs_)[k], k+1, priors_[k]);
    (*distrs_)[k]->printDistr();
  }
}

void ContinuousMultimixtureDistr::writeDistr(FILE* fout) {
  fprintf(fout, "\n");
  for(int k = 0; k < distrs_->length(); k++) {
    fprintf(fout, "            %s prior: %0.2f ", (*defs_)[k], priors_[k]);
    (*distrs_)[k]->writeDistr(fout);
  }
}

float ContinuousMultimixtureDistr::slice(float num, float left_val, float right_val) {
  float prob = 0.0;
  for(int k = 0; k < distrs_->length(); k++)
    if(strcmp((*distrtypes_)[k], "gaussian") == 0)
      prob += ((GaussianDistribution*)(*distrs_)[k])->slice(num, left_val, right_val) * priors_[k];
    else if(strcmp((*distrtypes_)[k], "gamma") == 0)
      prob += ((GammaDistribution*)(*distrs_)[k])->slice(num, left_val, right_val) * priors_[k];
    //else
      //std::cerr << "error in slice" << std::endl;

  return prob;
}

float ContinuousMultimixtureDistr::slice(float left_val, float right_val) {
  return slice(totalwts_, left_val, right_val);
}
