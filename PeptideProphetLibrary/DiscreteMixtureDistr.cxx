#include "DiscreteMixtureDistr.h"

/*

Program       : DiscreteMixtureDistr for PeptideProphet                                                       
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


DiscreteMixtureDistr::DiscreteMixtureDistr(int charge, int numbins, char* name, char* tag) : MixtureDistr(charge, name, tag) {
  numbins_ = numbins;
  ////std::cout << "...with " << numbins_ << " bins" << std::endl;
  maxdiff_ = 0.002;
  priors_ = NULL;
  numpos_priors_ = 2.0;
  numneg_priors_ = 2.0;
  bindefs_ = NULL ; 
}

Boolean2 DiscreteMixtureDistr::update(Array<float>* probs) {
  if(priors_ == NULL && intvals_->length() > 0) 
  {
    priors_ = new float[numbins_];
    for(int k = 0; k < numbins_; k++) 
	{
      priors_[k] = 0.0;
    }
	int numVals = intvals_->length() ; 
    for(int k = 0; k < numVals ; k++) 
	{
		int currentBinNum = (*intvals_)[k] ; 
		priors_[currentBinNum] += 1.0;
    }
    for(int k = 0; k < numbins_; k++) 
	{
      priors_[k] /= numVals ;
    }
    ((DiscreteDistribution*)(posdistr_))->setPriors(priors_, numpos_priors_);
    ((DiscreteDistribution*)(negdistr_))->setPriors(priors_, numneg_priors_);
  } // if no priors yet
  return MixtureDistr::update(probs);
}


void DiscreteMixtureDistr::initializeBinDefs(char** bindefs) 
{
	if (bindefs_ != NULL)
		delete bindefs_ ; 
	bindefs_ = new Array<char*>;
	for(int k = 0; k < numbins_; k++) 
	{
		bindefs_->insertAtEnd(bindefs[k]);
	}
}

void DiscreteMixtureDistr::init(char** bindefs) {
  if(bindefs != NULL)
    initializeBinDefs(bindefs);
  if(intvals_ != NULL)
  {
	  delete intvals_ ;
  }
  intvals_ = new Array<int>;

  if (posdistr_ != NULL)
  {
	  delete posdistr_ ; 
  }
  if (negdistr_ != NULL)
  {
	  delete negdistr_ ; 
  }

  posdistr_ = new DiscreteDistribution(numbins_, maxdiff_);
  negdistr_ = new DiscreteDistribution(numbins_, maxdiff_);
}


void DiscreteMixtureDistr::printDistr() {
  //std::cout << name_ << std::endl;
  //std::cout << "\tpos: ";
  printPosDistribution();
  //std::cout << "\tneg: ";
  printNegDistribution();
}

void DiscreteMixtureDistr::writeDistr(FILE* fout) {
  fprintf(fout, "%s\n", name_);
  fprintf(fout, "\tpos: ");
  fprintf(fout, "(");
  for(int k = 0; k < numbins_; k++) {
    fprintf(fout, "%0.3f %s", posdistr_->getProb(k), (*bindefs_)[k]);
    if(k < numbins_ - 1) {
      fprintf(fout, ", ");
    }
  }
  fprintf(fout, ")\n");
  fprintf(fout, "\tneg: ");
  fprintf(fout, "(");
  for(int k = 0; k < numbins_; k++) {
    fprintf(fout, "%0.3f %s", negdistr_->getProb(k), (*bindefs_)[k]);
    if(k < numbins_ - 1) {
      fprintf(fout, ", ");
    }
  }
  fprintf(fout, ")\n");
}



void DiscreteMixtureDistr::printPosDistribution() {
  printf("(");
  for(int k = 0; k < numbins_; k++) {
    printf("%0.3f %s", posdistr_->getProb(k), (*bindefs_)[k]);
    if(k < numbins_ - 1) {
      printf(", ");
    }
  }
  printf(")\n");
}

void DiscreteMixtureDistr::printNegDistribution() {
  printf("(");
  for(int k = 0; k < numbins_; k++) {
    printf("%0.3f %s", negdistr_->getProb(k), (*bindefs_)[k]);
    if(k < numbins_ - 1) {
      printf(", ");
    }
  }
  printf(")\n");
}


