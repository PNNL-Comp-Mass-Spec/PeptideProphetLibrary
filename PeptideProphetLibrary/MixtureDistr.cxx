#include "MixtureDistr.h"
#include <sstream>

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


MixtureDistr::MixtureDistr() 
{ 
    posdistr_ = NULL ; 
    negdistr_ = NULL ; 
}

MixtureDistr::MixtureDistr(int charge, char* name, char* tag) 
{
    posdistr_ = NULL ; 
    negdistr_ = NULL ; 
    initializeDistr(charge, name, tag);
}

MixtureDistr::~MixtureDistr() { }

void MixtureDistr::initializeDistr(int charge, char* name, char* tag) {
  charge_ = charge;
  intvals_ = NULL;
  floatvals_ = NULL;
  if(name == NULL) {
    //std::cerr << "null name for MixtureDistr" << std::endl;
      throw gcnew System::Exception("null name for MixtureDistr");
    //exit(1);
  }
  if(tag == NULL) {
    //std::cerr << "null tag for MixtureDistr" << std::endl;
      throw gcnew System::Exception("null tag for MixtureDistr");
    //exit(1);
  }

  name_ = new char[strlen(name)+1];
  strcpy(name_, name);
  name_[strlen(name)] = 0;
  tag_ = new char[strlen(tag)+1];
  strcpy(tag_, tag);
  tag_[strlen(tag)] = 0;
  pospriors_ = NULL;
  negpriors_ = NULL;
  negOnly_ = False;
}

int MixtureDistr::getNumVals() { 
  if(intvals_ != NULL) {
    return intvals_->length();
  }
  return floatvals_->length();
}

void MixtureDistr::enter(int index, char* val) {
  
  if(intvals_ != NULL) {
    enter(index, inttranslate(val));
  }
  else {
    enter(index, floattranslate(val));
  }
}


Boolean2 MixtureDistr::isValue(int index, int val) {
  if(intvals_ == NULL) {
    return False;
  }
  return (*intvals_)[index] == val;
}




int MixtureDistr::inttranslate(char* val) 
{
  return atoi(val);
}

float MixtureDistr::floattranslate(char* val) 
{
  return atof(val);
}

float MixtureDistr::getPosProb(int index) {
  if(intvals_ != NULL) {
    if(index < 0 || index >= intvals_->length()) {
      //std::cerr << "violation of index " << index << " for " << intvals_->length() << std::endl;
      //exit(1);
        std::stringstream str;
        str << "violation of index " << index << " for " << intvals_->length() << std::endl;
        throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
    }
    return posdistr_->getProb((*intvals_)[index]);
  }
  else {
    return posdistr_->getProb((*floatvals_)[index]);
  }
}


float MixtureDistr::getNegProb(int index) {
  if(intvals_ != NULL) {
    return negdistr_->getProb((*intvals_)[index]);
  }
  else {
    return negdistr_->getProb((*floatvals_)[index]);
  }
}

Boolean2 MixtureDistr::negOnly() {
  return negOnly_;
}

Distribution* MixtureDistr::getPosDistr() {
  return posdistr_;
}

void MixtureDistr::setPosDistr(MixtureDistr* distr) {
  Distribution* reject = distr->getPosDistr();
  if(posdistr_ != NULL)//Xiuxia, should it be reject here, instead of posdistr_?
    delete posdistr_;

  posdistr_ = distr->getPosDistr();//Xiuxia, why not just equal to reject?

}

    
void MixtureDistr::enter(int index, float val) {
  floatvals_->insertAtEnd(val);
}


void MixtureDistr::enter(int index, int val) {//Xiuxia, index is not used here
  intvals_->insertAtEnd(val);

  //std::cerr << intvals_[index] << std::endl;
}
    
Boolean2 MixtureDistr::update(Array<float>* probs) {
    //Xiuxia, get new mean and new stdev for both the Gaussian and the Gamma distribution
  
    //Xiuxia, initialize newtot_[0] = 0, newtotsq_ = 0, newtotwt_ = 0
    posdistr_->initUpdate(NULL);      
    negdistr_->initUpdate(NULL);

    for(int k = 0; k < probs->length(); k++) {

        if(intvals_ != NULL) {
            posdistr_->addVal((*probs)[k], (*intvals_)[k]);
            negdistr_->addVal(1.0 - (*probs)[k], (*intvals_)[k]);
        }
        else {
            posdistr_->addVal((*probs)[k], (*floatvals_)[k]);
            negdistr_->addVal(1.0 - (*probs)[k],(*floatvals_)[k]);
        }

    }
  
  Boolean2 output = False;
  
  if(posdistr_->update()) {
    output = True;
  }
  if(negdistr_->update()) {
    output = True;
  }
  return output;    //Xiuxia, negdistr_ takes precedence over posdistr_
}

char* MixtureDistr::getName() { return name_; }

char* MixtureDistr::getTag() { return tag_; }

void MixtureDistr::printVal(int index) {
  
  //std::cout << name_ << "=";
  if(intvals_ != NULL) {
    //std::cout << (*intvals_)[index];
  }
  else {
    //std::cout << (*floatvals_)[index];
  }
  
}


void MixtureDistr::printDistr() {
  //std::cout << name_ << std::endl;
  //std::cout << "\tpos: ";
  posdistr_->printDistr();
  //std::cout << "\tneg: ";
  negdistr_->printDistr();
}

void MixtureDistr::writeDistr(FILE* fout) {
  fprintf(fout, "%s\n", name_);
  fprintf(fout, "\tpos: ");
  posdistr_->writeDistr(fout);
  fprintf(fout, "\tneg: ");
  negdistr_->writeDistr(fout);
}
