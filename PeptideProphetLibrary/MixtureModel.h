#ifndef MIX_MODEL_H
#define MIX_MODEL_H

#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "SequestResult.h"
#include "NMCMixtureDistr.h"
#include "NTTMixtureDistr.h"
#include "GlycMixtureDistr.h"
#include "DiscrimValMixtureDistr.h"
#include "MascotDiscrimValMixtureDistr.h"
#include "ICATMixtureDistr.h"
#include "MixtureDistr.h"
#include "Array.h"
#include "Spectrum.h"
#include "OrderedResult.h"
#include "MassDifferenceDiscrMixtureDistr.h"
#include "DecayContinuousMultimixtureDistr.h"
#include "VariableOffsetMassDiffDiscrMixtureDistr.h"
#include "OutputContent.h"
#include "DatasetNumMap.h"


/*

Program       : MixtureModel for PeptideProphet                                                       
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 

Primary data object holding all mixture distributions for each precursor ion charge

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


int comp_nums(const void* num1, const void* num2);
int comp_specs(const void* num1, const void* num2);
int comp_ords(const void* num1, const void* num2);

class MixtureModel {

 public:

  MixtureModel(std::vector<SequestResult> &results, int maxnumiters, Boolean2 icat, Boolean2 glyc, Boolean2 massd, Boolean2 mascot, Boolean2 qtof, char *enzyme);
  Boolean2 updateDistr(char* name, int charge);
  Boolean2 iterate(int counter, int charge);
  Boolean2 iterate(int counter);
  float computeProb(int charge, int index);
  Boolean2 updatePriors(int charge);
  Boolean2 updateProbs(int counter, int charge); //Xiuxia
  float getProb(int charge, int index);
 
  virtual Boolean2 enterDistribution(int charge, int index, const char* tag, char* value);  
//  void writeResults(std::vector<SequestResult> &sequestResults, char* filename);
  void writeResults(char* filename);    //Xiuxia, removed sequestResults from the argument
  float getAdjDoublyTriplyProb(float prob_2_adj, float prob_of_partner);
  void printDistr();
  void readData(std::vector<SequestResult> &results);
  void setNegativeDiscrimValDistrs();
  void writeDiscrimValDistr(char* filename);
  void computeEstimatedSensAndError(char* filename);
  void computeAdjDoubleTriplySpectraProbs();
  Boolean2 isAdjusted(int index);
  void deriveModel(int maxnumiters);
  Boolean2 negOnly(int charge);
  void writeDistr(char* filename);
  void writeModelfile(char* modelfile);
  int getNegOnlyCharge(int charge);
  float getTotalProb(int charge);
  Boolean2 isNumber(char c);
  Boolean2 getMaldiInfo(char* filename);
  float computeProbWithNTT(int charge, int origntt, float orgiprob, int ntt, int index);
  void setMixtureDistributionNames(char* discrim, char* ntt);
  void validateModel();
//  void writeResultsInOrder(std::vector<SequestResult> &sequestResults, char* filename);
  void writeResultsInOrder(char* filename);    //Xiuxia, removed sequestResults from the argument
  char* getTagValue(char* data, char*tag);
  char* getEnzyme(char* filename);
  void writeResultsOrdered(const char* filename, std::vector<DatasetNumMap> &vecDatasetMap) ;
  static bool SortOutputResultsByHitnum(OutputContent &a, OutputContent &b) ; //Xiuxia, 08/02/2006


 protected:

    bool abort;

 MixtureDistr* getMixtureDistr(char* name, int charge);

 Array<Array<char*>*>* spectra_;  // spectra by precursor ion charge (0:1+, 1:2+, 2:3+)
 Array<Array<float>*>* probs_;    // probs by charge
 Array<Array<MixtureDistr*>*>* mixtureDistrs_;  // mixture distributions by charge
 Array<Array<int>*>* inds_;

 Array<Array<int>*>* dataset_num_All_ ; //Xiuxia
 Array<Array<int>*>* scanNumberAll_ ;    //Xiuxia
 Array<Array<float>*>* xcorrAll_ ;    //Xiuxia
 Array<Array<float>*>* fvalAll_ ; //Xiuxia
 Array<Array<float>*>* deltaCn2All_ ; //Xiuxia

 float* priors_; // for each charge
 int* numspectra_;
 int* done_;  // iterations to model termination
 Boolean2 icat_;  // use icat peptide cys information for probs
 Boolean2 glyc_;  // use peptide N-glycosylation motif for probs
 Boolean2 use_adj_probs_;  // constraint for 2+ and 3+ interpretations of same spectrum
 Boolean2 gamma_;  // use gamma distribution for fval distribution among incorrect search results
 int num_pairs_;
 Spectrum* pairedspecs_;  
 Array<int>* adjusted_;
 int max_num_iters_;
 Boolean2* negOnly_;  // no model results, use alternative means to crudely estimate prob ('0' or '-charge')
 Boolean2 maldi_;  // maldi data
 Boolean2 maldi_set_;
 Boolean2* pseudonegs_;
 int min_num_specs_;
 char* discrim_name_;
 char* ntt_name_;
 char* enzyme_;
 Boolean2 qtof_;
};

#endif
