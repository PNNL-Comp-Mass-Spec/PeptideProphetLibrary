#include "MixtureModel.h"
#include "GlobalVariable.h"
#include <iostream>
#include <sstream>
#include <algorithm>


using namespace std ;

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
/*
	DJ: Deep Jaitly. Tried to keep annotating what I was changing, but at some points, I ran out of patience
	because of the amount I was changing.
*/

MixtureModel::MixtureModel(std::vector<SequestResult> &results, int max_num_iters, Boolean2 icat, Boolean2 glyc, Boolean2 massd, Boolean2 mascot, Boolean2 qtof, char *enzyme)
{

  max_num_iters_ = max_num_iters;
  icat_ = icat;
  glyc_ = glyc;
  qtof_ = qtof;
  /*priors_ = new float[3];*/
  //numspectra_ = new int[3];
  //done_ = new int[3];
  priors_ = new float[numCharge];//Xiuxia
  numspectra_ = new int[numCharge];//Xiuxia
  done_ = new int[numCharge];//Xiuxia
  spectra_ = new Array<Array<char*>*>;
  probs_ = new Array<Array<float>*>;
  inds_ = new Array<Array<int>*>;
  mixtureDistrs_ = new Array<Array<MixtureDistr*>*>;

  dataset_num_All_ = new Array<Array<int>*> ; //Xiuxia
  scanNumberAll_ = new Array<Array<int>*> ;	//Xiuxia
  xcorrAll_ = new Array<Array<float>*> ;	//Xiuxia
  fvalAll_ = new Array<Array<float>*> ;		//Xiuxia
  deltaCn2All_ = new Array<Array<float>*> ;		//Xiuxia

  gamma_ = True;
  //use_adj_probs_ = True;
  use_adj_probs_ = False;
  pairedspecs_ = NULL;
  adjusted_ = NULL;
  negOnly_ = new Boolean2[numCharge];//Xiuxia
  maldi_ = False;
  maldi_set_ = True;
  pseudonegs_ = new Boolean2[numCharge];//Xiuxia
  min_num_specs_ = 35; //50;
  discrim_name_ = NULL;
  ntt_name_ = NULL;
  enzyme_ = NULL;

  //std::cerr << "\n PeptideProphet v. 1.0 A.Keller 11.7.02 ISB \n" << std::endl;
  //std::cerr << " Modified by Xiuxia Du and Deep Jaitly, June 21, 2006" << std::endl << std::endl ;

  enzyme_ = NULL ;	//Xiuxia, redundant assignment
  if (enzyme != NULL)
  {
	  enzyme_ = new char [strlen(enzyme)+1] ;
	  strcpy(enzyme_, enzyme) ;
	  enzyme_[strlen(enzyme)] = '\0' ;	//Xiuxia, terminating character
  }

  char enz_dig_distr_name[100];
  strcpy(enz_dig_distr_name, "no. ");
  if(enzyme_ == NULL)
    strcat(enz_dig_distr_name, "tryptic");
  else strcat(enz_dig_distr_name, enzyme_);
  strcat(enz_dig_distr_name, " term.");


  setMixtureDistributionNames("discrim score", enz_dig_distr_name); // fval and ntt mixture model names

  for(int charge = 0; charge < numCharge; charge++) {//Xiuxia

    numspectra_[charge] = 0;
    done_[charge] = -1;
    priors_[charge] = -1.0;
    negOnly_[charge] = False;
    pseudonegs_[charge] = False;
    Array<char*>* nextspec = new Array<char*>;
    spectra_->insertAtEnd(nextspec);
    Array<float>* nextprob = new Array<float>();
    probs_->insertAtEnd(nextprob);
    Array<int>* nextind = new Array<int>;
    inds_->insertAtEnd(nextind);

	//Xiuxia
	Array<int>* nextDatasetNum = new Array<int>() ;
	dataset_num_All_ -> insertAtEnd(nextDatasetNum) ;
	Array<int>* nextScanNumber = new Array<int>();
	scanNumberAll_ ->insertAtEnd(nextScanNumber) ;
	Array<float>* nextXCorr = new Array<float>();
	xcorrAll_ ->insertAtEnd(nextXCorr) ;
	Array<float>* nextFval = new Array<float>() ;
	fvalAll_->insertAtEnd(nextFval) ;
	Array<float>* nextDeltaCn2 = new Array<float>() ;
	deltaCn2All_ ->insertAtEnd(nextDeltaCn2) ;

    Array<MixtureDistr*>* next = new Array<MixtureDistr*>;	//Xiuxia, for each charge
	//                    add new distributions to model here

    /////////////////////////////////////////////////////////////////////////////////////////////
    if(mascot)
      next->insertAtEnd(new MascotDiscrimValMixtureDistr(charge, discrim_name_, "fval", maldi_, qtof_));
    else
      next->insertAtEnd(new DiscrimValMixtureDistr(charge, discrim_name_, "fval", gamma_, maldi_, qtof_));

    next->insertAtEnd(new NTTMixtureDistr(charge, ntt_name_, "ntt"));
    next->insertAtEnd(new NMCMixtureDistr(charge, "no. missed cl", "nmc"));

    if(massd) // now only available for charge 2/3
      //next->insertAtEnd(new MassDifferenceDiscrMixtureDistr(charge, "mass diff", "massd", 5.0, 1.0));
      next->insertAtEnd(new VariableOffsetMassDiffDiscrMixtureDistr(charge, "var offset mass diff", "massd", 5.0, 1.0, 0.0));
    if(icat_) {
      next->insertAtEnd(new ICATMixtureDistr(charge, "icat cys", "pep"));
    }
    if(glyc_) {
      next->insertAtEnd(new GlycMixtureDistr(charge, "N glyc motif", "pep"));
    }
    /////////////////////////////////////////////////////////////////////////////////////////////
    mixtureDistrs_->insertAtEnd(next);
  }


  readData(results);

  validateModel();

  deriveModel(max_num_iters_);

}

void MixtureModel::setMixtureDistributionNames(char* discrim, char* ntt) {
  discrim_name_ = new char[strlen(discrim)+1];
  strcpy(discrim_name_, discrim);
  discrim_name_[strlen(discrim)] = 0;
  ntt_name_ = new char[strlen(ntt)+1];
  strcpy(ntt_name_, ntt);
  ntt_name_[strlen(ntt)] = 0;
}

// after data entry, make sure have same number of data points in each mixture distr
void MixtureModel::validateModel() {
  int tot = 0;

  //for(int charge = 0; charge < 3; charge++) {
  for(int charge = 0; charge < numCharge; charge++) {//Xiuxia
    for(int k = 0; k < (*mixtureDistrs_)[charge]->length(); k++) {
      tot += (*spectra_)[charge]->length();
      if((*spectra_)[charge]->length() != (*(*mixtureDistrs_)[charge])[k]->getNumVals()) {
			//std::cerr << "have " << (*spectra_)[charge]->length() << " spectra and " << (*(*mixtureDistrs_)[charge])[k]->getNumVals() << " values for " << (*(*mixtureDistrs_)[charge])[k]->getName() << " charge " << (charge+1) << std::endl;
		  std::stringstream str;
		  str <<  "have " << (*spectra_)[charge]->length() << " spectra and " << (*(*mixtureDistrs_)[charge])[k]->getNumVals() << " values for " << (*(*mixtureDistrs_)[charge])[k]->getName() << " charge " << (charge+1) << std::endl;
		  throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
			//	exit(1);
      }
    } // next distr
  } // next charge
  if(tot == 0) { // nothing to do...
    //std::cerr << " read in no data" << std::endl;
    //exit(1);
	  std::stringstream str("");
	  str << " read in no data" << std::endl;
	  throw gcnew System::Exception(gcnew System::String(str.str().c_str()));

  }
}


// initialize fval negative distributions with 0 tryptic termini data, if available
void MixtureModel::setNegativeDiscrimValDistrs()
{
  for(int charge = 0; charge < numCharge; charge++) //Xiuxia, there is an outer loop wrt charge, why loop again wrt charge here?
  {
	  MixtureDistr* discrimDistr = getMixtureDistr(discrim_name_, charge) ;
	  MixtureDistr* nttDistr = getMixtureDistr(ntt_name_, charge) ;

	  if(discrimDistr != NULL && nttDistr != NULL)
	  {
  		  DiscrimValMixtureDistr *discrimMixDistr = (DiscrimValMixtureDistr*) discrimDistr ; //Xiuxia, upgrade to daughter class directly from the mother class
		  NTTMixtureDistr *nttMixDistr = (NTTMixtureDistr*) nttDistr ;
		  pseudonegs_[charge] = discrimMixDistr->initializeNegDistribution(nttMixDistr); //Xiuxia, initialized to False
	  }
  }

}


MixtureDistr* MixtureModel::getMixtureDistr(char* name, int charge)
{
  int numDistr = (*mixtureDistrs_)[charge]->length() ;

  for(int k = 0; k < numDistr ; k++)
  {
    MixtureDistr* distr = (*(*mixtureDistrs_)[charge])[k] ;
	char *distrName = distr->getName() ;
    if(strcmp(distrName, name) == 0) {
      return distr;	//Xiuxia, get the pointer to the mixture distribution. This is initialized earlier.
    }
  }
  return NULL;
}


float MixtureModel::getTotalProb(int charge) {
  if(negOnly_[charge]) {
    return 0.0;
  }
  float tot = 0.0;
  for(int k = 0; k < (*probs_)[charge]->length(); k++) {
    tot += (*(*probs_)[charge])[k];
  }
  return tot;
}

Boolean2 MixtureModel::updateDistr(char* name, int charge) {
  for(int k = 0; k < (*mixtureDistrs_)[charge]->length(); k++) {
    if(strcmp((*(*mixtureDistrs_)[charge])[k]->getName(), name) == 0) {
      return ((*(*mixtureDistrs_)[charge])[k]->update((*probs_)[charge]));
    }
  }

  return False;
}

void MixtureModel::deriveModel(int maxnumiters) {

  //setNegativeDiscrimValDistrs();

  // here if use mixture positives....then initialize with NTTMixtureDistr....
  if(maldi_ && getMixtureDistr(discrim_name_, 0) != NULL && getMixtureDistr(ntt_name_, 0) != NULL &&
     getMixtureDistr(discrim_name_, 0)->getPosDistr() != NULL) {
    ((DecayContinuousMultimixtureDistr*)(getMixtureDistr(discrim_name_, 0)->getPosDistr()))->initWithNTTs((NTTMixtureDistr*)(getMixtureDistr(ntt_name_, 0)));
    //std::cerr << "number of distrs: " << ((DecayContinuousMultimixtureDistr*)(getMixtureDistr(discrim_name_, 0)->getPosDistr()))->getNumDistributions() << std::endl;


  } // maldi
  if(qtof_)
    for(int ch = 0; ch < numCharge; ch++)
	{
		MixtureDistr* discrimDistr = getMixtureDistr(discrim_name_, ch) ;
		MixtureDistr* nttDistr = getMixtureDistr(ntt_name_, ch) ;
		if(discrimDistr != NULL &&  nttDistr != NULL)
		{
			Distribution *discrimPosDistr = discrimDistr->getPosDistr() ;
			if (discrimPosDistr != NULL)
			{
				DecayContinuousMultimixtureDistr* decayDistr = (DecayContinuousMultimixtureDistr*) discrimPosDistr ;
				decayDistr->initWithNTTs((NTTMixtureDistr*)nttDistr);
			}
		}
	}


  setNegativeDiscrimValDistrs();

  int counter = 1;

  while(counter <= maxnumiters && iterate(counter)) {
 //   if(counter == 1) {
 //     //std::cout << " iteration ";
 //   }
 //   else {
 //     std::cout << '\b';
 //     std::cout << '\b';
 //   }
 //   if(counter > 10) {
 //     std::cout << '\b';
 //   }
 //   if(counter > 100) {
 //     std::cout << '\b';
 //   }
	//std::cout << counter << " ";
 //   std::cout.flush();

	std::cout << "Iteration " << counter << "\n";
    counter++;
  }
  /*
  if(maldi_ || qtof_) {

      if(maldi_ && getMixtureDistr(discrim_name_, 0) != NULL && getMixtureDistr(ntt_name_, 0) != NULL &&
     getMixtureDistr(discrim_name_, 0)->getPosDistr() != NULL) {
	((DecayContinuousMultimixtureDistr*)(getMixtureDistr(discrim_name_, 0)->getPosDistr()))->commence();

      } // maldi
      if(qtof_)
	for(int ch = 0; ch < 3; ch++)
	  if(getMixtureDistr(discrim_name_, ch) != NULL && getMixtureDistr(ntt_name_, ch) != NULL &&
	     getMixtureDistr(discrim_name_, ch)->getPosDistr() != NULL)
	    ((DecayContinuousMultimixtureDistr*)(getMixtureDistr(discrim_name_, ch)->getPosDistr()))->commence();

      for(int ch = 0; ch < 3; ch++)
	done_[ch] = -1; // set back



      while(counter < maxnumiters && iterate(counter)) {
	if(counter == 1) {
	  //std::cout << " iteration ";
	}
	else {
	  //std::cout << '\b';
	  //std::cout << '\b';
	}
	if(counter > 10) {
	  //std::cout << '\b';
	}
	if(counter > 100) {
	  //std::cout << '\b';
	}
	//std::cout << counter << " ";
	//std::cout.flush();
	counter++;
      }
  } // if maldi or qtof
  */

  //for(int k = 0; k < 12; k++) {
  //  std::cout << '\b';
  //}
  //if(counter > 10) {
  //  std::cout << '\b';
  //}
  //if(counter > 100) {
  //  std::cout << '\b';
  //}

  std::cout << std::endl << std::endl << " Model complete after " << counter-1 << " iterations" << std::endl << std::endl;

	//std::cout << std::endl << " " << counter << " iterations for other charge states not listed above." << std::endl << std::endl;

  // if not both satisfactorily modeled
  for(int charge = 0; charge < numCharge; charge++) {

    if((*spectra_)[charge]->length() > 0 &&
       ((*spectra_)[charge]->length() < min_num_specs_ ||
       (priors_[charge] == 0.0 && (*spectra_)[charge]->length() > 0)))
	{
		//std::cerr<<"\t Setting negative only values"<<std::endl << std::endl;

 		negOnly_[charge] = negOnly(charge); //Xiuxia, update the probability in the called function and set negOnly_[charge] = 1 for this charge.
    }
  }

}


Boolean2 MixtureModel::iterate(int counter) {
	Boolean2 output = False;

	for(int charge = 0; charge < numCharge; charge++) {

		if(done_[charge] < 0) {
			if(! iterate(counter, charge)) {
				done_[charge] = counter;

				//std::cerr << std::endl << " " << counter << " iterations for " << charge+1 << "+" ;
			}
		}
		if(done_[charge] < 0) {
			output = True;
		}
	}

	return output;
}




Boolean2 MixtureModel::iterate(int counter, int charge) {//Xiuxia, counter is not used in this function
	Boolean2 output = False;

	if(priors_[charge] == 0.0)	//Xiuxia, initialized to -1.0 at the creation of the MixtureModel
		return False;

	if((*spectra_)[charge]->length() == 0) {
		updatePriors(charge);
		return False; // done
	}
	else if((*spectra_)[charge]->length() < min_num_specs_) {
		return False; // done, and will handle via negonly later
	}

	int maxcount = 100;

	//std::cout << "counter = " << counter << ", charge = " << charge << std::endl ;
	//std::cout << "priors_[" << charge << "]" << priors_[charge] << std::endl ;

	//std::cout << (*(*probs_)[1])[10] << std::endl ;

	if(updateProbs(counter, charge)) {
		output = True;
	}

	if(updateDistr(discrim_name_, charge)) {
		output = True;
	}

	if(updateDistr(ntt_name_, charge)) {
		output = True;
	}

	if(updatePriors(charge)) {
		output = True;
	}

  // do the rest here
	for(int k = 2; k < (*mixtureDistrs_)[charge]->length(); k++) {//Xiuxia, NMC
		if((*(*mixtureDistrs_)[charge])[k]->update((*probs_)[charge])) {
			output = True;
		}
	}

	return output;
}




Boolean2 MixtureModel::updateProbs(int counter, int charge) {
	Boolean2 output = False;
	float maxdiff = 0.1;

	Array<float>* newprobs = new Array<float>() ;

	for(int k = 0; k < (*probs_)[charge]->length(); k++)
	{
		float probValue = computeProb(charge, k) ;

		//std::cout << (*(*dataset_num_All_)[0])[k] << "\t" << (*(*fvalAll_)[0])[k] << "\t" << (*(*probs_)[0])[k] << std::endl ;

		newprobs->insertAtEnd(probValue);

		if(abs((*newprobs)[k] - (*(*probs_)[charge])[k]) >= maxdiff)
		{
			output = True;
		}
	}

	if(output)
	{
		assert(newprobs->length() == (*probs_)[charge]->length());
		delete (*probs_)[charge] ;
		probs_->remove(charge);
		probs_->insert(charge, newprobs);

		char fpath[80], fname[80], string_temp[10] ;
		strcpy(fpath, "out_Intermediate") ;

		/*
		if (counter > 200 && counter <= 210 && charge == 1) {

			(*(*mixtureDistrs_)[charge])[0]->printDistr() ;

			strcpy(fname, fpath) ;
			_itoa(counter, string_temp, 10);
			strcat(fname, string_temp) ;
			strcat(fname, "_charge") ;
			_itoa(charge, string_temp, 10) ;
			strcat(fname, string_temp);
			strcat(fname, ".txt");

			std::ofstream fout(fname);

			for (int k = 0 ; k < newprobs->length() ; k++) {
				fout << (*(*fvalAll_)[charge])[k] << "\t" << (*newprobs)[k] << "\n";
			}
			fout.close() ;
		}

		if (counter > 400 && counter <= 410 && charge == 1) {

			(*(*mixtureDistrs_)[charge])[0]->printDistr() ;

			strcpy(fname, fpath) ;
			_itoa(counter, string_temp, 10);
			strcat(fname, string_temp) ;
			strcat(fname, "_charge") ;
			_itoa(charge, string_temp, 10) ;
			strcat(fname, string_temp);
			strcat(fname, ".txt") ;
			std::ofstream fout(fname);

			for (int k = 0 ; k < newprobs->length() ; k++) {
				fout << (*(*fvalAll_)[charge])[k] << "\t" << (*newprobs)[k] << "\n";
			}
			fout.close() ;
		}
		*/
	}
	else
	{
		delete newprobs;
	}

	return output;
}


Boolean2 MixtureModel::updatePriors(int charge) {

  Boolean2 output = False;
  float maxdiff = 0.002;
  float next = 0;
  int tot = 0;
  for(int k = 0; k < (*probs_)[charge]->length(); k++) {
    if((*(*probs_)[charge])[k] >= 0) {
      next += (*(*probs_)[charge])[k];
	tot++;
    }
  }
  if(tot > 0) {
    next /= tot;
  }
  float diff = next - priors_[charge];
  if(priors_[charge] == -1 || abs(diff) >= maxdiff) {
    output = True;
  }
  priors_[charge] = next;
  return output;
}


int MixtureModel::getNegOnlyCharge(int charge) {
  //int alternatives[] = {1, 2, 1};
  //int secondalts[] = {2, 0, 0};
  int alternatives[] = {1, 2, 1, 1, 1};//Xiuxia, use setting for charge 3
  int secondalts[] = {2, 0, 0, 0, 0};//Xiuxia, use setting for charge 3
  if(! negOnly_[alternatives[charge]] && priors_[alternatives[charge]] != 0.0 && (*spectra_)[alternatives[charge]]->length() >= min_num_specs_) {
    return alternatives[charge];
  }
  if(! negOnly_[secondalts[charge]] && priors_[secondalts[charge]] != 0.0 && (*spectra_)[secondalts[charge]]->length() >= min_num_specs_) {
    return secondalts[charge];
  }
  return -1;
}



Boolean2 MixtureModel::isNumber(char c) {
  return c >= '0' && c <= '9';
}

Boolean2 MixtureModel::negOnly(int charge) {
	float min_prob = 0.9; // below that, call a '0'
	Boolean2 output = False;
	float negonlyprior = 0.1; // default setting

	int alt = getNegOnlyCharge(charge);

	if(alt < 0) {

		// use training data only....
		((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->reset(); // use training settings
	}
	else {

		if((*spectra_)[charge]->length() < 15) {
			// use training data anyway.
			((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->reset(); // use training settings
		}

		for(int k = 0; k < (*mixtureDistrs_)[charge]->length(); k++) {
			if((*(*mixtureDistrs_)[charge])[k]->negOnly()) {
				(*(*mixtureDistrs_)[charge])[k]->setPosDistr((*(*mixtureDistrs_)[alt])[k]);

			}
		} // next distr

		negonlyprior = priors_[alt];
	}

	output = True;

  // now compute rough estimated probs using the new distribution and negonlypriors....
	float posprob;
	float negprob;
	float prob;

	for(int d = 0; d < (*probs_)[charge]->length(); d++) {
		prob = ((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->getRightCumulativeNegProb(1.0, d, 10.0);
		posprob = (1.0 - prob) * negonlyprior;
		negprob = prob * (1.0 - negonlyprior);

		if(alt >= 0) {
			for(int k = 0; k < (*mixtureDistrs_)[charge]->length(); k++) {
				if((*(*mixtureDistrs_)[charge])[k]->negOnly()) {
					posprob *= (*(*mixtureDistrs_)[charge])[k]->getPosProb(d);
					negprob *= (*(*mixtureDistrs_)[charge])[k]->getNegProb(d);
				}
			}
		}

		if(posprob == 0.0) {
			(*probs_)[charge]->replace(d, 0.0);
		}
		else {
			prob = posprob / (posprob + negprob);

			if(prob >= min_prob) {
				prob = -1.0 * (charge) - 1.0;
			}
			else {
				prob = 0.0;
			}

			(*probs_)[charge]->replace(d, prob);
		}
	} // next data member

	priors_[charge] = 0.0; // set this to 0

	((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->resetTot(); // set prior to 0

	if(charge == 1 || charge == 2) {
		use_adj_probs_ = False;
	}

	return output; //Xiuxia, output = True, always, in this piece of code
}

float MixtureModel::computeProb(int charge, int index) {	// Xiuxia, index is the index to the spectra for this particular charge

	assert(charge < probs_->length());	//Xiuxia, only 3 charges are considered and probs is an array of 3 pointers
	assert(index < (*probs_)[charge]->length());	//Xiuxia, length() gives the number of spectra for this charge

	float posprob = 1.0;
	float negprob = 1.0;
	if(priors_ != NULL && priors_[charge] >= 0) { //Xiuxia, priors_ is initializaed to -1.0
		posprob *= priors_[charge];
		negprob *= (1.0 - priors_[charge]);
	}
	for(int k = 0; k < (*mixtureDistrs_)[charge]->length(); k++) {

		assert(index < (*(*mixtureDistrs_)[charge])[k]->getNumVals());
		if(index < 0) {
		}
		posprob *= (*(*mixtureDistrs_)[charge])[k]->getPosProb(index);
		negprob *= (*(*mixtureDistrs_)[charge])[k]->getNegProb(index);
	}

	if(posprob + negprob == 0.0) {
		return 0.0;
	}
	return (posprob / (posprob + negprob));	// Xiuxia, Bayesian formula
}

// this must be modified for maldi case (must substitute DiscrimValMixtureDistr values for positive as well)
float MixtureModel::computeProbWithNTT(int charge, int orig_ntt, float orig_prob, int ntt, int index = -1) {
	// index is only used when considering maldi case. But Why is tt set equal to -1 here?

  assert(charge < probs_->length());

  if(orig_ntt == ntt || orig_prob <= 0.0 || negOnly_[charge]) {
    return orig_prob;
  }

  float posntt1 = ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_, charge)))->getNTTPosFraction(orig_ntt);
  float negntt1 = ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_, charge)))->getNTTNegFraction(orig_ntt);
  if(maldi_) {
    posntt1 *= ((DecayContinuousMultimixtureDistr*)(((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->getPosDistr()))->getMixtureProbWithNTT(orig_ntt, index);
    negntt1 *= ((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->getNegProb(index);
  } // maldi case
  float ratio1;
  if(posntt1 == 0.0) {
    ratio1 = 999.0;
  }
  else {
    ratio1 = negntt1 / posntt1;
  }
  float posntt2 = ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_, charge)))->getNTTPosFraction(ntt);
  float negntt2 = ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_, charge)))->getNTTNegFraction(ntt);
  if(maldi_) {
    posntt2 *= ((DecayContinuousMultimixtureDistr*)(((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->getPosDistr()))->getMixtureProbWithNTT(ntt, index);
    negntt2 *= ((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->getNegProb(index);
  } // maldi case
  float ratio2;
  if(posntt2 == 0.0 && negntt2 == 0.0) {
    return 0.0;
  }

  if(posntt2 == 0.0) {
    ratio2 = 999.0;
  }
  else {
    ratio2 = negntt2 / posntt2;
  }

  float denom;
  if(ratio1 == 0.0) {
    denom = 999.0;
  }
  else {
    denom  = ((1.0 - orig_prob) * ratio2 / (orig_prob * ratio1)) + 1.0;
  }
  return 1 / denom;
}

float MixtureModel::getProb(int charge, int index) {
  return (*(*probs_)[charge])[index];
}


Boolean2 MixtureModel::enterDistribution(int charge, int index, const char* tag, char* value) {
  Boolean2 output = False;
  int numDistributions = (*mixtureDistrs_)[charge]->length() ;
  for(int k = 0; k < numDistributions ; k++)
  {
	MixtureDistr *mixDistr = (*(*mixtureDistrs_)[charge])[k] ;
	char *distrTag = mixDistr->getTag() ;
    if(strcmp(distrTag, tag) == 0)
	{
      mixDistr->enter(index, value); //index is not used in the function
      output = True;
    }
  }
  return output;
}


char* MixtureModel::getTagValue(char* data, char*tag) {
  char* output = NULL;
  char* result = strstr(data, tag);
  if(result == NULL)
    return output;
  int k = strlen(tag);
  while(k < strlen(result) && result[k] != ' ' && result[k] != '\t' &&
	result[k] != '\n')
    k++;
  output = new char[k - strlen(tag) + 1];
  strncpy(output, result + strlen(tag), k - strlen(tag));
  output[k - strlen(tag)] = 0;
  return output;
}


void MixtureModel::readData(std::vector<SequestResult> &results)
{
  char tag[75];
  char value[100];
  char spec[100];
  int charge;
  int index = 0; //-1;
  int nextpair;
  char next_line[200];


  int numResults = results.size() ;
  char *decoySpecName = "spec" ;
  char *nextSpec ;

  const char *nmcTag = "nmc" ;
  const char *discrTag = "fval" ;
  const char *nttTag = "ntt" ;
  char tagValue[128] ;

  for (int resultNum = 0 ; resultNum < numResults ; resultNum++)
  {
	  SequestResult sequestResult = results[resultNum] ;
	  // for our purposes we are just going to put the same charater in the spectra field everywhere.
	  // since we won't link to a html file.
	  nextSpec = new char[strlen(decoySpecName)+1];
	  strcpy(nextSpec, decoySpecName);
	  nextSpec[strlen(decoySpecName)] = 0;

	  int currentCharge = sequestResult.charge_ - 1;
	  (*spectra_)[currentCharge]->insertAtEnd(nextSpec);
	  (*probs_)[currentCharge]->insertAtEnd(0.5);
	  (*inds_)[currentCharge]->insertAtEnd(index++);



	  // add values for fval, ntt, nmc (F, # tryptic termini and # missed cleavages)
	  // while I am passing the wrong index (i.e. the index of the next value), I am doing so
	  // because that is how the original idiots that wrote PeptideProphet wrote it.
	  int currentNtt = sequestResult.mint_numTT ;
	  int currentNmc = sequestResult.mint_numMissedCleavages ;

	  int currentDatasetNum = sequestResult.dataset_num_ ; //Xiuxia
	  float currentFval = (float) sequestResult.mdbl_discriminantScore ; //Xiuxia
	  int currentScanNumber = sequestResult.ScanNumber ;	//Xiuxia
	  float currentXCorr = sequestResult.xcorr_ ;	//Xiuxia
	  float currentDeltaCn2 = sequestResult.delta_ ; // Xiuxia

	  (*dataset_num_All_)[currentCharge]->insertAtEnd(currentDatasetNum) ; //Xiuxia
	  (*scanNumberAll_)[currentCharge]->insertAtEnd(currentScanNumber);		//Xiuxia
	  (*xcorrAll_)[currentCharge]->insertAtEnd(currentXCorr);	//Xiuxia
	  (*fvalAll_)[currentCharge]->insertAtEnd(currentFval) ;	//Xiuxia
	  (*deltaCn2All_)[currentCharge]->insertAtEnd(currentDeltaCn2) ; //Xiuxia

	  _itoa(currentNtt, tagValue,10) ;
	  enterDistribution(currentCharge, index, nttTag, tagValue);

	  _itoa(currentNmc, tagValue,10) ;
	  enterDistribution(currentCharge, index, nmcTag, tagValue);

	  sprintf(tagValue, "%f", currentFval) ;
	  enterDistribution(currentCharge, index, discrTag, tagValue);

  }

  //std::cerr << " Processed: " << std::endl ;
	for (int charge = 0; charge < numCharge; charge++){
		//std::cerr << " " << (*spectra_)[charge]->length() << "\t" << charge+1 <<  "+ " << std::endl ;
  }
}

void MixtureModel::writeResultsInOrder(char* filename) {
	if(use_adj_probs_ && pairedspecs_ == NULL) {
		computeAdjDoubleTriplySpectraProbs();
	}

	int ordered_ind = 0;
	int charge;
	OrderedResult* nextresult;
	int total = (*spectra_)[0]->length() + (*spectra_)[1]->length() + (*spectra_)[2]->length();
	OrderedResult* ordered = new OrderedResult[total];

  // first the 1+
	for(int k = 0; k < (*spectra_)[0]->length(); k++) {
		nextresult = new OrderedResult(0, k, (*(*inds_)[0])[k]);
		ordered[ordered_ind++] = *nextresult;
	}
	if(use_adj_probs_) {
		for(int k = 0; k < num_pairs_; k++) {
			charge = atoi(pairedspecs_[k].name_ + strlen(pairedspecs_[k].name_)-1);
			nextresult = new OrderedResult(charge-1, k, pairedspecs_[k].ind_);
			ordered[ordered_ind++] = *nextresult;
		}
	}
	else {
		for(int ch = 1; ch < numCharge; ch++) {
			for(int k = 0; k < (*spectra_)[ch]->length(); k++) {
				nextresult = new OrderedResult(ch, k, (*(*inds_)[ch])[k]);
				ordered[ordered_ind++] = *nextresult;
			}
		}
	}

  // now order them
	qsort(ordered, total, sizeof(OrderedResult), (int(*)(const void*, const void*)) comp_ords);

	std::ofstream fout(filename);
	if(! fout) {
		//std::cerr << "could not open filename" << std::endl;
		std::stringstream str;
		str << "could not open filename" << std::endl;
		throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
		//exit(1);
	}

	for(int k = 0; k < total; k++)
	{
		int originalIndex = ordered[k].index_ ;
		int currentCharge = ordered[k].charge_ ;

		int temp = (*(*inds_)[currentCharge])[30000] ;
		int originalVectorIndex = (*(*inds_)[currentCharge])[originalIndex] ;
		//Xiuxia, error here. For charge = 2, (*(*inds_))[2]->length() = 30843, but originalIndex = 37755.

		float peptideProphetPValue = 0 ;
		bool negativeOnly = false ;
		bool currentSpectraAdjusted = false ;
		bool adjusted = false ;

		if(currentCharge == 0)
		{
			if(! negOnly_[0] && priors_[0] == 0.0)
			{
			  // no data
				peptideProphetPValue = 0 ;
			}
			else
				peptideProphetPValue = (*(*probs_)[0])[originalIndex];
			if(negOnly_[0])
			{
				negativeOnly = true ;
			}
		}
		else
		{
			if(use_adj_probs_)
			{
			  //remember that if paired, the paired information would be in pairedspecs_[originalIndex]
				peptideProphetPValue = pairedspecs_[originalIndex].prob_;
				if(isAdjusted(originalIndex))
				{
					adjusted = true ;
				}
			}
			else
			{
				if(! negOnly_[currentCharge] && priors_[currentCharge] == 0.0)
				{ // no data
					peptideProphetPValue = 0 ;
				}
				else
				{
					peptideProphetPValue = (*(*probs_)[currentCharge])[originalIndex];
				}
				if(negOnly_[currentCharge])
				{
					negativeOnly = true ;
				}
			}
		}
//	  fout<<currentCharge<<","<<sequestResults[originalVectorIndex].mdbl_discriminantScore<<","<<sequestResults[originalVectorIndex].massdiff_<<","<<sequestResults[originalVectorIndex].xcorr_<<","<<sequestResults[originalVectorIndex].delta_<<","<<sequestResults[originalVectorIndex].peptide_<<","<<sequestResults[originalVectorIndex].protein_<<std::endl ;
	} // charge 2,3

	fout.close();
}

// write out probabilities to output file (input to mixture_aft.pl)
//void MixtureModel::writeResults(std::vector<SequestResult> &results, char* filename)
void MixtureModel::writeResults(char* filename)
{
  if(
	  use_adj_probs_ && pairedspecs_ == NULL)
  {
	computeAdjDoubleTriplySpectraProbs();
  }

  std::ofstream fout(filename);
  if(! fout) {
    //std::cerr << "could not open filename" << std::endl;
	  std::stringstream str;
	  str << "could not open filename" << std::endl;
	  throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
    //exit(1);
  }


  int ind1, ind2 ;
  // header
  //fout << "HitNum" << "\t" << "FScore" << "\t" << "charge" << "\t" << "NTT" << "\t" << "NMC" << "\t" << "Probability" << "\t" << "negOnly" << std::endl ;
  fout << "HitNum" << "\t" << "FScore" << "\t" << "Probability" << "\t" << "negOnly" << std::endl ;

  // first the 1+
	//std::cout << " negOnly_[" << 0 << "] = " << negOnly_[0] << std::endl ;

	for(int k = 0; k < (*spectra_)[0]->length(); k++) {
		int temp = 0;

		if(! negOnly_[0] && priors_[0] == 0.0) { // definite 0
			fout << "0 \t nodata";
		}
		else {

			int temp = 1 ;	//charge

			//fout << (*(*dataset_num_All_)[0])[k] << "\t" << (*(*scanNumberAll_)[0])[k] << "\t" << temp << "\t" << (*(*xcorrAll_)[0])[k] << "\t" << (*(*deltaCn2All_)[0])[k] << "\t" << (*(*fvalAll_)[0])[k] << "\t" << ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_,0)))->getNTTValue(k) << "\t" << ((NMCMixtureDistr*)(getMixtureDistr("no. missed cl",0)))->getNMCValue(k) << "\t" << (*(*probs_)[0])[k] ; //Xiuxia<< (*(*probs_)[0])[k] ; //Xiuxia
			//fout << (*(*dataset_num_All_)[0])[k] << "\t" << (*(*fvalAll_)[0])[k] << "\t" << temp << "\t" << ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_,0)))->getNTTValue(k) << "\t" << ((NMCMixtureDistr*)(getMixtureDistr("no. missed cl",0)))->getNMCValue(k) << "\t" << (*(*probs_)[0])[k] ;
			fout << (*(*dataset_num_All_)[0])[k] << "\t" << (*(*fvalAll_)[0])[k] << "\t" << (*(*probs_)[0])[k] ; //Xiuxia<< (*(*probs_)[0])[k] ; //Xiuxia
		}

		// for probabilities with other NTT, commented out for the pipeline
		//fout <<  "\t(";
		//fout <<  "\t";//Xiuxia
		//for(int n = 0; n < 3; n++) {
		//	fout << computeProbWithNTT(0, ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_, 0)))->getNTTValue(k), (*(*probs_)[0])[k], n, k);
		//	if(n < 2) {
		//	fout << "\t";
		//	}
		//}
		//fout << ")";//Xiuxia, commented out

		if(negOnly_[0]) {
			//fout << "\t negonly";
			ind1 = 1 ;
			fout << "\t" << ind1 ; //Xiuxia
		}
		else {
			ind1 = 0 ;
			fout << "\t" << ind1 ;
		}
		fout << "\n";
	}

	if(use_adj_probs_) {
		for(int k = 0; k < num_pairs_; k++) {
			fout << pairedspecs_[k].name_ << "\t" << pairedspecs_[k].prob_;

			fout << "\t(";
			for(int n = 0; n < 3; n++) {
				fout << computeProbWithNTT(
				   atoi(pairedspecs_[k].name_+strlen(pairedspecs_[k].name_)-1)-1,
				   pairedspecs_[k].ntt_,
				   pairedspecs_[k].prob_, n, pairedspecs_[k].ind_);
				if(n < 2) {
					fout << ",";
				}
			}
			fout << ")";
			if(isAdjusted(k)) {
				fout << "\tadj";
			}
			fout << std::endl;
		}
	}
	else {
		for(int charge = 1; charge < numCharge; charge++) {

			//std::cout << " negOnly_[" << charge << "] = " << negOnly_[charge] << std::endl ;

			for(int k = 0; k < (*spectra_)[charge]->length(); k++) {

				if(! negOnly_[charge] && priors_[charge] == 0.0) { // definite 0
					fout << "0 \t nodata";
				}
				else {
					//fout << (*(*dataset_num_All_)[charge])[k] << "\t" << (*(*scanNumberAll_)[charge])[k] << "\t" << charge+1 << "\t" << (*(*xcorrAll_)[charge])[k] << "\t" << (*(*deltaCn2All_)[charge])[k] << "\t" << (*(*fvalAll_)[charge])[k] << "\t" << ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_,charge)))->getNTTValue(k) << "\t" << ((NMCMixtureDistr*)(getMixtureDistr("no. missed cl",charge)))->getNMCValue(k) << "\t" << (*(*probs_)[charge])[k] ; //Xiuxia
					//fout << (*(*dataset_num_All_)[charge])[k] << "\t" << (*(*fvalAll_)[charge])[k] << "\t" << charge+1 << "\t" << ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_,charge)))->getNTTValue(k) << "\t" << ((NMCMixtureDistr*)(getMixtureDistr("no. missed cl",charge)))->getNMCValue(k) << "\t" << (*(*probs_)[charge])[k] ; //Xiuxia
					fout << (*(*dataset_num_All_)[charge])[k] << "\t" << (*(*fvalAll_)[charge])[k] << "\t" << (*(*probs_)[charge])[k] ; //Xiuxia

				}

				// for probabilities with other NTT, commented out for the pipeline
				//fout <<  "\t(";
				//fout <<  "\t";

				//for(int n = 0; n < 3; n++) {//Xiuxia, n = ntt

				//	fout << computeProbWithNTT(charge, ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_, charge)))->getNTTValue(k), (*(*probs_)[charge])[k], n, k);

				//	if(n < 2) {
				//		//fout << ",";
				//		fout << "\t";
				//	}
				//}
				//fout << ")";//Xiuxia, commented out

				if(negOnly_[charge]) {
					//fout << "\t negonly";
					ind2 = 1 ;
					fout << "\t" << ind2 ;
				}
				else {
					ind2 = 0 ;
					fout << "\t" << ind2 ;
				}
				fout << "\n";
			}	// next spectrum

		} // next charge
	} // else
	fout.close();
}







// Xiuxia, added this function on 08/02/2006
void MixtureModel::writeResultsOrdered(const char* filename, std::vector<DatasetNumMap> &vecDatasetMap)
{

	if(use_adj_probs_ && pairedspecs_ == NULL)
	{
		computeAdjDoubleTriplySpectraProbs();
	}

	OutputContent oneOutputUnit ;
	std::vector<OutputContent> vecOutputContent ;

	for(int k = 0; k < (*spectra_)[0]->length(); k++) {
		int temp = 0;

		if(! negOnly_[0] && priors_[0] == 0.0) { // definite 0
			oneOutputUnit.HitNum = (*(*dataset_num_All_)[0])[k] ;
			oneOutputUnit.fscore = 0.0 ;
			oneOutputUnit.prob = 0.0 ;
			oneOutputUnit.negonly = 1 ;

			vecOutputContent.push_back(oneOutputUnit) ;
		}
		else {
			oneOutputUnit.HitNum = (*(*dataset_num_All_)[0])[k] ;
			oneOutputUnit.fscore = (*(*fvalAll_)[0])[k] ;
			oneOutputUnit.prob = (*(*probs_)[0])[k] ;

			if(negOnly_[0]) {
				//fout << "\t negonly";
				oneOutputUnit.negonly = 1 ;
			}
			else {
				oneOutputUnit.negonly = 0;
			}

			vecOutputContent.push_back (oneOutputUnit) ;
		}
	}


	for(int charge = 1; charge < numCharge; charge++) {

		for(int k = 0; k < (*spectra_)[charge]->length(); k++) {

			if(! negOnly_[charge] && priors_[charge] == 0.0) { // definite 0
				oneOutputUnit.HitNum = (*(*dataset_num_All_)[charge])[k] ;
				oneOutputUnit.fscore = 0.0 ;
				oneOutputUnit.prob = 0.0 ;
				oneOutputUnit.negonly = 1 ;

				vecOutputContent.push_back(oneOutputUnit) ;
			}
			else {
				oneOutputUnit.HitNum = (*(*dataset_num_All_)[charge])[k] ;
				oneOutputUnit.fscore = (*(*fvalAll_)[charge])[k] ;
				oneOutputUnit.prob = (*(*probs_)[charge])[k] ;

				if(negOnly_[0]) {
					oneOutputUnit.negonly = 1 ;
				}
				else {
					oneOutputUnit.negonly = 0;
				}

				vecOutputContent.push_back (oneOutputUnit) ;
			}
		}	// next spectrum

	} // next charge

	sort(vecOutputContent.begin(), vecOutputContent.end(), MixtureModel::SortOutputResultsByHitnum) ;

//	std::ofstream fout(filename);
	FILE *fp = fopen(filename, "w") ;

	if(fp == NULL) {
		std::cerr << "could not open filename" << std::endl;
		exit(1);
	}

//	fout << "HitNum" << "\t" << "FScore" << "\t" << "Probability" << "\t" << "negOnly" << std::endl ;
	fprintf(fp, "HitNum\tFScore\tProbability\tnegOnly\n") ;

	int indexMap = 0 ;
	int totaloutrow = 0 ;
	const int MAXSIZE = 16 ;
	char tempBuff[MAXSIZE];

	for (int index = 0; index < vecOutputContent.size(); index++ )
	{
		while (vecOutputContent[index].HitNum == vecDatasetMap[indexMap].dataset_num_start)
		{
			//fout << vecDatasetMap[indexMap].dataset_num_other << "\t" << vecOutputContent[index].fscore << "\t" << vecOutputContent[index].prob << "\t" << vecOutputContent[index].negonly << "\n" ;
			int datasetNum = vecDatasetMap[indexMap].dataset_num_other ;
			float fScore = vecOutputContent[index].fscore ;
			float prob = vecOutputContent[index].prob ;
			int negOnly = vecOutputContent[index].negonly ;

			sprintf(tempBuff, "%.4f", prob);

			if (strcmp(tempBuff, "0.0000") == 0 || strcmp(tempBuff, "1.0000") == 0)
				fprintf(fp, "%d\t%.4f\t%.0f\t%d\n", datasetNum, fScore, prob, negOnly) ;
			else
				fprintf(fp, "%d\t%.4f\t%.4f\t%d\n", datasetNum, fScore, prob, negOnly) ;

			indexMap++ ;
			totaloutrow++ ;
		}

	}

	fclose(fp) ;
//	fout.close();
	std::cout << "  ... wrote out " << totaloutrow << " lines" << std::endl ;
}






bool MixtureModel::SortOutputResultsByHitnum(OutputContent &a, OutputContent &b)
{
	if (a.HitNum < b.HitNum)
		return true ;
	if (a.HitNum > b.HitNum)
		return false ;

	return false ;
}








Boolean2 MixtureModel::isAdjusted(int index) {
  if(adjusted_ == NULL) {
    return False;
  }
  for(int k = 0; k < adjusted_->length(); k++) {
    if((*adjusted_)[k] == index) {
      return True;
    }
    else if((*adjusted_)[k] > index) {
      return False;
    }
  }
  return False;
}

void MixtureModel::writeModelfile(char* modelfile) {
  writeDistr(modelfile);
  computeEstimatedSensAndError(modelfile);
  writeDiscrimValDistr(modelfile);
}

// model learned distribution summaries
void MixtureModel::writeDiscrimValDistr(char* filename) {
  float minval = -5.0;
  float maxval = 10.0;
  float window = 0.2;
  FILE* fout;
  if((fout = fopen(filename, "a")) == NULL) {
    //std::cerr << "could not open " << filename << std::endl;
	  std::stringstream str;
	  str  << "could not open " << filename << std::endl;
	  throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
    //exit(1);
  }
  fprintf(fout, "\n\nOBSERVED AND MODEL DISTRIBUTIONS\n");
  fprintf(fout, "#Fval\t1+ distr\t1+ pos\t1+neg\t2+ distr\t2+ pos\t2+ neg\t3+ distr\t3+ pos\t3+ neg\n");
  int numwins = int ((maxval - minval)/window);
  for(int k = 0; k < numwins; k++) {
    fprintf(fout, "%0.2f\t", (k+0.5)*window + minval);
    for(int charge = 0; charge < 3; charge++) {
      fprintf(fout, "%d\t%0.2f\t%0.2f", ((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->slice(k*window + minval, (k+1)*window + minval),((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->posSlice(k*window + minval, (k+1)*window + minval), ((DiscrimValMixtureDistr*)(
		  (discrim_name_, charge)))->negSlice(k*window + minval, (k+1)*window + minval));
      if(charge < 2) {
	fprintf(fout, "\t");
      }
    } // next charge
    fprintf(fout, "\n");
  }

  if(maldi_ || qtof_) {
    fprintf(fout, "\n\nMODEL DISTRIBUTIONS BY NTT VALUE\n");
    if(maldi_)
    fprintf(fout, "#Fval\tntt0 1+ pos\tntt1 1+ pos\tntt2 1+ pos\n");
    else
      fprintf(fout, "#Fval\tntt0 1+ pos\tntt0 2+ pos\tntt0 3+ pos\tntt1 1+ pos\tntt1 2+ pos\tntt1 3+ pos\tntt2 1+ pos\tntt2 2+ pos\tntt2 3+ pos\n");
    for(int k = 0; k < numwins; k++) {
      fprintf(fout, "%0.2f", (k+0.5)*window + minval);
      for(int ntt = 0; ntt < 3; ntt++) {
	for(int charge = 0; charge < 3; charge++) {
	  if(qtof_ || (maldi_ && charge == 0))
	    fprintf(fout, "\t%0.2f", ((DiscrimValMixtureDistr*)(getMixtureDistr(discrim_name_, charge)))->posSliceWithNTT(k*window + minval, (k+1)*window + minval, ntt));
	  //if(qtof_ && charge < 2) {
	  //fprintf(fout, "\t");
	  //}
	} // next charge
	//fprintf(fout, "\n");
      } // next ntt
      fprintf(fout, "\n");
    } // next win

  } // if maldi or qtof


  fclose(fout);
}

void MixtureModel::printDistr() {
  for(int charge = 0; charge < 3; charge++) {
    printf("charge %d after %d iterations:\n", charge+1, done_[charge]);
    printf("\tprior: %0.3f\n", priors_[charge]);
    for(int k = 0; k < (*mixtureDistrs_)[charge]->length(); k++) {
      (*(*mixtureDistrs_)[charge])[k]->printDistr();
    }
    printf("\n\n");
  }
}

void MixtureModel::writeDistr(char* filename) {
  FILE* fout;
  if((fout = fopen(filename, "w")) == NULL) {
    //std::cerr << "cannot open " << filename << std::endl;
    //exit(1);
	  std::stringstream str;
	  str << "cannot open " << filename << std::endl;
	  throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
  }
  fprintf(fout, "PeptideProphet v. 1.0  Andrew Keller  ISB\n\n\n");


  for(int charge = 0; charge < 3; charge++) {
    fprintf(fout, "FINAL %d+ MODEL after %d iterations:\n", charge+1, done_[charge]);
    fprintf(fout, "number of spectra: %d\n", (*spectra_)[charge]->length());
    if(negOnly_[charge]) {
      if(getNegOnlyCharge(charge) < 0)
	fprintf(fout, "using training data positive distributions to identify candidates ('-%d') above background ('0')\n", charge+1);
      else
	fprintf(fout, "using %d+ positive distributions to identify candidates ('-%d') above background ('0')\n", getNegOnlyCharge(charge)+1, charge+1);
    }
    else if(pseudonegs_[charge]) {
      fprintf(fout, "using %s 0 data as pseudonegatives\n", ntt_name_);
    }
    else {
      fprintf(fout, "using training data negative distributions\n");
    }

    fprintf(fout, "\tprior: %0.3f, total: %0.1f\n", priors_[charge], getTotalProb(charge));
    for(int k = 0; k < (*mixtureDistrs_)[charge]->length(); k++) {
      (*(*mixtureDistrs_)[charge])[k]->writeDistr(fout);
    }
   fprintf(fout, "\n");
  }
  fclose(fout);
}


// add constraint that final probs for 2+ and 3+ precursor ion interpretations
// of same spectrum cannot total more than unity
void MixtureModel::computeAdjDoubleTriplySpectraProbs() {
	//std::cerr << negOnly_[1] << "\t" << negOnly_[2] << std::endl ;

  if(negOnly_[1] || negOnly_[2]) {
    return;
  }
  num_pairs_ = 0;
  for(int charge = 1; charge < 3; charge++) {
    num_pairs_ += (*spectra_)[charge]->length();
  }


  pairedspecs_ = new Spectrum[num_pairs_];
  adjusted_ = new Array<int>;
  int total = 0;
  char* next;
  for(int charge = 1; charge < 3; charge++) {
      for(int k = 0; k < (*probs_)[charge]->length(); k++) {
	char* next = new char[strlen((*(*spectra_)[charge])[k]) + 3];
	strcpy(next, (*(*spectra_)[charge])[k]);
	strcat(next, ".");
	next[strlen((*(*spectra_)[charge])[k]) + 1] = charge + 49;
	next[strlen((*(*spectra_)[charge])[k]) + 2] = 0;

	Spectrum* spec = new Spectrum(next, (*(*probs_)[charge])[k], ((NTTMixtureDistr*)(getMixtureDistr(ntt_name_, charge)))->getNTTValue(k), (*(*inds_)[charge])[k]);
	pairedspecs_[total++] = *spec;
      }
  }
  // now order it
  qsort(pairedspecs_, num_pairs_, sizeof(Spectrum), (int(*)(const void*, const void*)) comp_specs);
  // now walk down list and make 2/3 adjustments....
  for(int k = 0; k < num_pairs_-1; k++) {
    if(strncmp(pairedspecs_[k].name_, pairedspecs_[k+1].name_, strlen(pairedspecs_[k].name_)-2) == 0) {
      // make adjustment

      float prob1 = getAdjDoublyTriplyProb(pairedspecs_[k].prob_, pairedspecs_[k+1].prob_);
      float prob2 = getAdjDoublyTriplyProb(pairedspecs_[k+1].prob_, pairedspecs_[k].prob_);
      if(pairedspecs_[k].prob_ >= 0.5 && pairedspecs_[k+1].prob_ >= 0.5) {
	adjusted_->insertAtEnd(k);
	adjusted_->insertAtEnd(k+1);
      }
      pairedspecs_[k].prob_ = prob1;
      pairedspecs_[k+1].prob_ = prob2;
      k++;
    }
  }

}

// writes out sensitivity and error for various minimum probs
void MixtureModel::computeEstimatedSensAndError(char* filename) {
  if(use_adj_probs_ && pairedspecs_ == NULL) {
    computeAdjDoubleTriplySpectraProbs();
  }
  int total = 0;
  for(int charge = 0; charge < 3; charge++) {
    if(! negOnly_[charge]) {
      total += (*probs_)[charge]->length();
    }
  }

  float* combinedprobs = new float[total];
  total = 0; // index
  for(int k = 0; k < (*probs_)[0]->length(); k++) {
    if(! negOnly_[0]) {
      combinedprobs[total++] = (*(*probs_)[0])[k];
    }
  }
  if(use_adj_probs_) {
    for(int k = 0; k < num_pairs_; k++) {
      if(! negOnly_[atoi(pairedspecs_[k].name_ + strlen(pairedspecs_[k].name_)-1)-1]) {
	combinedprobs[total++] = pairedspecs_[k].prob_;
      }
    }
  }
  else {
    for(int charge = 1; charge < 3; charge++) {
      for(int k = 0; k < (*probs_)[charge]->length(); k++) {
	if(! negOnly_[charge]) {
	  combinedprobs[total++] = (*(*probs_)[charge])[k];
	}
      }
    } // next
  } // next charge

  qsort(combinedprobs, total, sizeof(float), (int(*)(const void*, const void*)) comp_nums);

     // now sens and error as walk down list
  float thresh[] = {0.99, 0.98, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.0};
  int threshind = 0;
  float correct = 0.0;
  float incorrect = 0.0;
  float totcorrect = 0.0;
  FILE* fout;
  for(int k = 0; k < total; k++) {
    totcorrect += combinedprobs[k];
  }

  if( (fout = fopen(filename, "a")) == NULL) {
    //std::cerr << "cannot open " << filename << std::endl;
    //exit(1);
	  std::stringstream str;
	  str << "cannot open " << filename << std::endl;
	  throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
  }


  fprintf(fout, "\n\nJOINT 1+/2+/3+ MODEL SENSITIVITY/ERROR PREDICTIONS (%d tot correct)\n", (int)(getTotalProb(0) + getTotalProb(1) + getTotalProb(2)));
  fprintf(fout, "#min_prob\tsens\terr\n");
  for(int k = 0; k < total; k++) {
    if(combinedprobs[k] >= thresh[threshind]) {
      correct += combinedprobs[k];
      incorrect += 1 - combinedprobs[k];
    }
    else {
      fprintf(fout, "%0.2f\t\t%0.4f\t%0.4f\n", thresh[threshind], correct/totcorrect, incorrect/(correct + incorrect));
      threshind++;
      k--; // assay again using next threshold
    }
  } // next member
  fprintf(fout, "%0.2f\t\t%0.4f\t%0.4f\n", thresh[threshind], correct/totcorrect, incorrect/(correct + incorrect));
  fprintf(fout, "\n\n");

  fprintf(fout, "\n\nJOINT 1+/2+/3+ MIN PROB THRESHOLDS FOR SPECIFIED ERROR RATES\n");
  fprintf(fout, "#err\t\t\tmin_prob\n");
  float error_rates[] = {0.0, 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15};
  correct = 0.0;
  incorrect = 0.0;
  threshind = 0;
  for(int k = 0; k < total; k++) {
    correct += combinedprobs[k];
    incorrect += 1 - combinedprobs[k];
    if(incorrect / (correct + incorrect) > error_rates[threshind]) {
      if(k == 0) {
	fprintf(fout, "%0.3f\t\t\t%0.2f\n", error_rates[threshind], combinedprobs[k]);
      }
      else {
	fprintf(fout, "%0.3f\t\t\t%0.2f\n", error_rates[threshind], combinedprobs[k-1]);
      }
      threshind++;
      k--;
    }
    if(threshind == (sizeof(error_rates)/sizeof(float)) - 1) {
      k = total; // done
    }
  } // next ordered prob
  fprintf(fout, "\n\n");

  fclose(fout);
}

int comp_nums(const void* num1, const void* num2) {
  float* num1_ = (float*)num1;
  float* num2_ = (float*)num2;
  if(*num1_ < *num2_) return 1;
  if(*num1_ == *num2_) return 0;
  if(*num1_ > *num2_) return -1;
  return 0;
}

int comp_specs(const void* num1, const void* num2) {
  Spectrum* spec1 = (Spectrum*)num1;
  Spectrum* spec2 = (Spectrum*)num2;
  return strcmp((*spec1).name_, (*spec2).name_);
}

int comp_ords(const void* num1, const void* num2) {
  OrderedResult* res1 = (OrderedResult*)num1;
  OrderedResult* res2 = (OrderedResult*)num2;
  if((*res1).input_index_ < (*res2).input_index_) return -1;
  if((*res1).input_index_ > (*res2).input_index_) return 1;
  return 0;
}


// adjustment imposing constraint that final probs for 2+ and 3+ precursor ion interpretations
// of same spectrum do not total more than unity
float MixtureModel::getAdjDoublyTriplyProb(float prob_2_adj, float prob_of_partner) {
  if(prob_2_adj == 0 || prob_2_adj + prob_of_partner == 0) {
    return 0.0;
  }
  return (prob_2_adj * (prob_2_adj + prob_of_partner - prob_2_adj * prob_of_partner) / (prob_2_adj + prob_of_partner));
}

