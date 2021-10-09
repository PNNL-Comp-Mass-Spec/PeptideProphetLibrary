#include "DiscrimValMixtureDistr.h"

/*

Program       : DiscrimValMixtureDistr for PeptideProphet                                                       
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


DiscrimValMixtureDistr::DiscrimValMixtureDistr(int charge, char* name, char* tag, Boolean2 gamma, Boolean2 maldi, Boolean2 qtof) : MixtureDistr(charge, name, tag) {  
  all_negs_ = False;
  gamma_ = gamma;
  maldi_ = maldi;
  maxdiff_ = 0.002;
  qtof_ = qtof;
  floatvals_ = new Array<float>() ;
 
  if(qtof_ || (maldi_ && charge_ == 0)) {
    //std::cerr << "instantiating decaycontmultimix distr" << std::endl;
    posdistr_ = new DecayContinuousMultimixtureDistr(maxdiff_ * 2.0, qtof_);
  }
  else
    posdistr_ = new GaussianDistribution(maxdiff_);


  if(gamma_) {
    negdistr_ = new GammaDistribution(maxdiff_);
  }
  else {
    negdistr_ = new GaussianDistribution(maxdiff_);
  }
  
  float singlyposprior[] = {2.0, 0.4};
  float doublyposprior[] = {4.102, 1.64};
  float triplyposprior[] = {4.563, 1.84};
  float quadruplyposprior[] = {4.563, 1.84}; //Xiuxia
  float quintuplyposprior[] = {4.563, 1.84};//Xiuxia

  float singlynegprior[] = {-2.0, 0.9};
  float doublynegprior[] = {-1.945, 0.87};
  float triplynegprior[] = {-1.563, 0.78};
  float quadruplynegprior[] = {-1.563, 0.78};//Xiuxia
  float quintuplynegprior[] = {-1.563, 0.78};//Xiuxia

  float singlyneggammprior[] = {5.0, 26.0, -5.0};
  float doublyneggammprior[] = {6.16, 39, -5.0};
  float triplyneggammprior[] = {6.06, 37, -5.0};
  float quadruplyneggammprior[] = {6.06, 37, -5.0};//Xiuxia
  float quintuplyneggammprior[] = {6.06, 37, -5.0};//Xiuxia

  float singlyposmaldiprior[] = {5.0, 0.6};
  float singlynegmaldiprior[] = {1.0, 1.5};
  float singlyneggammmaldiprior[] = {6.5, 25.6, -5.52};

  //float minvals[] = {-2.0, -5.0, -5.0 };
  float minvals[] = {-2.0, -5.0, -5.0, -5.0, -5.0 };//Xiuxia

  negmean_ = 0.0;
  MIN_NUM_PSEUDOS_ = 50;
  ZERO_SET_ = 100; // ?
  NUM_DEVS_ = 6;
  USE_TR_NEG_DISTR_ = False;
  posinit_ = NULL;
  neginit_ = NULL;

  if(charge == 0) {
    if(maldi_ || qtof_) {
      ; //posinit_ = copy(singlyposmaldiprior, 2);
    }
    else {
      posinit_ = copy(singlyposprior, 2);
    }
    if(gamma) {
      if(maldi_) {
		neginit_ = copy(singlyneggammmaldiprior, 3);
      }
      else {
		neginit_ = copy(singlyneggammprior, 3);
      }
    }
    else { // not gamma
      if(maldi_) {
		neginit_ = copy(singlynegmaldiprior, 2);
      }
      else {
		neginit_ = copy(singlynegprior, 2);
      }
    }
  }
  else if(charge == 1) {
    if(! qtof_)
      posinit_ = copy(doublyposprior, 2);
    if(gamma) {
      neginit_ = copy(doublyneggammprior, 3);
    }
    else {
      neginit_ = copy(doublynegprior, 2);
    }
  }
  else if(charge == 2) {
    if(! qtof_)
      posinit_ = copy(triplyposprior, 2);
    if(gamma) {
      neginit_ = copy(triplyneggammprior, 3);
    }
    else {
      neginit_ = copy(triplynegprior, 2);
    }
  }
  else if(charge == 3) {
    if(! qtof_)
      posinit_ = copy(quadruplyposprior, 2);
    if(gamma) {
      neginit_ = copy(quadruplyneggammprior, 3);
    }
    else {
      neginit_ = copy(quadruplynegprior, 2);
    }
  }
  else if(charge == 4) {
    if(! qtof_)
      posinit_ = copy(quintuplyposprior, 2);
    if(gamma) {
      neginit_ = copy(quintuplyneggammprior, 3);
    }
    else {
      neginit_ = copy(quintuplynegprior, 2);
    }
  }

  reset();

  if(gamma_) {
    minval_ = minvals[charge];
    negdistr_->setMinVal(minval_);
  }

}


float* DiscrimValMixtureDistr::copy(float* init, int num) {
  float* output = new float[num];
  for(int k = 0; k < num; k++) {
    output[k] = init[k];
  }
  return output;
}

void DiscrimValMixtureDistr::reset() {
  if(qtof_ || (maldi_ && charge_ == 0)) {
    ((DecayContinuousMultimixtureDistr*)(posdistr_))->reset();
  }
  else {
    posdistr_->init(posinit_);
  }
  negdistr_->init(neginit_);
} 

void DiscrimValMixtureDistr::resetTot() {
  ((ContinuousDistribution*)(posdistr_))->resetTot();
  ((ContinuousDistribution*)(negdistr_))->resetTot();
} 

Boolean2 DiscrimValMixtureDistr::noDistr() {
  float num_stdevs = 0.45; //0.15;
  Boolean2 first;
  Boolean2 second;
  if(qtof_ || (charge_ == 0 && maldi_)) {
    return False;
    first = gamma_ && ((DecayContinuousMultimixtureDistr*)(posdistr_))->getMean() + ((DecayContinuousMultimixtureDistr*)(posdistr_))->getStdev() < ((ContinuousDistribution*)(negdistr_))->getMean() + ((GammaDistribution*)(negdistr_))->getZero() + num_stdevs * sqrt(((ContinuousDistribution*)(negdistr_))->getStdev());
  Boolean2 second = !gamma_ && ((DecayContinuousMultimixtureDistr*)(posdistr_))->getMean() + ((DecayContinuousMultimixtureDistr*)(posdistr_))->getStdev() < ((ContinuousDistribution*)(negdistr_))->getMean() + ((ContinuousDistribution*)(negdistr_))->getStdev();
  }
  else {
    first = gamma_ && ((ContinuousDistribution*)(posdistr_))->getMean() + ((ContinuousDistribution*)(posdistr_))->getStdev() < ((ContinuousDistribution*)(negdistr_))->getMean() + ((GammaDistribution*)(negdistr_))->getZero() + num_stdevs * sqrt(((ContinuousDistribution*)(negdistr_))->getStdev());
    second = !gamma_ && ((ContinuousDistribution*)(posdistr_))->getMean() + ((ContinuousDistribution*)(posdistr_))->getStdev() < ((ContinuousDistribution*)(negdistr_))->getMean() + ((ContinuousDistribution*)(negdistr_))->getStdev();

	//std::cerr << "mean = " << ((ContinuousDistribution*)(posdistr_))->getMean() << "\t" << "stdev = " << ((ContinuousDistribution*)(posdistr_))->getStdev() << std::endl ;
  }


  return (first || second);


			       /*
  Boolean2 first = gamma_ && ((ContinuousDistribution*)(posdistr_))->getMean() + ((ContinuousDistribution*)(posdistr_))->getStdev() < ((ContinuousDistribution*)(negdistr_))->getMean() + ((GammaDistribution*)(negdistr_))->getZero() + num_stdevs * sqrt(((ContinuousDistribution*)(negdistr_))->getStdev());
  Boolean2 second = !gamma_ && ((ContinuousDistribution*)(posdistr_))->getMean() + ((ContinuousDistribution*)(posdistr_))->getStdev() < ((ContinuousDistribution*)(negdistr_))->getMean() + ((ContinuousDistribution*)(negdistr_))->getStdev();
  return (first || second);
			       */
}
//float DiscrimValMixtureDistr::getPosProb(int index) {
float DiscrimValMixtureDistr::getPosProb(int index) {//Xiuxia
	//std::cerr << ((GaussianDistribution*)(posdistr_))->mean_ << std::endl ;//Xiuxia

	float temp; //Xiuxia

	//DiscrimValMixtureDistr::printDistr() ;
	//std::cerr << "fval = " << (*floatvals_)[index] << "\t" << "prob = " ;

	if(qtof_ || (maldi_ && charge_ == 0)) {
    //std::cerr << "pos prob: " << ((DecayContinuousMultimixtureDistr*)(posdistr_))->getMixtureProb(index, (*floatvals_)[index]) << std::endl;
    //if(((DecayContinuousMultimixtureDistr*)(posdistr_))->oneProb((*floatvals_)[index])) {
    // return 1.0;
    //}

    return ((DecayContinuousMultimixtureDistr*)(posdistr_))->getMixtureProb(index, (*floatvals_)[index]);
   
    //std::cerr << "MALDI" << std::endl;
    //std::cerr << "number of distrs in DISCRIMVAL: " << ((DecayContinuousMultimixtureDistr*)(posdistr_))->getNumDistributions() << std::endl;
	}

	if(all_negs_ || (*floatvals_)[index] < negmean_) {//Xiuxia, all_negs_ is initialized to False
		temp = 0.0;
		//std::cerr << temp << std::endl ;

		return temp ;//0.0;
	}
	else if(! maldi_ && ((GaussianDistribution*)(posdistr_))->oneProb((*floatvals_)[index])) { //Xiuxia, check if fval > mean_ + stdev_
		if(charge_ == 0) {
		;//std::cout << "one probs for " << (*floatvals_)[index] << std::endl;
		}
		temp = 1.0 ;
		//std::cerr << temp << std::endl ;
		return temp; //1.0;
	}
	else if(gamma_ && ((GammaDistribution*)(negdistr_))->zeroProb((*floatvals_)[index])) {
		temp = 0.0;
		//std::cerr << temp << std::endl ;
		return temp;//0.0;
	}
	temp = MixtureDistr::getPosProb(index);
	//std::cerr << temp << std::endl ;
	return temp;
}



float DiscrimValMixtureDistr::getNegProb(int index) {
  if(all_negs_ || (*floatvals_)[index] < negmean_) {
    return 1.0;
  }
  else if(! maldi_ && ! qtof_ && ((GaussianDistribution*)(posdistr_))->oneProb((*floatvals_)[index])) {
    return 0.0;
  }
  //else if(((maldi_ && charge_ == 0) || qtof_) && ((DecayContinuousMultimixtureDistr*)(posdistr_))->oneProb((*floatvals_)[index])) {
  //return 0.0;
  //}
  else if(gamma_ && ((GammaDistribution*)(negdistr_))->zeroProb((*floatvals_)[index])) {
    return 1.0;
  }
  //std::cerr << "neg prob: " << MixtureDistr::getNegProb(index) << std::endl;

  return MixtureDistr::getNegProb(index);
}


void DiscrimValMixtureDistr::enter(int index, float val) {
  MixtureDistr::enter(index, val);
  if(index == 0 || val < min_dataval_) {
    min_dataval_ = val;
  }
}


float DiscrimValMixtureDistr::getRightCumulativeNegProb(int index, float right_val) {
  return negSlice((*floatvals_)[index], right_val);
}

float DiscrimValMixtureDistr::getRightCumulativeNegProb(float total, int index, float right_val) {
  return ((ContinuousDistribution*)(negdistr_))->slice(total, (*floatvals_)[index], right_val);
}

int DiscrimValMixtureDistr::slice(float left_val, float right_val) {
  int tot = 0;
  for(int k = 0; k < floatvals_->length(); k++) {
    if((*floatvals_)[k] > left_val && (*floatvals_)[k] <= right_val) {
      tot++;
    }
  }
  return tot;
}


Boolean2 DiscrimValMixtureDistr::decayMultimixture() { return qtof_ || (charge_ == 0 && maldi_); }

float DiscrimValMixtureDistr::posSliceWithNTT(float left_val, float right_val, int ntt) {
  if(! decayMultimixture()) {
    //std::cerr << "error in posSliceWithNTT" << std::endl;
	  throw gcnew System::Exception("error in posSliceWithNTT");
    //exit(1);
  }
  return ((DecayContinuousMultimixtureDistr*)(posdistr_))->sliceWithNTT(left_val, right_val, ntt);
}

float DiscrimValMixtureDistr::posSlice(float left_val, float right_val) {
  if(decayMultimixture())
    return ((DecayContinuousMultimixtureDistr*)(posdistr_))->slice(left_val, right_val);
  return ((ContinuousDistribution*)(posdistr_))->slice(left_val, right_val);
}

float DiscrimValMixtureDistr::negSlice(float left_val, float right_val) {
  return ((ContinuousDistribution*)(negdistr_))->slice(left_val, right_val);
}

void DiscrimValMixtureDistr::setNegativeDistr(float mean, float stdev, float zero) {
  float* next = new float[3];
  next[0] = mean;
  next[1] = stdev;
  next[2] = zero;

  negdistr_->init(next);
  delete [] next;
}
void DiscrimValMixtureDistr::setPositiveDistr(float mean, float stdev) {
  float* next = new float[2];
  next[0] = mean;
  next[1] = stdev;
  posdistr_->init(next);
  delete [] next;
}
void DiscrimValMixtureDistr::printDistr() {
  MixtureDistr::printDistr();
  printf("\tnegmean: %0.2f\n", negmean_);
}

void DiscrimValMixtureDistr::writeDistr(FILE* fout) {
  MixtureDistr::writeDistr(fout);
  fprintf(fout, "\tnegmean: %0.2f\n", negmean_);
}

Boolean2 DiscrimValMixtureDistr::update(Array<float>* probs) {
  Boolean2 result = MixtureDistr::update(probs);
  if(noDistr()) {
    all_negs_ = True;
    result = True;
  }
  return result; 
}

// use data with 0 tryptic termini to initialize negative distributions (when available)
Boolean2 DiscrimValMixtureDistr::initializeNegDistribution(NTTMixtureDistr* nttdistr) {
  assert(nttdistr->getNumVals() == getNumVals());
  float mean = 0.0;
  float stdev = 0.0;
  int tot = 0;
  float totsq = 0.0;
  float zero;
  /*float posmean[] = { 2.0, 4.102, 4.563 };*/ 
  //loat posstdev[] = { 0.4, 1.64, 1.84 };
  float posmean[] = { 2.0, 4.102, 4.563, 4.563, 4.563 };//Xiuxia
  float posstdev[] = { 0.4, 1.64, 1.84, 1.84, 1.84 };//Xiuxia
  float MAX_SINGLY_NEGMEAN = 1.0;

  //float negmean_num_stds[] = {-1.0, 0.5, 0.5}; // by charge
  float negmean_num_stds[] = {-1.0, 0.5, 0.5, 0.5, 0.5};//Xiuxia
  float maldi_negmean_num_stds = 0.5; //-0.1; //-1.0; //0.5; //-0.1; //0.25; //0.1;
  float qtof_negmean_num_stds = 0.5; //-0.1; //0.5;

  float min_singly_fval = -2;
  if(maldi_)
    min_singly_fval = -6;

   
  GammaDistribution* negGammaDistr = NULL ; 
  if (gamma_)
  {
	  negGammaDistr = (GammaDistribution*) negdistr_ ; 
	  //Xiuxia, negdist_ is already a GammaDistribution by initialization, why do type cast here?
  }

  if(gamma_ && charge_ == 0) 
  {
    negGammaDistr->setDistrMinval(min_singly_fval);
  }

  int numVals = getNumVals() ; 
  for(int k = 0; k < numVals ; k++) 
  {
	  float value = (*floatvals_)[k] ; 
	  if(nttdistr->isValue(k, 0) && (! gamma_ || negGammaDistr->aboveMin(value))) //Xiuxia, only consider the ntt=0 peptides. These are taken as incorrectly identified peptides.
	  {
		mean += value ;
		totsq += value * value ;
		tot++;
	  }
  }

  // what to do if not enough values above the min value of the negative distribution. 
  if(tot < MIN_NUM_PSEUDOS_) 
  { 
	  if(! gamma_) 
	  {
		  zero = 0.0; // don't need this
	  }
	  else 
	  {
		  if(charge_ == 0) 
		  {
			  zero = min_dataval_ - 3.0;
			  if(zero < min_singly_fval) 
			  {
				  zero = min_singly_fval;
			  }
		  }
		  else 
		  {
			  zero = min_dataval_ - 0.1;
		  }
		  if(zero < minval_) 
		  {    
			  zero = minval_;
		  }
	  } // gamma

	  tot = 0; // restart
	  mean = 0.0;
	  totsq = 0.0;

	  for(int k = 0; k < numVals; k++) 
	  {
		  float value = (*floatvals_)[k] ; 
		  if(value > zero && value < posmean[charge_] - posstdev[charge_]) 
		  {
			  mean += value - zero ;
			  totsq += (value - zero) * (value - zero);
			  tot++;
		  }
	  } // next

	  if(tot > 0) 
	  {
		  mean /= tot ; 
		  // calculate standar deviation
		  if(gamma_) 
		  {
			  stdev = totsq / tot;
		  }
		  else
		  {
			  stdev = totsq / tot - mean * mean;
		  }
	  }

	  float* newsettings = new float[3];
	  newsettings[0] = mean;
	  newsettings[1] = stdev;

	  if(!gamma_) 
	  {
		  newsettings[1] = sqrt(newsettings[1]);
	  }

	  newsettings[2] = zero;
	  negdistr_->init(newsettings);
	  delete [] newsettings;

	  if(! gamma_) 
	  {
		  negmean_ = mean - negmean_num_stds[charge_] * stdev;
	  }
	  else 
	  {
		  negmean_ = zero + mean - negmean_num_stds[charge_] * sqrt(stdev);
	  }
	  if(negmean_ > MAX_SINGLY_NEGMEAN) 
	  {
		  negmean_ = MAX_SINGLY_NEGMEAN;
	  }

	  //std::cerr << "setting negmean to " << negmean_ << " for charge " << (charge_+1) << " with zero " << zero << std::endl;

	  USE_TR_NEG_DISTR_ = True;
	  return False; // done
  } // if not enough pseudos
    
  mean /= tot;
  stdev = (totsq / tot) - mean * mean;
  stdev = sqrt(stdev);  // delete this later ?

  if(gamma_) 
  {
	  zero = mean - NUM_DEVS_ * stdev;
	  if(min_dataval_ - ZERO_SET_ * stdev > zero) 
	  {
		  zero = min_dataval_ - ZERO_SET_ * stdev;
	  }
  }
  else 
  {
    zero = 0.0;
  }


  
  if(maldi_) 
  {
	  negmean_ = mean - maldi_negmean_num_stds * stdev;
	  //std::cerr << "mean: " << mean << ", negmeanstdevs: " << negmean_num_stds[0] << ", stdev: " << stdev << std::endl;
	  //negmean_ = -1.4;
  }
  else if(qtof_) 
  {
	  negmean_ = mean - qtof_negmean_num_stds * stdev;
  }
  else
  {
	  negmean_ = mean - negmean_num_stds[charge_] * stdev;
  }


  // now recompute real mean and stdev
  mean = 0.0;
  totsq = 0.0;
  tot = 0;
    
  for(int k = 0; k < numVals; k++) 
  {
	  float value = (*floatvals_)[k] ; 
	  if(nttdistr->isValue(k, 0) && value > zero && (! gamma_ || negGammaDistr->aboveMin(value))) 
	  {
		  mean += value - zero;
		  totsq += (value - zero) * (value - zero);
		  tot++;
	  }
  } // next
  
  if(tot > 0)
  {
	  mean /= tot;
	  if(gamma_) 
	  {
		  stdev = totsq / tot;
	  }
	  else 
	  {
		  stdev = (totsq / tot) - mean * mean;
		  stdev = sqrt(stdev);
	  }
  }

  Boolean2 reset = False;
  float newposmean = posmean[charge_];
  float newposstdev = posstdev[charge_];

  while(! maldi_ && noDistr()) 
  {
    newposmean += 0.5;
    newposstdev += 0.95;
    setPositiveDistr(newposmean, newposstdev);
    //reset = True;
  }
  //if(reset) {
  //setPositiveDistr(newposmean, newposstdev);
  //}

  setNegativeDistr(mean, stdev, zero);
  return True;
}

float DiscrimValMixtureDistr::getnegmean() {return negmean_ ;}