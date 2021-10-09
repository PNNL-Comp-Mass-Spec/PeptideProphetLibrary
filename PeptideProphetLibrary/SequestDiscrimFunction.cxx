#include "SequestDiscrimFunction.h"
#include <math.h>
#include <sstream>

using namespace std ;
using namespace System;


/*

Program       : DiscriminantFunction for discr_calc of PeptideProphet
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


SequestDiscrimFunction::SequestDiscrimFunction(int charge) : DiscriminantFunction(charge) {

  //float consts[] = {0.646, -0.959, -1.460, -0.959, -0.959};	//Xiuxia, added the last two numbers for charge = 4, 5
  //float xcorrs[] = {5.49, 8.362, 9.933, 8.362, 8.362};//Xiuxia, added the last two numbers for charge = 4, 5
  //float deltas[] = {4.643, 7.386, 11.149, 7.386, 7.386};//Xiuxia, added the last two numbers for charge = 4, 5
  //float ranks[] = {-0.455, -0.194, -0.201, -0.194, -0.194};//Xiuxia, added the last two numbers for charge = 4, 5
  //float massdiffs[] =  {-0.84, -0.314, -0.277, -0.314, -0.314};//Xiuxia, added the last two numbers for charge = 4, 5
  //int max_pep_lens[] = {100, 15, 25, 50, 50}; //Xiuxia, ?
  //int num_frags[] = {2, 2, 4, 6, 6};	//Xiuxia, ?

	consts[0] = 0.646 ;
	consts[1] = -0.959 ;
	consts[2] = -1.460 ;
	consts[3] = -0.959 ;
	consts[4] = -0.959 ;

	xcorrs[0] = 5.49 ;
	xcorrs[1] = 8.362 ;
	xcorrs[2] = 9.933 ;
	xcorrs[3] = 8.362 ;
	xcorrs[4] = 8.362 ;

	deltas[0] = 4.643 ;
	deltas[1] = 7.386 ;
	deltas[2] = 11.149 ;
	deltas[3] = 7.386 ;
	deltas[4] = 7.386 ;

	ranks[0] = -0.455 ;
	ranks[1] = -0.194 ;
	ranks[2] = -0.201 ;
	ranks[3] = -0.194 ;
	ranks[4] = -0.194 ;

	massdiffs[0] = -0.84 ;
	massdiffs[1] = -0.314 ;
	massdiffs[2] = -0.277 ;
	massdiffs[3] = -0.314 ;
	massdiffs[4] = -0.314 ;

	max_pep_lens[0] = 100 ;
	max_pep_lens[1] = 15 ;
	max_pep_lens[2] = 25 ;
	max_pep_lens[3] = 50 ;
	max_pep_lens[4] = 50 ;

	num_frags[0] = 2 ;
	num_frags[1] = 2 ;
	num_frags[2] = 4 ;
	num_frags[3] = 6 ;
	num_frags[4] = 6 ;

	const_ = consts[charge_];
	xcorr_p_wt_ = xcorrs[charge_];
	delta_wt_ = deltas[charge_];
	log_rank_wt_ = ranks[charge_];
	abs_massd_wt_ = massdiffs[charge_];
	max_pep_len_ = max_pep_lens[charge_];
	num_frags_ = num_frags[charge_];
}



Boolean2 SequestDiscrimFunction::isComputable(SequestResult &result) {
  return result.xcorr_ > 0.0 ;
}

float SequestDiscrimFunction::getDiscriminantScore(SequestResult &seqresult)
{
  float tot = const_;
  tot += xcorr_p_wt_ * getXcorrP(seqresult.xcorr_, getPepLen(seqresult.peptide_));
  tot += delta_wt_ * seqresult.delta_;
  float lg_val = log((float)(1.0*seqresult.rank_)) ;
  tot += log_rank_wt_ * lg_val;
  tot += abs_massd_wt_ * abs(seqresult.massdiff_);
  Boolean2 writeMultivariate = False;

  if(writeMultivariate) {
	  std::ofstream fMulti("multivariate.txt", std::ios::app);
    if(! fMulti) {
      //std::cerr << "cannot append multivariate info for " << seqresult.spectrum_.c_str() << std::endl;
      //exit(1);
	  std::stringstream str;
	  str << "cannot append multivariate info for " << seqresult.spectrum_.c_str() << std::endl;
	  throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
    }
	double lg_val = log(1.0*seqresult.rank_) ;
    fMulti << seqresult.spectrum_.c_str() << "\t" << getXcorrP(seqresult.xcorr_, getPepLen(seqresult.peptide_)) << "\t" << seqresult.delta_ << "\t" << lg_val<< "\t" << abs(seqresult.massdiff_) << "\t" << tot << std::endl;
    fMulti.close();
  }

   return tot;
}


int SequestDiscrimFunction::getPepLen(std::string peptide)
{
	const char *pep = peptide.c_str() ;

	if (peptide.length() == 0)
		return 0 ;
	int start = 0;
	int end = peptide.length() ;

	if(end > 4 && pep[1] == '.' && pep[end-2] == '.')
	{
		start = 2;
		end = strlen(pep) - 2;
	}
	int tot = 0;
	for(int k = start; k < end; k++)
		if(pep[k] >= 'A' && pep[k] <= 'Z')
			tot++;
	return tot;
}

float SequestDiscrimFunction::getXcorrP(float xcorr, int peplen) {
  int eff_pep_len = peplen;
  if(eff_pep_len > max_pep_len_)
    eff_pep_len = max_pep_len_;
  float lg_xcorr = log(xcorr) ;
  float lg_eff_len = log((float)(1.0*eff_pep_len * num_frags_)) ;
  return lg_xcorr / lg_eff_len;
}

void SequestDiscrimFunction::writeCoef(char *d_output_file_param) {

  std::ofstream fout(d_output_file_param);
  if(! fout) {
    //std::cerr << "could not open filename" << std::endl;
	  std::stringstream str;
	  str << "could not open filename" << std::endl;
	  throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
    //exit(1);
  }

  fout << "CS" << "\t" << "Xcorr" << "\t" << "DeltaCn2" << "\t" << "RankSp" << "\t" << "DelM" << "\t" << "Const" << "\n" ;

  for (int charge =0; charge < numCharge; charge++) {
	  fout << charge+1 << "\t" << xcorrs[charge] << "\t" << deltas[charge] << "\t" << ranks[charge] << "\t" << massdiffs[charge] << "\t" << consts[charge] << "\n" ;
  }

  fout.close() ;
}
