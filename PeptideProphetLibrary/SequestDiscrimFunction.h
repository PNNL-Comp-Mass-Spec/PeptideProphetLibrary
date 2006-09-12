#ifndef SEQDISCR_FUN_H
#define SEQDISCR_FUN_H

#include <string.h>
#include <iostream>
#include <fstream>

#include "SequestResult.h"
#include "DiscriminantFunction.h"
#include "GlobalVariable.h"

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


class SequestDiscrimFunction : public DiscriminantFunction {

 public:

  SequestDiscrimFunction(int charge);
  Boolean2 isComputable(SequestResult &result);
  float getDiscriminantScore(SequestResult &result);
  float getXcorrP(float xcorr, int peplen);
  int getPepLen(std::string peptide) ;

  void writeCoef(char *d_output_file_param) ;

 protected:

	 float consts[numCharge] ;
	 float xcorrs[numCharge] ;
	 float deltas[numCharge] ;
	 float ranks[numCharge] ;
	 float massdiffs[numCharge] ;
	 float max_pep_lens[numCharge] ;
	 float num_frags[numCharge] ;

	//float consts[] = {0.646, -0.959, -1.460, -0.959, -0.959};	//Xiuxia, added the last two numbers for charge = 4, 5
	//float xcorrs[] = {5.49, 8.362, 9.933, 8.362, 8.362};//Xiuxia, added the last two numbers for charge = 4, 5
	//float deltas[] = {4.643, 7.386, 11.149, 7.386, 7.386};//Xiuxia, added the last two numbers for charge = 4, 5
	//float ranks[] = {-0.455, -0.194, -0.201, -0.194, -0.194};//Xiuxia, added the last two numbers for charge = 4, 5
	//float massdiffs[] =  {-0.84, -0.314, -0.277, -0.314, -0.314};//Xiuxia, added the last two numbers for charge = 4, 5
	//int max_pep_lens[] = {100, 15, 25, 50, 50}; //Xiuxia, ?
	//int num_frags[] = {2, 2, 4, 6, 6};	//Xiuxia, ?

	float xcorr_p_wt_;
	float delta_wt_;
	float log_rank_wt_;
	float abs_massd_wt_;
	int max_pep_len_;
	int num_frags_;

}; // class

#endif
