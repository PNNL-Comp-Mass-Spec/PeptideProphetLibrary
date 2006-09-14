#include "SequestDiscrimScoreCalculator.h"
#include "GlobalVariable.h"

/*

Program       : DiscriminantCalculator for discr_calc of PeptideProphet 
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

SequestDiscrimScoreCalculator::SequestDiscrimScoreCalculator() : ScoreCalculator() { }


SequestDiscrimScoreCalculator::SequestDiscrimScoreCalculator(std::vector<SequestResult> &results, Boolean2 exclude_deltastars, Boolean2 windows, Boolean2 massd, Boolean2 modify_deltas, Boolean2 maldi, char* enz) : ScoreCalculator(enz) {
  modify_deltas_ = modify_deltas;
  maldi_ = maldi;
  maldi_set_ = True;
  init(results, exclude_deltastars, windows, massd);
}

SequestDiscrimScoreCalculator::SequestDiscrimScoreCalculator(std::vector<SequestResult> &results, Boolean2 exclude_deltastars, Boolean2 windows, Boolean2 massd, Boolean2 modify_deltas, Boolean2 maldi) : ScoreCalculator(NULL) {
  modify_deltas_ = modify_deltas;
  maldi_ = maldi;
  maldi_set_ = True;
  init(results, exclude_deltastars, windows, massd);
}

void SequestDiscrimScoreCalculator::init(std::vector<SequestResult> &results, Boolean2 exclude_deltastars, Boolean2 windows, Boolean2 massd) 
{
  abort = false;
  set_deltastars_zero_ = True;
  windows_ = windows;

  if(exclude_deltastars) {
    exclude_deltastars_ = True;
    set_deltastars_zero_ = False;
  }

  discr_funcs_ = new SequestDiscrimFunction* [numCharge];

  for(int ch = 0; ch < numCharge; ch++)		//Xiuxia, to consider charge = 4, 5
    if(massd)
      discr_funcs_[ch] = new AbbrevSequestDiscrimFunction(ch);
    else if(maldi_ && ch == 0)
      discr_funcs_[ch] = new ToftofAbbrevSequestDiscrimFunction(); 
    else 
      discr_funcs_[ch] = new SequestDiscrimFunction(ch);


  processData(results);
}

void SequestDiscrimScoreCalculator::processData(std::vector<SequestResult> &results)
{

  //std::cerr << " Computing SEQUEST discriminant scores";
  //if(modify_deltas_)
    //std::cerr << " (leaving original delta* values)";
  //else if(exclude_deltastars_)
  //  //std::cerr << " (excluding delta* entries)";
  //else if(set_deltastars_zero_)
  //  //std::cerr << " (setting all delta* to zero)";
  //if(enzyme_ != NULL)
    //std::cerr << " (" << enzyme_ << ")" << std::endl;
  //std::cerr << std::endl;
	

  float discr_score;
  int ntt;
  int nmc;
  std::string stripped ;
  int degen;

  int format = -1;

  int numResults = results.size() ; 
  for (int resultNum = 0 ; resultNum < numResults ; resultNum++)
  {
	  if(abort)
	  {
		  throw new System::Exception("Process aborted in SequestDiscrimScoreCalculator::processData");
	  }

	  if(results[resultNum].xcorr_ > 0.0) 
	  {
		  std::string currentPeptide = results[resultNum].peptide_ ; 
		  int charge = results[resultNum].charge_ - 1; 
		  // Xiuxia, 05/16/2006

		  //if (charge <= 2)
		  if (charge < numCharge) //Xiuxia
		  {
		    stripped = strip(currentPeptide, 1); // remove all modifications also
			nmc = numMissedCleavages(stripped);
			// numMissedCleavages() is a function in the base class ScoreCalculator

			ntt = numCompatibleTermini(currentPeptide);

			results[resultNum].mint_numMissedCleavages = nmc;
			results[resultNum].mint_numTT = ntt;
			
			results[resultNum].mdbl_discriminantScore = discr_funcs_[charge]->getDiscriminantScore(results[resultNum]);	
		  }
      } // if process
    } // if real data line
}

