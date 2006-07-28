#include "ScoreCalculator.h"
#include <sstream>

/*

Program       : ScoreCalculator for discr_calc of PeptideProphet 
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

ScoreCalculator::ScoreCalculator() { }

ScoreCalculator::ScoreCalculator(char* enz) {
  // in the future, want to specify the enzyme specificity through the c'tor	

  // for now
  if(enz == NULL)
    enzyme_ = NULL;
  else {
    enzyme_ = new char[strlen(enz)+1];
    strcpy(enzyme_, enz);
    enzyme_[strlen(enz)] = 0;	
  }
  enzyme_digestion_ = (new EnzymeSpecificity())->getEnzymeDigestion(enzyme_);

  if(enzyme_digestion_ == NULL) {
    //std::cerr << "could not get enzyme digestion for " << enzyme_ << std::endl;
    //exit(1);   

	  std::stringstream str;
	  str << "could not get enzyme digestion for " << enzyme_ << std::endl;

	  throw new System::Exception(str.str().c_str());	  
  }

}

int ScoreCalculator::numMissedCleavages(std::string peptide) {
  return enzyme_digestion_->numMissedCleavages(peptide);
}

int ScoreCalculator::numCompatibleTermini(std::string peptide) {
  return enzyme_digestion_->numCompatibleTermini(peptide);
}


// take off all special symbols, as well as preceding and following aa
std::string ScoreCalculator::strip(std::string peptide, Boolean2 remove_mods) 
{
  int start = 0;
  int stop = peptide.length() ;

  const char *pep = peptide.c_str() ; 
  char* output = NULL;
  if(stop > 4 && pep[1] == '.')
	  start = 2;
  if(strlen(pep) > 4 && pep[stop-2] == '.')
	  stop = stop-3;

  if(! remove_mods) 
  {
    output = new char[stop - start + 2];
    strncpy(output, pep+start, stop - start + 1);
    output[stop - start + 1] = 0;
	std::string outStr = output ; 
	delete [] output ; 
    return outStr;
  }

  int length = 0;
  for(int k = start; k <= stop; k++) 
    if(pep[k] >= 'A' && pep[k] <= 'Z')
      length++;
  output = new char[length+1];
  length = 0;
  for(int k = start; k <= stop; k++) 
    if(pep[k] >= 'A' && pep[k] <= 'Z')
      output[length++] = pep[k];
  output[length] = 0;

  std::string outStr = output ; 
  delete [] output ; 
  return outStr;
}

char* ScoreCalculator::parseSpectrum(char* spec) {
  if(spec == NULL) {
    //std::cerr << "null spectrum" << std::endl;
    //exit(1);
	  std::stringstream str;
	  str << "null spectrum" << std::endl;
	  throw new System::Exception(str.str().c_str());
  }
  char* output = NULL;
  if(strlen(spec) > 2 && spec[strlen(spec)-2] == '.') {
    output = new char[strlen(spec)-1];
    strncpy(output, spec, strlen(spec)-2);
    output[strlen(spec)-2] = 0;
  }
  else {
    //std::cerr << "cannot parse spec " << spec << std::endl;
	  std::stringstream str;
	  str << "cannot parse spec " << spec << std::endl;
	  throw new System::Exception(str.str().c_str());
  }
  return output;
}

Boolean2 ScoreCalculator::isNumber(char c) {
  return c >= '0' && c <= '9';
}
