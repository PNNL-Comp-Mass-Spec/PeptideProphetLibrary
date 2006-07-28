#include "MassDifferenceDiscrMixtureDistr.h"

/*

Program       : MassDifferenceDiscrMixtureDistr for PeptideProphet                                                       
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

MassDifferenceDiscrMixtureDistr::MassDifferenceDiscrMixtureDistr(int charge, char* name, char* tag, float range, float window) : DiscreteMixtureDistr(charge, 2 * ((int)(range/window)) + 1, name, tag) {

  bindefs_ = new Array<char*>;
  range_ = range;
  window_ = window;

  float bin;
  int num = 0;
  char* next;
  for(int k = 0; k < numbins_; k++) {
    bin = -1.0 * range_ + (k * window_);
    if(bin < 0.0)
      num++;
    next = new char[num + 10];
    sprintf(next, "massd=%0.1f", bin, window_/2.0);
    next[num+9] = 0;
    bindefs_->insertAtEnd(next);
  }

  maxdiff_ = 0.001;
  negOnly_ = True;
  DiscreteMixtureDistr::init(NULL);

}

int MassDifferenceDiscrMixtureDistr::inttranslate(char* val) {
  float value = atof(val);
  for(int k = 0; k < numbins_; k++)
    if(value <= (-1.0 * range_ + (k * window_)) + window_ / 2.0)
      return k;
  return numbins_ - 1;
}

Boolean2 MassDifferenceDiscrMixtureDistr::haveDataWithValue(int bin) {
  for(int k = 0; k < intvals_->length(); k++)
    if((*intvals_)[k] == bin)
      return True;
  return False;
}

void MassDifferenceDiscrMixtureDistr::writeDistr(FILE* fout) {
  
  fprintf(fout, "%s\n", name_);
  fprintf(fout, "\tpos: ");
  fprintf(fout, "(");
  int next;
  int counter = 0;
  int column_width = 4;
  for(int k = 0; k < numbins_; k++) {
    if(haveDataWithValue(k)) {
      counter++;
      fprintf(fout, "%0.2f %s", posdistr_->getProb(k), (*bindefs_)[k]);
      if(k < numbins_ - 1 && haveDataWithValue(k+1)) {
	fprintf(fout, ", ");
      }
      if(counter%column_width == 0)
	fprintf(fout, "\n\t\t");
    } // if have value
  }
  fprintf(fout, ")\n");
  
  fprintf(fout, "\tneg: ");
  fprintf(fout, "(");
  counter = 0;
  for(int k = 0; k < numbins_; k++) {
    if(haveDataWithValue(k)) {
      counter++;
      fprintf(fout, "%0.2f %s", negdistr_->getProb(k), (*bindefs_)[k]);
      if(k < numbins_ - 1 && haveDataWithValue(k+1)) {
	fprintf(fout, ", ");
      }
      if(counter%column_width == 0)
	fprintf(fout, "\n\t\t");
    } // if have value
  }
  
  fprintf(fout, ")\n");
}
