#include "VariableOffsetMassDiffDiscrMixtureDistr.h"

/*

Program       : VariableOffsetMassDiffDiscrMixtureDistr for PeptideProphet
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


VariableOffsetMassDiffDiscrMixtureDistr::VariableOffsetMassDiffDiscrMixtureDistr(int charge, char* name, char* tag, float range, float window, float orig) : MassDifferenceDiscrMixtureDistr(charge, name, tag, range, window) {
  offset_init_ = orig;
  offset_ = offset_init_;
  vals_ = new Array<float>;
  offset_set_ = False;
  update_ctr_ = 0;
  min_ctr_ = 2; // min value for offset update
}

void VariableOffsetMassDiffDiscrMixtureDistr::enter(int index, char* val) {
  assert(intvals_ != NULL);
  intvals_->insertAtEnd(inttranslate(val));
  vals_->insertAtEnd(atof(val)); // put the float value here
}

int VariableOffsetMassDiffDiscrMixtureDistr::getIntegralValue(float val) {
  for(int k = 0; k < numbins_; k++)
    if(val <= (-1.0 * range_ + (k * window_) + offset_) + window_ / 2.0)
      return k;
  return numbins_ - 1;

}

Boolean2 VariableOffsetMassDiffDiscrMixtureDistr::update(Array<float>* probs) {
  Boolean2 output = False;
  if(! offset_set_ && update_ctr_ >= min_ctr_) {
    //std::cerr << "here1" << std::endl;
    float new_offset = getMode(0.1, probs);
    if(new_offset - offset_ > maxdiff_ || offset_ - new_offset > maxdiff_) { // update
      output = True;
      offset_ = new_offset;
      assert(vals_->length() == intvals_->length());
      for(int k = 0; k < vals_->length(); k++)
        intvals_->replace(k, getIntegralValue((*vals_)[k]));
    }
    else {
      offset_set_ = True; // done
    }
  } // if update offset
  //std::cerr << "here2" << std::endl;
  if(! offset_set_)
    update_ctr_++;
  if(DiscreteMixtureDistr::update(probs) || output)
    return True;

  return False;

}

float VariableOffsetMassDiffDiscrMixtureDistr::getMode(float window, Array<float>* probs) {
  float min_tot = 5.0;
  int num_windows = (int)(2 * range_ / window) + 1;
  if(probs == NULL)
    return offset_init_;

  float* win = new float[num_windows];
  for(int k = 0; k < num_windows; k++)
    win[k] = 0.0;
  float max = 0.0;
  int max_ind = -1;
  float tot = 0.0;

  assert(probs->length() == vals_->length());


  for(int k = 0; k < probs->length(); k++) {
    int next = (int)(((*vals_)[k] + range_) / window);
    //std::cerr << k << " out of " << probs->length() << ": " << next << " vs " << num_windows << std::endl;
    if(next < 0)
      next = 0;
    if(next >= num_windows)
      next = num_windows - 1;
    win[next] += (*probs)[k];
    tot += (*probs)[k];
  }

  // now find the max
  if(tot < min_tot) {
    delete[] win;
    return offset_init_;
  }

  for(int k = 0; k < num_windows; k++)
    if(win[k] > max) {
      max = win[k];
      max_ind = k;
    }
  delete[] win;

  if(max_ind == -1)
    return offset_init_;
  //std::cerr << update_ctr_ << ": " << max_ind << " " << max << std::endl;
  return ((float)max_ind * window - range_);
}

void VariableOffsetMassDiffDiscrMixtureDistr::writeDistr(FILE* fout) {

  fprintf(fout, "%s (offset: %0.2f)\n", name_, offset_);
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
