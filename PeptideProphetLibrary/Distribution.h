#ifndef DISTR_H
#define DISTR_H


#include <iostream>
#include <stdio.h>
#include <math.h>

#include "sysdepend.h"

/*

Program       : Distribution for PeptideProphet                                                       
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

class Distribution {

public:
  Distribution();
  Distribution(float maxdiff);
  ~Distribution(); 

  virtual float getProb(float val) = 0;
  virtual float getProb(int val) = 0;
  virtual void init(float* priors) = 0;
  virtual void addVal(float wt, float val) = 0;
  virtual void addVal(float wt, int val) = 0;
  virtual Boolean2 update() = 0;
  virtual void initUpdate(float* prior) = 0;
  virtual void printDistr() = 0;
  virtual void writeDistr(FILE* fout) = 0;
  virtual void setMinVal(float min); 

  protected: 

float* tot_; 
float totwt_;
float* newtot_; // for update 
float newtotwt_; // for update Boolean2 set_;
float maxdiff_;
 Boolean2 set_;
 float minval_;
};

















#endif
