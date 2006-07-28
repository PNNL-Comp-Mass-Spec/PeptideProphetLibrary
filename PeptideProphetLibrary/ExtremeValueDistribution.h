#ifndef EXTREME_H
#define EXTREME_H

#include <math.h>

#include "ContinuousDistribution.h"

/*

Program       : ExtremeValueDistribution for PeptideProphet                                                       
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

class ExtremeValueDistribution : public ContinuousDistribution {

 public:
  ExtremeValueDistribution(float maxdiff);
  void init(float* prior);
  void initUpdate(float* prior);
  void addVal(float wt, float val);
  void addVal(float wt, int val);
  Boolean2 update();
  float getProb(int val);
  float getExtremeValueProb(float val, float beta, float mu);
  float getProb(float val);
  void computeMoments(float mean, float std);
  float cumulative(float x, float mu, float beta);
  float extremeValueSlice(float num, float left_val, float right_val, float mu, float beta);
  float slice(float num, float left_val, float right_val);
  void printDistr();
  void writeDistr(FILE* fout);

 protected:
  float beta_;
  float mu_;

};



#endif
