#ifndef CONT_DISTR_H
#define CONT_DISTR_H


#include "Distribution.h"

/*

Program       : ContinuousDistribution for PeptideProphet                                                       
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


class ContinuousDistribution : public Distribution {

  public:
  ContinuousDistribution();
  ContinuousDistribution(float maxdiff);
  ~ContinuousDistribution();

  virtual void init(float* prior);
  virtual float slice(float left_val, float right_val);
  virtual float slice(float num, float left_val, float right_val) = 0;
  virtual float getMean();
  virtual float getStdev();
  void resetTot();
  void setStdev(float stdev);
  void ignoreStdev();


  float mean_;//Xiuxia moved to public
  float stdev_;//Xiuxia moved to public

  protected:

  //float mean_;
  //float stdev_;
  float newtotsq_;
  Boolean2 use_stdev_;

};






# endif
