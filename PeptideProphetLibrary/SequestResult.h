#ifndef SEQRESULT_H
#define SEQRESULT_H

#include<string.h>
#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>


#include "SearchResult.h"

#define SIZE_BUF 8192
#define VAL_UNCERTAINTY 1

/*

Program       : SequestResult for discr_calc of PeptideProphet 
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


class SequestResult : public SearchResult {

 public:

  SequestResult();
  ~SequestResult();
  char* getName();
  int dataset_num_ ;// Xiuxia, 06/05/2006
  int ScanNumber;
  int charge_;
  float massdiff_;
  float xcorr_;
  float delta_;
  int mint_numMissedCleavages ; //DJ 5/14/2006
  int mint_numTT ; // DJ 5/14/2006
  double mdbl_discriminantScore ; 
  Boolean2 deltastar_;
  int rank_;
  int format_; // input format (index)

 protected:

};

#endif
