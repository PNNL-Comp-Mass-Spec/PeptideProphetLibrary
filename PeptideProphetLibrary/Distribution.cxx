
#include "Distribution.h"

/*

Program       : Distribution for PeptideProphet                                                       
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 

Distribution

This software will be available as open source in the future.  However, at this
time we are sharing it with you as a collaborator and ask that you do not
redistribute it.

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


Distribution::Distribution() 
{
	newtot_ = NULL ; 
}

Distribution::Distribution(float maxdiff) {
  maxdiff_ = maxdiff;
  set_ = False;
  newtot_ = NULL;
}


Distribution::~Distribution() { 
  if(tot_ != NULL) {
    delete [] tot_;
  }
  if(newtot_ != NULL) {
    delete [] newtot_;
  }

}


void Distribution::setMinVal(float min) { 
  minval_ = min; 
}

