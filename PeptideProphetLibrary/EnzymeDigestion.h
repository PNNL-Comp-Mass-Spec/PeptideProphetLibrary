#ifndef ENZ_DIG_H
#define ENZ_DIG_H

#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "sysdepend.h"

/*

Program       : EnzymeDigestion for PeptideProphet                                                       
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


class EnzymeDigestion {

 public:

  EnzymeDigestion(char* sites, char* term_not_following, char* term_not_preceding, 
          int min_edge_dist, int min_dist);
  ~EnzymeDigestion();
  virtual int numMissedCleavages(std::string pep);
  virtual int numCompatibleTermini(std::string peptide);
  Boolean2 isCompatibleTerminus(char c);
  Boolean2 termCanPrecede(char c);
  Boolean2 termCanFollow(char c);
  char* strCopy(char* orig);

 protected:

  char* recognition_sites_;
  char* term_not_following_;
  char* term_not_preceding_;
  int min_edge_dist_;
  int min_dist_;
};


#endif
