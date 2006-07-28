#ifndef SEARCHRESULT_H
#define SEARCHRESULT_H

#include <string.h>
#include <iostream>

#include "sysdepend.h"
#include "constants.h"



/*

Program       : SearchResult for discr_calc of PeptideProphet 
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


class SearchResult {

 public:

  SearchResult();

  virtual Boolean2 isMaldi(char* spec);
  char* strCopy(char* orig);
  Boolean2 isProcessed();

  int charge_;
  std::string spectrum_;
  std::string protein_;
  std::string peptide_;
  Boolean2 maldi_; // whether or not maldi data
  Boolean2 degen_; // whether or not degen
  virtual char* getName() = 0;

  virtual std::ostream& print(std::ostream& os);


 protected:

  virtual void init();
  Boolean2 processed_; // whether or not ok


};

#endif
