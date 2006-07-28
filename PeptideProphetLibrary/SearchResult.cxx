#include "SearchResult.h"

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


/*
format -1 uncertain
format 0 no mw column
format 1 has mw column
*/

SearchResult::SearchResult() { }


void SearchResult::init() {
  processed_ = False;
}

Boolean2 SearchResult::isProcessed() { return processed_; }



char* SearchResult::strCopy(char* orig) {
  char* output = new char[strlen(orig)+1];
  strcpy(output, orig);
  output[strlen(orig)] = 0;
  return output;
}


Boolean2 SearchResult::isMaldi(char* spec) {
  return False;
  // need to find a way in future....
}

std::ostream& SearchResult::print(std::ostream& os) 
{
  os << charge_ << " " << protein_.c_str() << " " << peptide_.c_str() << " " << spectrum_.c_str() << std::endl;
  return os;
}

