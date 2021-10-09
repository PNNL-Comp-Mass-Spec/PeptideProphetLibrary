#include <string>
using namespace std;

#include "EnzymeDigestion.h"
#include <sstream>

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

EnzymeDigestion::EnzymeDigestion(char* sites, char* term_not_following, char* term_not_preceding,
                 int min_edge_dist, int min_dist) {
  recognition_sites_ = strCopy(sites);
  term_not_following_ = strCopy(term_not_following);
  term_not_preceding_ = strCopy(term_not_preceding);
  min_edge_dist_ = min_edge_dist;
  min_dist_ = min_dist;
}

EnzymeDigestion::~EnzymeDigestion() {
  if(recognition_sites_ != NULL)
    delete [] recognition_sites_;
  if(term_not_following_ != NULL)
    delete [] term_not_following_;
  if(term_not_preceding_ != NULL)
    delete [] term_not_preceding_;
}



// pep seq must be stripped
int EnzymeDigestion::numMissedCleavages(std::string peptide)
{
  const char *pep = peptide.c_str() ;
  if(strlen(pep) > 4 && pep[1] == '.' && pep[strlen(pep)-2] != '.')
  {
    //std::cerr << "cannot parse peptide " << pep << std::endl;
      std::stringstream str;
      str << "cannot parse peptide " << pep << std::endl;
      throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
    //exit(1);
  }
  int counter = 0;
  int nmc = 0;
  Boolean2 found = False;

  for(int k = 0; k < strlen(pep); k++) {
    if(counter >= min_dist_ && k >= min_edge_dist_ && k < strlen(pep) - min_edge_dist_) { // position ok
      if((k > 0 || termCanFollow(pep[k-1])) && (k < strlen(pep) - 1 || termCanPrecede(pep[k+1])) &&
     isCompatibleTerminus(pep[k])) {
    nmc++;
    counter = 0;
      }
    } // ok posn

  } // next pep posn
  return nmc;

}

// pep of form: P.XXXXX.F
int EnzymeDigestion::numCompatibleTermini(std::string peptide)
{
    const char *pep1 = peptide.c_str() ;

    //Xiuxia, July 06, 2006
    if(pep1[0] == '.')
    {
        peptide = '-' + peptide ;
    }
    if(pep1[strlen(pep1)-1] == '.')
    {
        peptide = peptide + '-' ;
    }

    const char *pep = peptide.c_str() ;

    if(strlen(pep) < 4 || pep[1] != '.' || pep[strlen(pep)-2] != '.')
    {
        //std::cerr << strlen(pep) << "\t" << pep[1] << "\t" << pep[strlen(pep)-2] << std::endl << std::endl ;
        //std::cerr << "cannot parse peptide " << pep << std::endl;
        std::stringstream str;
        str << strlen(pep) << "\t" << pep[1] << "\t" << pep[strlen(pep)-2] << std::endl << std::endl ;
        str << "cannot parse peptide " << pep << std::endl;
        throw gcnew System::Exception(gcnew System::String(str.str().c_str()));
        //exit(1);
    }
    int nct = 0;

    if(pep[0] == '-' || pep[0] == '1')
    {
        nct++;
    }
    else if(strlen(pep) > 2 && termCanPrecede(pep[2]))
    {
        for(int s = 0; s < strlen(recognition_sites_); s++)
        {
            if(pep[0] == recognition_sites_[s])
            {
                nct++;
                s = strlen(recognition_sites_); // done
            }
        }
    }
    if(pep[strlen(pep)-3] == '1' || pep[strlen(pep)-1] == '-')
    {
        nct++;
    }
    else if(strlen(pep) > 3 && termCanFollow(pep[strlen(pep)-4]))
    {
        for(int s = 0; s < strlen(recognition_sites_); s++)
        {
            if(pep[strlen(pep)-3] == recognition_sites_[s])
            {
                nct++;
                s = strlen(recognition_sites_); // done
            }
        }
    }
    return nct;
}

Boolean2 EnzymeDigestion::isCompatibleTerminus(char c) {
  for(int k = 0; k < strlen(recognition_sites_); k++)
    if(recognition_sites_[k] == c)
      return True;
  return False;
}

Boolean2 EnzymeDigestion::termCanPrecede(char c) {
  for(int k = 0; k < strlen(term_not_preceding_); k++)
    if(term_not_preceding_[k] == c)
      return False;
  return True;
}

Boolean2 EnzymeDigestion::termCanFollow(char c) {
  for(int k = 0; k < strlen(term_not_following_); k++)
    if(term_not_following_[k] == c)
      return False;
  return True;
}

char* EnzymeDigestion::strCopy(char* orig) {
  char* output = new char[strlen(orig)+1];
  strcpy(output, orig);
  output[strlen(orig)] = 0;
  return output;
}
