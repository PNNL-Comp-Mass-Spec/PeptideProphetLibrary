#include "AspNEnzymeDigestion.h"
#include<sstream>

/*

Program       : AspNEnzymeDigestion for PeptideProphet                                                       
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

AspNEnzymeDigestion::AspNEnzymeDigestion() : EnzymeDigestion("D", "", "", 1, 0) { }

int AspNEnzymeDigestion::numCompatibleTermini(std::string peptide) 
{
	const char *pep = peptide.c_str() ; 

	if(strlen(pep) < 4 || pep[1] != '.' || pep[strlen(pep)-2] != '.') 
	{
		//std::cerr << "cannot parse peptide " << pep << std::endl;
		std::stringstream str;
		str << "cannot parse peptide " << pep << std::endl;
		throw new System::Exception(str.str().c_str());		
		//exit(1);
	}

	int nct = 0;
	if(pep[0] == '-' || pep[0] == '1')
	{
		nct++;
	}
	else if(strlen(pep) > 2)
	{
		for(int s = 0; s < strlen(recognition_sites_); s++) 
		{
			if(pep[2] == recognition_sites_[s]) 
			{
				nct++;
				s = strlen(recognition_sites_); // done
			}
			if(pep[strlen(pep)-3] == '1' || pep[strlen(pep)-1] == '-')
			{
				nct++;
			}
			else 
			{
				for(int s = 0; s < strlen(recognition_sites_); s++) 
				{
					if(pep[strlen(pep)-1] == recognition_sites_[s]) 
					{
						nct++;
						s = strlen(recognition_sites_); // done
					}
				}
			}
		}
	}
	return nct;
}



