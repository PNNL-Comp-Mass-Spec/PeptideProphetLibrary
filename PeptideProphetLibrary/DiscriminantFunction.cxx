#include "DiscriminantFunction.h"
#include "GlobalVariable.h"

/*

Program       : DiscriminantFunction for discr_calc of PeptideProphet
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


DiscriminantFunction::DiscriminantFunction(int charge) {

  charge_ = charge;
  error(charge_);
}


void DiscriminantFunction::error(int charge) {
  //if(charge < 0 || charge > 2) {
    if(charge < 0 || charge > numCharge-1) {    //Xiuxia
    //std::cerr << "illegal charge: " << charge << std::endl;
        throw gcnew System::Exception(gcnew System::String("illegal charge: " + charge));
    exit(1);
  }
}
