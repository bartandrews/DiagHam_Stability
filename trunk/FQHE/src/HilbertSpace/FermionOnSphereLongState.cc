////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of state of fermion on a sphere without any restriction       //
//    on the number of fermion nor the number of states that can be reached   //
//                                                                            //
//                        last modification : 02/03/2004                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////



#include "config.h"
#include "HilbertSpace/FermionOnSphereLongState.h"


// default constructor
// 

FermionOnSphereLongState::FermionOnSphereLongState()
{
  this->StateDescription = 0;
}
  
// copy constructor
// 
// state = reference on the state to copy
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1

FermionOnSphereLongState::FermionOnSphereLongState(FermionOnSphereLongState& state, const int& reducedNbrState)
{
  this->StateDescription = new unsigned long [reducedNbrState + 1];
  for (int i = 0; i <= reducedNbrState; ++i)
    this->StateDescription[i] = state.StateDescription[i];
}
  
  
// basic constructor
// 
// reducedNbrState = reduced number of state (aka the number of unsigned long per state) minus 1

FermionOnSphereLongState::FermionOnSphereLongState(const int& reducedNbrState)
{
  this->StateDescription = new unsigned long [reducedNbrState + 1];
  for (int i = 0; i <= reducedNbrState; ++i)
    this->StateDescription[i] = (unsigned long) 0;
}
  
  
// destructor
// 

FermionOnSphereLongState::~FermionOnSphereLongState()
{
  if (this->StateDescription != 0)
    delete[] this->StateDescription;
}
  
