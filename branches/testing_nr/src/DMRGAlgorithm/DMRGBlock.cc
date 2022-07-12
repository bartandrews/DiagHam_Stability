////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of elementary DMRG block                      //
//                                                                            //
//                        last modification : 31/05/2001                      //
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


#include "DMRGAlgorithm/DMRGBlock.h"


// constructor from datas
//
// hamiltonian = Hamiltonian associated to the block
// iterationID = identfifcator of iteration associated to the block

DMRGBlock::DMRGBlock (ExplicitHamiltonian* hamiltonian, int iterationID) 
{
  this->Hamiltonian = hamiltonian;
  this->IterationID = iterationID;
}

// destructor
//

DMRGBlock::~DMRGBlock () 
{
}
