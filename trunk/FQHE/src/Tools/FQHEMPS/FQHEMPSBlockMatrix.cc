////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of MPS matrix built as a block diagonal matrix         //
//                                                                            //
//                        last modification : 17/02/2013                      //
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
#include "Tools/FQHEMPS/FQHEMPSBlockMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"

#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor 
//

FQHEMPSBlockMatrix::FQHEMPSBlockMatrix()
{
}

// constructor from two MPS matrices (the number of B matrices has to be identical for all of them)
//
// matrix1 = first MPS matrix
// matrix2 = second MPS matrix

FQHEMPSBlockMatrix::FQHEMPSBlockMatrix(AbstractFQHEMPSMatrix* matrix1, AbstractFQHEMPSMatrix* matrix2)
{
  this->NbrBMatrices = matrix1->GetNbrMatrices();
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  for (int i = 0; i < this->NbrBMatrices; ++i)
    this->RealBMatrices[i] = CreateBlockDiagonalMatrix(matrix1->GetMatrices()[i], matrix2->GetMatrices()[i]);
}

// create the B matrices for the block state
//

void FQHEMPSBlockMatrix::CreateBMatrices ()
{
}

