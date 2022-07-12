////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of MPO matrix for density operator                 //
//                                                                            //
//                        last modification : 29/07/2016                      //
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
#include "Tools/FQHEMPS/FQHEMPODensityOperator.h"
#include "Matrix/SparseRealMatrix.h"


// default constructor 
//

FQHEMPODensityOperator::FQHEMPODensityOperator()
{
}

// constructor 
//
// maxOccupation = maximum occupation for a single orbital
// prefactor = repfactor in front of the density operator

FQHEMPODensityOperator::FQHEMPODensityOperator(int maxOccupation, double prefactor)
{
  this->NbrOMatrices =  maxOccupation + 1;
  this->BondDimension = 2;
  this->UpperPhysicalIndices = new unsigned long[this->NbrOMatrices];
  this->LowerPhysicalIndices = new unsigned long[this->NbrOMatrices];
  this->RealOMatrices = new SparseRealMatrix[this->NbrOMatrices];
  int* TmpNbrElementPerRow = new int[2];
  TmpNbrElementPerRow[0] = 2;
  TmpNbrElementPerRow[1] = 1;
  for (int i = 0; i < this->NbrOMatrices; ++i)
    {
      this->UpperPhysicalIndices[i] = (0x1ul << i) - 0x1ul;
      this->LowerPhysicalIndices[i] = (0x1ul << i) - 0x1ul;
      this->RealOMatrices[i] = SparseRealMatrix(2, 2, TmpNbrElementPerRow);
      this->RealOMatrices[i].SetMatrixElement(0, 0, 1.0);
      this->RealOMatrices[i].SetMatrixElement(0, 1, prefactor * ((double) i));
      this->RealOMatrices[i].SetMatrixElement(1, 1, 1.0);
      this->RealOMatrices[i].LockMatrix();
    }
}

// destructor
//

FQHEMPODensityOperator::~FQHEMPODensityOperator()
{
} 

// get the boundary indices of the MPO representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index

void FQHEMPODensityOperator::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex)
{
  rowIndex = 0;
  columnIndex = 1;
}

// get the name describing the O matrices 
// 
// return value = name 

char* FQHEMPODensityOperator::GetName()
{
  char* TmpString = new char[16];
  sprintf(TmpString, "density");
  return TmpString;
}

// set a new prefactor in front of the operator
//
// prefactor = prefactor in front of the operator

void FQHEMPODensityOperator::SetPrefactor(double prefactor)
{
  for (int i = 0; i < this->NbrOMatrices; ++i)
    {
      this->RealOMatrices[i].SetMatrixElement(0, 1, prefactor * ((double) i));
    }
}

