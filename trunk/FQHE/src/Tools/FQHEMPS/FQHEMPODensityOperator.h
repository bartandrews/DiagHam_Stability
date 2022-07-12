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


#ifndef FQHEMPODENSITYOPERATOR_H
#define FQHEMPODENSITYOPERATOR_H


#include "config.h"
#include "Tools/FQHEMPS/AbstractFQHEMPOMatrix.h"


class FQHEMPODensityOperator : public AbstractFQHEMPOMatrix
{

 protected:

 public:
  
  // default constructor 
  //
  FQHEMPODensityOperator();

  // constructor 
  //
  // maxOccupation = maximum occupation for a single orbital
  // prefactor = prefactor in front of the density operator
  FQHEMPODensityOperator(int maxOccupation, double prefactor);

  // destructor
  //
  ~FQHEMPODensityOperator();
  
  // get the boundary indices of the MPO representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex);

  // get the name describing the O matrices 
  // 
  // return value = name 
  virtual char* GetName();

  // get the maximum occupation per orbital
  //
  // return value = aximum occupation per orbital
  virtual int GetMaximumOccupation();

  // set a new prefactor in front of the operator
  //
  // prefactor = prefactor in front of the operator
  virtual void SetPrefactor(double prefactor);

 protected:

};

// get the maximum occupation per orbital
//
// return value = aximum occupation per orbital

inline int FQHEMPODensityOperator::GetMaximumOccupation()
{
  return (this->NbrOMatrices - 1);
}

#endif
