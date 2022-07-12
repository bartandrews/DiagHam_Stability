////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of factorial coefficient                      //
//                                                                            //
//                        last modification : 04/06/2002                      //
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


#ifndef NSPHEREPARAMETERS_H
#define NSPHEREPARAMETERS_H


#include "config.h"

#include <iostream>
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"


class NSphereParameters
{

  // dimension
  int Dimension;

  // complex flag
  bool IsComplex;

  // number of real parameters
  int NbrParameters;

  // vector to hold complex coordinates
  ComplexVector ComplexCoordinates;

  // vector to hold real coordinates
  RealVector RealCoordinates;

  // vector to hold parameters
  RealVector Parameters;

  // table holding temporary data
  double *CosTable;
  double *SinTable;
  
 public:

  // default constructor 
  //
  NSphereParameters();

  // constructor from an integer
  //
  // dim = dimension of sphere
  // isComplex = flag indicating whether complex N-Sphere is used
  NSphereParameters(int dim, bool isComplex);

  // destructor
  //
  ~NSphereParameters();

  // assignement
  //
  // factorial = factorial coefficient to assign
  // return value = reference on current factorial coefficient
  NSphereParameters& operator = (const NSphereParameters& sphere);

  // get number of parameters
  int GetNbrParameters(){return this->NbrParameters;}

  // set parameters
  void SetParameters(double *parameters);
  
  // get complex coordinates for last parameter set
  ComplexVector& GetComplexCoordinates();
  
  // get complex coordinates for last parameter set
  RealVector& GetRealCoordinates();
  

 private:


};

// get complex coordinates for last parameter set
inline ComplexVector& NSphereParameters::GetComplexCoordinates()
{
  return this->ComplexCoordinates;
}
  
// get complex coordinates for last parameter set
inline RealVector& NSphereParameters::GetRealCoordinates()
{
  return this->RealCoordinates;
}


#endif


