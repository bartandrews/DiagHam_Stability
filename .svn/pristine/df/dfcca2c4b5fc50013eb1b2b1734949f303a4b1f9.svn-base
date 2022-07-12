////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of abstract cluster architecture                  //
//                                                                            //
//                        last modification : 15/06/2004                      //
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


#ifndef ABSTRACTCLUSTERARCHITECTURE_H
#define ABSTRACTCLUSTERARCHITECTURE_H


#include "config.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/Complex.h"


class AbstractClusterArchitecture : public AbstractArchitecture
{

 private:


 public:
  
  // destructor
  //
  virtual ~AbstractClusterArchitecture();
  
  // get id of the local process
  // 
  // return value = local process id
  virtual int GetProcessId() = 0;

  // get a vector element from a real vector
  //
  // vectorId = id of vector where the lement is located
  // index = element index in the vector
  // return value = corresponding vector element
  virtual double RequestRealVectorElement(int vectorId, int index) = 0;
  
  // set a vector element into a real vector
  //
  // component = value of the vector element
  // vectorId = id of vector where the element is located
  // index = element index in the vector
  virtual void SetRealVectorElement(const double& component, int vectorId, int index) = 0;
  
  // get a vector element from a complex vector
  //
  // vectorId = id of vector where the lement is located
  // index = element index in the vector
  // return value = corresponding vector element
  virtual Complex RequestComplexVectorElement(int vectorId, int index) = 0;
  
  // set a vector element into a real vector
  //
  // component = value of the vector element
  // vectorId = id of vector where the lement is located
  // index = element index in the vector
  virtual void SetComplexVectorElement(const Complex& component, int vectorId, int index) = 0;
  
  // get a real vector 
  // 
  // vectorId = id of vector to get
  // return value = corresponding vector
  virtual RealVector GetRealVector(int vectorId) = 0;
  
  // set a real vector 
  // 
  // vector = reference on the vector which contains thedatas 
  // vectorId = id of vector where to store the datas
  virtual void SetRealVector(RealVector& vector, int vectorId) = 0;
  
  // get square norm of a given vector
  // 
  // vectorId = vector id 
  // return value = square norm
  //  virtual double GetVectorSqrNorm(int vectorId) = 0;
  
  // normalize a given vector
  // 
  // vectorId = if of the vector to normalize
  //  virtual void NormalizeVector(int vectorId) = 0;
  
};

#endif
