////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class for n dimensional real vector                     //
//            whose memory allocation is done only on one process             //
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


#ifndef DELOCALIZEDREALVECTOR_H
#define DELOCALIZEDREALVECTOR_H


#ifdef USE_CLUSTER_ARCHITECTURE

#include "config.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"
#include "Architecture/ClusterArchitecture/AbstractClusterArchitecture.h"

#include <iostream>
#include <fstream>


using std::ostream;
using std::ifstream;
using std::ofstream;


class Complex;
class ComplexVector;
#ifdef USE_OUTPUT
class MathematicaOutput;
#endif
class BlockDiagonalMatrix;
class Matrix;
class RealMatrix;
class RealAntisymmetricMatrix;
class RealDiagonalMatrix;
class RealSymmetricMatrix;
class RealTriDiagonalSymmetricMatrix;
class ComplexVector;


class DelocalizedRealVector : public RealVector
{

 protected:
  
  // id of the process where the vector is localized
  int LocalizationId;
  // id of the local process
  int LocalId;
  // pointer to the cluster architecture in use
  AbstractClusterArchitecture* Architecture;

  // a dummy element used when passing cevtor element as reference
  double DummyElement;
  // pointer to an index of the component correponding to the dummy element
  int* DummyElementPosition;

 public:

  // default constructor
  //
  DelocalizedRealVector();

  // constructor for an empty real vector 
  //
  // size = Vector Dimension 
  // architecture = pointer to the cluster architecture in use
  // vectorId = id of the vector
  // localizationId = id of the process where the vector is localized
  // zeroFlag = true if all coordinates have to be set to zero
  DelocalizedRealVector(int size, AbstractClusterArchitecture* architecture, int vectorId, int localizationId, bool zeroFlag = false);

  // copy constructor
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  DelocalizedRealVector(const DelocalizedRealVector& vector, bool duplicateFlag = false);

  // copy constructor from a real vector
  //
  // vector = vector to copy
  // duplicateFlag = true if datas have to be duplicated
  DelocalizedRealVector(const RealVector& vector, AbstractClusterArchitecture* architecture, bool duplicateFlag = false);

  // copy constructor from a complex vector (keep only real part and datas are duplicated)
  //
  // vector = vector to copy
  DelocalizedRealVector(const ComplexVector& vector, AbstractClusterArchitecture* architecture);

  // copy constructor from a vector (duplicate datas if necessary)
  //
  // vector = vector to copy
  DelocalizedRealVector(const Vector& vector, AbstractClusterArchitecture* architecture);

  // destructor
  //
  ~DelocalizedRealVector ();

  // assignement
  //
  // vector = vector to assign
  // return value = reference on current vector
  DelocalizedRealVector& operator = (const DelocalizedRealVector& vector);

  // assignement from a real vector
  //
  // vector = vector to assign
  // return value = reference on current vector
  DelocalizedRealVector& operator = (const RealVector& vector);

  // assignement from a complex vector (keep only real part and datas are duplicated)
  //
  // vector = vector to assign
  // return value = reference on current vector
  DelocalizedRealVector& operator = (const ComplexVector& vector);

  // assignement from a vector (duplicate datas if necessary)
  //
  // vector = vector to assign
  // return value = reference on current vector
  DelocalizedRealVector& operator = (const Vector& vector);

  // Resize vector
  //
  // dimension = new dimension
  void Resize (int dimension);

  // Resize vector and set to zero all components that have been added
  //
  // dimension = new dimension
  void ResizeAndClean (int dimension);

  // copy a vector into another
  //
  // vector = vector to copy
  // coefficient = optional coefficient which multiply source to copy
  // return value = reference on current vector
  DelocalizedRealVector& Copy (RealVector& vector, double coefficient = 1.0);

  // create a new vector with same size and same type but without duplicating datas
  //
  // zeroFlag = true if all coordinates have to be set to zero
  // return value = pointer to new vector 
  Vector* EmptyClone(bool zeroFlag = false);

  // return vector i-th coordinate (without testing if position is valid)
  //
  // i = coordinate position
  double& operator [] (int i);

  // read vector from a file 
  //
  // fileName = name of the file where the vector has to be read
  // return value = true if no error occurs
  bool ReadVector (char* fileName);

  // input file stream overload
  //
  // file = reference on output file stream
  // vector = reference on vector to load
  // return value = reference on output file stream
  friend ifstream& operator >> (ifstream& file, DelocalizedRealVector& vector);

  // localize the current vector to the current process
  // 
  void Localize();

  // delocalize the current vector from the current process
  // 
  // transfertFlag = indicates if the current vector datas have to sent to the vector real location
  void Delocalize(bool transfertFlag = false);

 private: 

  // flush dummy element to take into account last change
  //
  void FlushDummyElement();

};

// return vector i-th coordinate (without testing if position is valid)
//
// i = coordinate position

inline double& DelocalizedRealVector::operator [] (int i)
{
  if (this->LocalizationId == this->LocalId)
    {
      return this->Components[i];
    }
  else
    {
      if ((*(this->DummyElementPosition)) != -1)
	{
	  this->Architecture->SetRealVectorElement(this->DummyElement, this->VectorId, (*(this->DummyElementPosition)));	  
	}
      this->DummyElement = this->Architecture->RequestRealVectorElement(this->VectorId, i);
      (*(this->DummyElementPosition)) = i;
      return this->DummyElement;
    }
}
 


// flush dummy element to take into account last change
//

inline void DelocalizedRealVector::FlushDummyElement()
{
  if ((this->LocalId != this->LocalizationId) && ((*(this->DummyElementPosition)) != -1))
    {
      this->Architecture->SetRealVectorElement(this->DummyElement, this->VectorId, (*(this->DummyElementPosition)));
    }
  (*(this->DummyElementPosition)) = -1;
}

#endif

#endif

