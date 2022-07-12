////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class implementing the linear superposition of several Hamiltonians    //
//                           class author: Gunnar Möller                      //
//                                                                            //
//                        last modification : 02/03/2007                      //
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


#ifndef LINEARLYSUPERPOSEDQHEONSPHEREHAMILTONIAN_H
#define LINEARLYSUPERPOSEDQHEONSPHEREHAMILTONIAN_H


#include "AbstractQHEHamiltonian.h"
#include "AbstractQHEOnSphereHamiltonian.h"
#include "Vector/RealVector.h"


class LinearlySuperposedQHEOnSphereHamiltonian: public AbstractQHEOnSphereHamiltonian
{

 protected:
 
  // number of Hamiltonians
  int NumHamiltonians;
  // array that contains the cofficient in front of each Hamiltonian
  double *LinearCoefficients;
  //  array of hamiltonians
  AbstractQHEOnSphereHamiltonian **ListHamiltonians;

  // temporary vectors
  RealVector VectorSum;
  RealVector VectorSummand;

 public:

  // creates a Hamiltonian that consists of a superposition of Number Hamiltonians with an array of coefficients and
  // an array of Hamiltonians
  //
  // number = number of Hamiltonians
  // coefficients = array that contains the cofficient in front of each Hamiltonian
  // hamiltonians = array of hamiltonians
  LinearlySuperposedQHEOnSphereHamiltonian(int number, double *coefficients, AbstractQHEOnSphereHamiltonian **hamiltonians);
  
  virtual ~LinearlySuperposedQHEOnSphereHamiltonian();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  
  void ShiftHamiltonian (double shift);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored

  RealVector& LowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
						     int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored

  RealVector& LowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
					     int firstComponent, int nbrComponent);


  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored

  RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
					  int firstComponent, int nbrComponent);


 protected:

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored

  RealVector* LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							     int firstComponent, int nbrComponent);

  
  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  
  RealVector* LowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						     int firstComponent, int nbrComponent);

  // test the amount of memory needed for fast multiplication algorithm
  //
  // allowedMemory = amount of memory that cam be allocated for fast multiplication
  // return value = amount of memory needed

  long FastMultiplicationMemory(long allowedMemory);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  
  long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);
  

  // enable fast multiplication algorithm
  //

  void EnableFastMultiplication();

  // enable fast multiplication algorithm using on disk cache 
  //
  // fileName = prefix of the name of the file where temporary matrix elements will be stored
  
  void EnableFastMultiplicationWithDiskStorage(char* fileName);

  // save precalculations in a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be stored
  // return value = true if no error occurs
  
  bool SavePrecalculation (char* fileName);

  // load precalculations from a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be read
  // return value = true if no error occurs

  bool LoadPrecalculation (char* fileName);


  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();
    
};

#endif  // LINEARLYSUPERPOSEDQHEONSPHEREHAMILTONIAN_H
