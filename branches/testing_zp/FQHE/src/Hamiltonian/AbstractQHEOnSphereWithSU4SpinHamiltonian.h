////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract fractional quantum Hall hamiltonian          //
//              associated to particles with SU(4) spin on a sphere           //
//                                                                            //
//                        last modification : 28/11/2006                      //
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


#ifndef ABSTRACTQHEONSPHEREWITHSU4SPINHAMILTONIAN_H
#define ABSTRACTQHEONSPHEREWITHSU4SPINHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"

#include <iostream>


using std::ostream;


class AbstractQHEOnSphereWithSU4SpinHamiltonian : public AbstractQHEOnSphereHamiltonian
{

 protected:
  
  // number of index sum for the creation (or annhilation) operators for the intra spin/isopin sector
  int NbrIntraSectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the intra spin/isopin sector
  int* NbrIntraSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators for the intra spin/isopin sector
  int** IntraSectorIndicesPerSum;

  // number of index sum for the creation (or annhilation) operators for the inter spin/isopin sector
  int NbrInterSectorSums;
  // array containing the number of index group per index sum for the creation (or annhilation) operators for the inter spin/isopin sector
  int* NbrInterSectorIndicesPerSum;
  // array containing the (m1,m2) indices per index sum for the creation (or annhilation) operators for the inter spin/isopin sector
  int** InterSectorIndicesPerSum;

  // arrays containing all interaction factors, the first index correspond correspond to index sum for the creation (or annhilation) operators
  // the second index is a linearized index (m1,m2) + (n1,n2) * (nbr element in current index sum) (m for creation operators, n for annhilation operators)
  // array containing all interaction factors for spin up,isospin plus and spin up,isospin plus
  double** InteractionFactorsupup;
  // array containing all interaction factors for spin up,isospin plus and spin up,isospin minus
  double** InteractionFactorsupum;
  // array containing all interaction factors for spin up,isospin plus and spin down,isospin plus
  double** InteractionFactorsupdp;
  // array containing all interaction factors for spin up,isospin plus and spin down,isospin minus
  double** InteractionFactorsupdm;
  // array containing all interaction factors for spin up,isospin minus and spin up,isospin minus
  double** InteractionFactorsumum;
  // array containing all interaction factors for spin up,isospin minus and spin down,isospin plus
  double** InteractionFactorsumdp;
  // array containing all interaction factors for spin down,isospin minus and spin down,isospin minus
  double** InteractionFactorsumdm;
  // array containing all interaction factors for spin down,isospin plus and spin down,isospin plus
  double** InteractionFactorsdpdp;
  // array containing all interaction factors for spin down,isospin plus and spin down,isospin minus
  double** InteractionFactorsdpdm;
  // array containing all interaction factors for spin down,isospin minus and spin down,isospin minus
  double** InteractionFactorsdmdm;

  // array that contains all one-body interaction factors for particles with spin up/isospin plus
  double* OneBodyInteractionFactorsupup;
  // array that contains all one-body interaction factors for particles with spin up/isospin minus
  double* OneBodyInteractionFactorsumum;
  // array that contains all one-body interaction factors for particles with spin down/isospin plus
  double* OneBodyInteractionFactorsdpdp;
  // array that contains all one-body interaction factors for particles with spin down/isospin minus
  double* OneBodyInteractionFactorsdmdm;

 public:

  // destructor
  //
  virtual ~AbstractQHEOnSphereWithSU4SpinHamiltonian() = 0;

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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
  virtual  RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						   int firstComponent, int nbrComponent);

 protected:
 
  // multiply a et of vectors by the current hamiltonian for a given range of indices 
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

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication();

  // enable fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int lastComponent);

  // enable fast multiplication algorithm
  //
  // fileName = prefix of the name of the file where temporary matrix elements will be stored
  virtual void EnableFastMultiplicationWithDiskStorage(char* fileName);

};

#endif
