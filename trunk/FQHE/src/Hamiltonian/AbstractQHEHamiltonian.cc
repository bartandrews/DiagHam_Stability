////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of quantum Hall hamiltonian                     //
//                                                                            //
//                        last modification : 03/07/2003                      //
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
#include "Hamiltonian/AbstractQHEHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::cout;
using std::endl;


// destructor
//

AbstractQHEHamiltonian::~AbstractQHEHamiltonian()
{
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{  
  for (int i = 0; i < nbrVectors; ++i)
    {
      vDestinations[i].ClearVector();
    }
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors);
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							     int firstComponent, int nbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    {
      vDestinations[i].ClearVector();
    }
  return this->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vector sat which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->LowLevelAddMultiply(vSources[i], vDestinations[i]);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
								int firstComponent, int nbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    this->LowLevelAddMultiply(vSources[i], vDestinations[i], firstComponent, nbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

// ComplexVector* AbstractQHEHamiltonian::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
// {
//   for (int i = 0; i < nbrVectors; ++i)
//     this->LowLevelMultiply(vSources[i], vDestinations[i]);
//   return vDestinations;
// }

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

// ComplexVector* AbstractQHEHamiltonian::LowLevelMultipleMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
// 								int firstComponent, int nbrComponent)

// {
//   for (int i = 0; i < nbrVectors; ++i)
//     this->LowLevelMultiply(vSources[i], vDestinations[i], firstComponent, nbrComponent);
//   return vDestinations;
// }
  
// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored

// ComplexVector* AbstractQHEHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors)
// {
//   for (int i = 0; i < nbrVectors; ++i)
//     this->LowLevelAddMultiply(vSources[i], vDestinations[i]);
//   return vDestinations;
// }

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

//ComplexVector* AbstractQHEHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int n// brVectors, 
// 								   int firstComponent, int nbrComponent)
// {
//   for (int i = 0; i < nbrVectors; ++i)
//     this->LowLevelAddMultiply(vSources[i], vDestinations[i], firstComponent, nbrComponent);
//   return vDestinations;
// }


// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool AbstractQHEHamiltonian::SavePrecalculation (char* fileName)
{
  return false;
}

// evaluate all interaction factors
//   

void AbstractQHEHamiltonian::EvaluateInteractionFactors()
{
  return;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long AbstractQHEHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  return ((long) -1);
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  return ((long) -1);
}

// test the amount of memory needed for fast multiplication algorithm, and keep track of unique matrix entries (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element
long AbstractQHEHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent, 
					       SortedRealUniqueArray &realElements, 
					       SortedComplexUniqueArray &complexElements)
{
   return ((long) -1);
}
 
// enable fast multiplication algorithm
//

void AbstractQHEHamiltonian::EnableFastMultiplication()
{
  return;
}

// enable fast multiplication algorithm (partial evaluation)
//
// jobIndex = index of the job that proceeds part of the fast multiplication evaluation
// nbrJob = number of jobs that proceed the fast multiplication evaluation

void AbstractQHEHamiltonian::PartialEnableFastMultiplication(int jobIndex, int nbrJob)
{
  return;
}

// load precalculations from a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be read
// return value = true if no error occurs

bool AbstractQHEHamiltonian::LoadPrecalculation (char* fileName)
{
  return false;
}

