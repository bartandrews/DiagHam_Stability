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
#include "Hamiltonian/ProjectedQHEHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"
#include "LanczosAlgorithm/InternalReorthogonalizedLanczosAlgorithm.h"
#include "MathTools/Complex.h"

#include <iostream>

using std::ostream;  


// constructor
// hamiltonian = core Hamiltonian to be diagonalized
// projector = hamiltonian onto whose groundstate subspace should be projected
// architecture = architecture to use for precalculation and multiplications
// maxNbrVectorsProj = largest number of eigenvalues stored in memory for projector Lanczos
// maxIterProj = largest number of Lanczos iterations to be used in projection
// projectorPrecision = Lanczos precision required for projection
ProjectedQHEHamiltonian::ProjectedQHEHamiltonian(AbstractQHEHamiltonian *hamiltonian, AbstractQHEHamiltonian *projector,
						 AbstractArchitecture* architecture, int maxNbrVectorsProj,
						 int maxIterProj, double projectorPrecision)
{
  this->Hamiltonian=hamiltonian;
  this->Projector=projector;
  this->Flag.Initialize();
  this->LanczosAlgorithm = new InternalReorthogonalizedLanczosAlgorithm(architecture, 1, maxNbrVectorsProj, maxIterProj);
  this->LanczosAlgorithm->SetHamiltonian(this->Projector);
  if (projectorPrecision != 0.0)
	LanczosAlgorithm->SetEigenvaluePrecision(projectorPrecision);
}

// destructor
//
ProjectedQHEHamiltonian::~ProjectedQHEHamiltonian()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete this->LanczosAlgorithm;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use
void ProjectedQHEHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Hamiltonian->SetHilbertSpace(hilbertSpace);
  this->Projector->SetHilbertSpace(hilbertSpace);
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space
AbstractHilbertSpace* ProjectedQHEHamiltonian::GetHilbertSpace ()
{
  return this->Hamiltonian->GetHilbertSpace();
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension
int ProjectedQHEHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Hamiltonian->GetHilbertSpaceDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value
void ProjectedQHEHamiltonian::ShiftHamiltonian (double shift)
{
  this->Hamiltonian->ShiftHamiltonian(shift);
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element
Complex ProjectedQHEHamiltonian::MatrixElement (RealVector& V1, RealVector& V2)
{
  return Hamiltonian->MatrixElement(V1,V2);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element
Complex ProjectedQHEHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Hamiltonian->MatrixElement(V1,V2);
}


// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored
RealVector& ProjectedQHEHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination)
{
  this->Hamiltonian->LowLevelMultiply(vSource, vDestination);
  this->LanczosAlgorithm->ProjectVector(vDestination);
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
RealVector& ProjectedQHEHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination,
				       int firstComponent, int nbrComponent)
{
  this->Hamiltonian->LowLevelMultiply(vSource, vDestination, firstComponent, nbrComponent);
  this->LanczosAlgorithm->ProjectVector(vDestination);
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored
RealVector& ProjectedQHEHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
				       int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  this->Hamiltonian->LowLevelMultiply(vSource,vDestination, sourceStart, sourceStep, sourceShift, sourceNbrComponent,
					 destinationStart, destinationStep, destinationShift, destinationNbrComponent);
  this->LanczosAlgorithm->ProjectVector(vDestination);
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored
RealVector& ProjectedQHEHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  this->Hamiltonian->LowLevelAddMultiply(vSource, vDestination);
  this->LanczosAlgorithm->ProjectVector(vDestination);
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
RealVector& ProjectedQHEHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					  int firstComponent, int nbrComponent)
{
  this->Hamiltonian->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  this->LanczosAlgorithm->ProjectVector(vDestination);
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// sourceStart = source vector first index
// sourceStep = step to add to go to the following source vector index
// sourceShift = shift to apply when directly accessing source vector component (must be substracted to the real index)
// sourceNbrComponent = number of component to take into account in the source vector
// destinationStart = destination vector first index
// destinationStep = step to add to go to the following destination vector index
// destinationShift = shift to apply when directly accessing destination vector component (must be substracted to the real index)
// destinationNbrComponent = number of component to take into account in the destination vector
// return value = reference on vector where result has been stored
RealVector& ProjectedQHEHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
					int sourceStart, int sourceStep, int sourceShift, int sourceNbrComponent,
					int destinationStart, int destinationStep, int destinationShift, int destinationNbrComponent)
{
  this->Hamiltonian->LowLevelAddMultiply(vSource, vDestination, sourceStart, sourceStep, sourceShift, sourceNbrComponent,
					 destinationStart, destinationStep, destinationShift, destinationNbrComponent);
  this->LanczosAlgorithm->ProjectVector(vDestination);
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian and store result in another set of vectors
// low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored
RealVector* ProjectedQHEHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  this->Hamiltonian->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors);
  for (int i=0; i<nbrVectors; ++i)
    {
      this->LanczosAlgorithm->ProjectVector(vDestinations[i]);
    }
  return vDestinations;
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
RealVector* ProjectedQHEHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
					       int firstComponent, int nbrComponent)
{
  this->Hamiltonian->LowLevelMultipleMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
  for (int i=0; i<nbrVectors; ++i)
    {
      this->LanczosAlgorithm->ProjectVector(vDestinations[i]);
    }
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vector sat which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// return value = pointer to the array of vectors where result has been stored
RealVector* ProjectedQHEHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors)
{
  this->Hamiltonian->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors);
  for (int i=0; i<nbrVectors; ++i)
    {
      this->LanczosAlgorithm->ProjectVector(vDestinations[i]);
    }
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
RealVector* ProjectedQHEHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent)
{
  this->Hamiltonian->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
  for (int i=0; i<nbrVectors; ++i)
    {
      this->LanczosAlgorithm->ProjectVector(vDestinations[i]);
    }
  return vDestinations;
}

 
// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs
bool ProjectedQHEHamiltonian::SavePrecalculation (char* fileName)
{
  return false;
}

// evaluate all interaction factors
//   
void ProjectedQHEHamiltonian::EvaluateInteractionFactors()
{
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed
long ProjectedQHEHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  return 0l;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element
long ProjectedQHEHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  return 0l;
}

// enable fast multiplication algorithm
//
void ProjectedQHEHamiltonian::EnableFastMultiplication()
{  
}

// enable fast multiplication algorithm (partial evaluation)
//
// jobIndex = index of the job that proceeds part of the fast multiplication evaluation
// nbrJob = number of jobs that proceed the fast multiplication evaluation
void ProjectedQHEHamiltonian::PartialEnableFastMultiplication(int jobIndex, int nbrJob)
{
}

// load precalculations from a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be read
// return value = true if no error occurs
bool ProjectedQHEHamiltonian::LoadPrecalculation (char* fileName)
{
  return false;
}

