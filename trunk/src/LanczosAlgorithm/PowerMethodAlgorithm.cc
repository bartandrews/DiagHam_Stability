////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2003 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of power methof algorithm                      //
//                                                                            //
//                        last modification : 24/01/2013                      //
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


#include "LanczosAlgorithm/PowerMethodAlgorithm.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"

#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//

PowerMethodAlgorithm::PowerMethodAlgorithm()
{ 
}

// default constructor
//
// architecture = architecture to use for matrix operations
// maxIter = an approximation of maximal number of iteration

PowerMethodAlgorithm::PowerMethodAlgorithm(AbstractArchitecture* architecture, int maxIter)
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->MaximumNumberIteration = maxIter;
  this->Architecture = architecture;
  this->Flag.Initialize();
  this->PreviousEigenvalue = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->EigenvectorPrecision = 0.0;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

PowerMethodAlgorithm::PowerMethodAlgorithm(const PowerMethodAlgorithm& algorithm) 
{
  this->Index = algorithm.Index;
  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->CurrentVector = algorithm.CurrentVector;
  this->Flag = algorithm.Flag;
  this->Architecture = algorithm.Architecture;
  this->PreviousEigenvalue = algorithm.PreviousEigenvalue;
  this->CurrentEigenvalue = algorithm.CurrentEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->EigenvectorPrecision = algorithm.EigenvectorPrecision;
}

// destructor
//

PowerMethodAlgorithm::~PowerMethodAlgorithm() 
{
}

// initialize Lanczos algorithm with a random vector
//

void PowerMethodAlgorithm::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->CurrentVector = RealVector (Dimension);
  this->PreviousVector = RealVector (Dimension);
  for (int i = 0; i < Dimension; i++)
    {
      this->CurrentVector[i] = (rand() - 32767) * 0.5;
    }
  this->CurrentVector /= this->CurrentVector.Norm();
  this->Index = 0;
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void PowerMethodAlgorithm::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->PreviousVector = RealVector (Dimension);
  this->CurrentVector = vector;
  this->Index = 0;
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& PowerMethodAlgorithm::GetGroundState()
{
  return this->CurrentVector;
}

// get the n first eigenstates
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* PowerMethodAlgorithm::GetEigenstates(int nbrEigenstates)
{
  RealVector* Eigenstates = new RealVector [1];
  Eigenstates[0] = this->CurrentVector;
  return Eigenstates;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void PowerMethodAlgorithm::RunLanczosAlgorithm (int nbrIter) 
{
  if (this->Index == 0)
    {
      RealVector TmpVector;
      TmpVector = this->CurrentVector;
      this->CurrentVector = this->PreviousVector;
      this->PreviousVector = TmpVector;
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(this->PreviousVector), &(this->CurrentVector));
      Operation1.ApplyOperation(this->Architecture);      
      this->CurrentEigenvalue = (this->PreviousVector * this->CurrentVector);
      this->CurrentVector /= this->CurrentVector.Norm();
      this->Index = 1;
    }
  for (int i = 0; i < nbrIter; ++i)
    {
      RealVector TmpVector;
      TmpVector = this->CurrentVector;
      this->CurrentVector = this->PreviousVector;
      this->PreviousVector = TmpVector;
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(this->PreviousVector), &(this->CurrentVector));
      Operation1.ApplyOperation(this->Architecture);      
      this->PreviousEigenvalue = this->CurrentEigenvalue;
      this->CurrentEigenvalue = (this->PreviousVector * this->CurrentVector);
      this->CurrentVector /= this->CurrentVector.Norm();
      ++this->Index;
    }
}

  
// test if convergence has been reached
//
// return value = true if convergence has been reached

bool PowerMethodAlgorithm::TestConvergence ()
{
  cout << this->Index << " : " << this->CurrentEigenvalue << " " << fabs(this->CurrentEigenvalue) << " " << fabs(this->CurrentEigenvalue - this->PreviousEigenvalue) << endl;
  if (fabs(this->CurrentEigenvalue - this->PreviousEigenvalue) < (this->EigenvaluePrecision * fabs(this->CurrentEigenvalue)))
    {
      return true;
    }
  else
    return false;
}

// get the n first eigenvalues
//
// eigenvalues = reference on the array where the eigenvalues will be stored (allocation done by the method itself)
// nbrEigenstates = number of needed eigenvalues

void PowerMethodAlgorithm::GetEigenvalues (double*& eigenvalues, int nbrEigenvalues)
{
  eigenvalues = new double [1];
  eigenvalues[0] = this->CurrentEigenvalue;
}

// get the n first eigenvalues
//
// eigenvalues = reference on the array where the eigenvalues will be stored (allocation done by the method itself)
// nbrEigenstates = number of needed eigenvalues

void PowerMethodAlgorithm::GetEigenvalues (Complex*& eigenvalues, int nbrEigenvalues)
{
  eigenvalues = new Complex [1];
  eigenvalues[0] = this->CurrentEigenvalue;
}


