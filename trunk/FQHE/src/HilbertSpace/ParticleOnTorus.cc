////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of particle on a torus                       //
//                                                                            //
//                        last modification : 18/07/2002                      //
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
#include "HilbertSpace/ParticleOnTorus.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/FQHETorusApplyCNRotationOperation.h"
#include "Architecture/ArchitectureOperation/FQHETorusSymmetrizeU1U1StateOperation.h"


// virtual destructor
//

ParticleOnTorus::~ParticleOnTorus ()
{
}


// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space
void ParticleOnTorus::SetTargetSpace(ParticleOnTorus* targetSpace)
{
  printf("Attention - trying to call SetTargetSpace on generic class ParticleOnTorus\n");
}


// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space
void ParticleOnTorus::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->SetTargetSpace((ParticleOnTorus*) targetSpace);
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// kySector = Ky sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnTorus::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, RealVector& groundState)
{
  RealSymmetricMatrix TmpDensityMatrix;
  return TmpDensityMatrix;
}
// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// kySector = Ky sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix ParticleOnTorus::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, ComplexVector& groundState)
{
  HermitianMatrix TmpDensityMatrix;
  return TmpDensityMatrix;
}

// apply a magnetic translation along x to a given state
//
// index = state index 
// return value = translated state index

int ParticleOnTorus::ApplyXMagneticTranslation(int index)
{
  return -1;
}
  
// apply a magnetic translation along x to a given state
//
// index = state index 
// sign = additional sign due to the particle statistics
// return value = translated state index

int ParticleOnTorus::ApplyXMagneticTranslation(int index, double& sign)
{
  sign = 1.0;
  return this->ApplyXMagneticTranslation(index);
}

// remove part of each Fock state, discarding component if the Fock state does not a given pattern
//
// inputVector = state to truncate
// reducedSpace = Hilbert space where the truncated state will lie
// pattern = array describing the pattern 
// patternSize = pattern size
// patternShift = indicate where the pattern has to be applied
// return value = trucated state

RealVector ParticleOnTorus::TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnTorus* reducedSpace, int* pattern, int patternSize, int patternShift)
{
  RealVector Tmp;
  return Tmp;
}

// apply the C4 rotation to a given state assumin it is an eigenstate of both kx and ky
//
// inputState = reference on the state that has to be rotated
// inputSpace = Hilbert space associated to the input state
// architecture = pointer to the architecture
// clockwise = the rotation is done clockwise
// return value = rotated state

ComplexVector ParticleOnTorus::C4Rotation (ComplexVector& inputState, ParticleOnTorus* inputSpace, bool clockwise, AbstractArchitecture* architecture)
{
  ComplexVector OutputState (this->HilbertSpaceDimension, true);
  if (architecture != 0)
    {
      if (clockwise == false)
	{
	  FQHETorusApplyCNRotationOperation Operation(4, &inputState, &OutputState, inputSpace, this);
	  Operation.ApplyOperation(architecture);
	}
      else
	{
	  FQHETorusApplyCNRotationOperation Operation(-4, &inputState, &OutputState, inputSpace, this);
	  Operation.ApplyOperation(architecture);
	}
    }
  else
    {
      this->CoreC4Rotation(inputState, inputSpace, OutputState, 0 , this->HilbertSpaceDimension, clockwise);
    }
  return OutputState;
}

// core part of the C4 rotation
//
// inputState = reference on the state that has to be rotated
// inputSpace = Hilbert space associated to the input state
// outputState = reference on the rotated state
// minIndex = minimum index that has to be computed
// nbrIndices = number of indices that have to be computed
// clockwise = the rotation is done clockwise
// return value = reference on the rotated state

ComplexVector& ParticleOnTorus::CoreC4Rotation (ComplexVector& inputState, ParticleOnTorus* inputSpace, ComplexVector& outputState, int minIndex, int nbrIndices, bool clockwise)
{
  return outputState;
}

// symmetrized a product of two uncoupled states 
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

RealVector ParticleOnTorus::SymmetrizeU1U1State (RealVector& leftVector, RealVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, AbstractArchitecture* architecture)
{
  RealVector SymmetrizedVector (this->LargeHilbertSpaceDimension,true);
  FQHETorusSymmetrizeU1U1StateOperation Operation (this, leftSpace, rightSpace, &SymmetrizedVector, &leftVector, &rightVector);
  Operation.ApplyOperation(architecture);
  return SymmetrizedVector;
}
  
// symmetrized a product of two uncoupled states 
//
// outputVector = reference on the vector which will contain the symmetrozed state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

ComplexVector ParticleOnTorus::SymmetrizeU1U1State (ComplexVector& leftVector, ComplexVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, AbstractArchitecture* architecture)
{
  ComplexVector SymmetrizedVector(this->LargeHilbertSpaceDimension, true);
//   timeval TotalStartingTime;
//   gettimeofday (&TotalStartingTime, 0);
  FQHETorusSymmetrizeU1U1StateOperation Operation (this, leftSpace, rightSpace, &SymmetrizedVector, &leftVector, &rightVector);
  Operation.ApplyOperation(architecture);
//   unsigned long firstComponent = 0ul;
//   unsigned long nbrComponent = leftSpace->GetLargeHilbertSpaceDimension();
//   this->SymmetrizeU1U1StateCore (SymmetrizedVector, leftVector, rightVector,  leftSpace,  rightSpace, unnormalizedBasisFlag, firstComponent, nbrComponent);
//   timeval TotalEndingTime;
//   gettimeofday (&TotalEndingTime, 0);
//   double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
//   cout << this->FirstComponent << " " <<  this->NbrComponent << " : " << Dt << "s" << endl;
//   cout << SymmetrizedVector.Norm() << endl;
//   if (unnormalizedBasisFlag == false)
//     SymmetrizedVector /= SymmetrizedVector.Norm();
  return SymmetrizedVector;
}
  
// symmetrized a product of several uncoupled states 
//
// inputStates = states which will be symmetrized
// inputSpaces = Hilbert spaces attached to each states
// nbrStates = number of states to symmetrize
// architecture = pointer to the architecture 
// return value = symmetrized state

RealVector ParticleOnTorus::SymmetrizeU1U1State (RealVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, double precision, AbstractArchitecture* architecture)
{
  ComplexVector* TmpVectors = new ComplexVector[nbrStates];
  for (int i = 0 ; i < nbrStates; ++i)
    {
      TmpVectors[i] = ComplexVector(inputStates[i]);
    }
  ComplexVector TmpOutputState (this->LargeHilbertSpaceDimension, true); 
  unsigned long FirstComponent = 0;
  unsigned long NbrComponents = inputSpaces[0]->GetHilbertSpaceDimension();
  this->SymmetrizeU1U1StateCore(TmpOutputState, TmpVectors, inputSpaces, nbrStates, FirstComponent, NbrComponents);
  delete[] TmpVectors;
  return RealVector(TmpOutputState);
}

// symmetrize a product of several uncoupled states 
//
// inputStates = states which will be symmetrized
// inputSpaces = Hilbert spaces attached to each states
// nbrStates = number of states to symmetrize
// architecture = pointer to the architecture 
// return value = symmetrized state

ComplexVector ParticleOnTorus::SymmetrizeU1U1State (ComplexVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, double precision,  AbstractArchitecture* architecture)
{
  ComplexVector SymmetrizedVector (this->LargeHilbertSpaceDimension, true);
  unsigned long FirstComponent = 0;
  unsigned long NbrComponents = inputSpaces[0]->GetHilbertSpaceDimension();
  this->SymmetrizeU1U1StateCore (SymmetrizedVector, inputStates, inputSpaces, nbrStates, FirstComponent, NbrComponents);
  double TmpNorm = SymmetrizedVector.Norm();
  cout << "The norm of the symmetrized vector is " << TmpNorm << endl;
  if (TmpNorm > precision)
    SymmetrizedVector /= TmpNorm;
  return SymmetrizedVector;
}

// symmetrize a vector by grouping neighbouring orbitals
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
// kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
// architecture = pointer to the architecture
// return value = symmetrized state

int ParticleOnTorus::SymmetrizeSingleStateGroupingNeighbouringOrbitals (RealVector& inputVector, int nbrOrbitals, RealVector*& symmetrizedVectors, int*& kySectors, AbstractArchitecture* architecture, double precision)
{
  int TargetSpaceNbrOrbitals = this->GetNbrOrbitals() / nbrOrbitals;
  ComplexVector* TmpVectors = new ComplexVector[TargetSpaceNbrOrbitals];
  ComplexVector TmpInputVector(inputVector, true);
  this->SymmetrizeSingleStateGroupingNeighbouringOrbitalsCore(TmpInputVector, TmpVectors, nbrOrbitals, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  ++NbrGeneratedSectors;
	}     
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new RealVector[NbrGeneratedSectors];
  kySectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  double TmpNorm = TmpVectors[i].Norm();
	  cout << "The norm of the symmetrized vector is " << TmpNorm << endl;
	  if (TmpNorm > precision)
	    {
	      TmpVectors[i] /= TmpNorm;
	      symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i];
	      kySectors[NbrGeneratedSectors] = i;
	      ++NbrGeneratedSectors;	  
	    }
	}     
    }  
  delete[] TmpVectors;
  return NbrGeneratedSectors;
}

// symmetrize a vector by grouping neighbouring orbitals
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
// kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
// architecture = pointer to the architecture
// return value = symmetrized state

int ParticleOnTorus::SymmetrizeSingleStateGroupingNeighbouringOrbitals (ComplexVector& inputVector, int nbrOrbitals, ComplexVector*& symmetrizedVectors, int*& kySectors, AbstractArchitecture* architecture, double precision)
{
  int TargetSpaceNbrOrbitals = this->GetNbrOrbitals() / nbrOrbitals;
  ComplexVector* TmpVectors = new ComplexVector[TargetSpaceNbrOrbitals];
  this->SymmetrizeSingleStateGroupingNeighbouringOrbitalsCore(inputVector, TmpVectors, nbrOrbitals, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  ++NbrGeneratedSectors;
	}     
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new ComplexVector[NbrGeneratedSectors];
  kySectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  double TmpNorm = TmpVectors[i].Norm();
	  cout << "The norm of the symmetrized vector is " << TmpNorm << endl;
	  if (TmpNorm > precision)
	    {
	      TmpVectors[i] /= TmpNorm;
	      symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i];
	      kySectors[NbrGeneratedSectors] = i;
	      ++NbrGeneratedSectors;	  
	    }
	}     
    }  
  return NbrGeneratedSectors;
}

// symmetrize a vector by grouping distant and equally separated orbitals
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
// kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
// architecture = pointer to the architecture
// return value = symmetrized state

int ParticleOnTorus::SymmetrizeSingleStateGroupingDistantOrbitals (RealVector& inputVector, int nbrOrbitals, RealVector*& symmetrizedVectors, int*& kySectors, AbstractArchitecture* architecture, double precision)
{
  int TargetSpaceNbrOrbitals = this->GetNbrOrbitals() / nbrOrbitals;
  ComplexVector* TmpVectors = new ComplexVector[TargetSpaceNbrOrbitals];
  ComplexVector TmpInputVector(inputVector, true);
//   for (int j = 0; j < TmpInputVector.GetVectorDimension(); ++j)
//     cout << TmpInputVector[j] << endl;
  this->SymmetrizeSingleStateGroupingDistantOrbitalsCore(TmpInputVector, TmpVectors, nbrOrbitals, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  ++NbrGeneratedSectors;
	}     
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new RealVector[NbrGeneratedSectors];
  kySectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)
	{
	  double TmpNorm = TmpVectors[i].Norm();
// 	  for (int j = 0; j < TmpVectors[i].GetVectorDimension(); ++j)
// 	    cout << TmpVectors[i][j] << " " ;
// 	  cout << endl;
	  cout << "The norm of the symmetrized vector is " << TmpNorm << endl;
	  if (TmpNorm > precision)
	    {
	      TmpVectors[i] /= TmpNorm;
	      symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i];
	      kySectors[NbrGeneratedSectors] = i;
	      ++NbrGeneratedSectors;	  
	    }
	}     
    }  
  delete[] TmpVectors;
  return NbrGeneratedSectors;
}

// symmetrize a vector by grouping distant and equally separated orbitals
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
// kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
// architecture = pointer to the architecture
// return value = symmetrized state

int ParticleOnTorus::SymmetrizeSingleStateGroupingDistantOrbitals (ComplexVector& inputVector, int nbrOrbitals, ComplexVector*& symmetrizedVectors, int*& kySectors, AbstractArchitecture* architecture, double precision, bool twistedTorus)
{
  int TargetSpaceNbrOrbitals = this->GetNbrOrbitals() / nbrOrbitals;
  ComplexVector* TmpVectors = new ComplexVector[TargetSpaceNbrOrbitals];
  this->SymmetrizeSingleStateGroupingDistantOrbitalsCore(inputVector, TmpVectors, nbrOrbitals, 0ul, this->LargeHilbertSpaceDimension, twistedTorus);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  ++NbrGeneratedSectors;
	}     
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new ComplexVector[NbrGeneratedSectors];
  kySectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)
	{
	  double TmpNorm = TmpVectors[i].Norm();
	  cout << "The norm of the symmetrized vector is " << TmpNorm << endl;
	  if (TmpNorm > precision)
	    {
	      TmpVectors[i] /= TmpNorm;
	      symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i];
	      kySectors[NbrGeneratedSectors] = i;
	      ++NbrGeneratedSectors;	  
	    }
	}     
    }  
  return NbrGeneratedSectors;
}


// symmetrize a vector by keeping only a subset of equally separated orbitals
//
// inputVector = reference on the vector to symmetrize
// firstOrbitalIndex = index of the first orbital to keep
// periodicity = momentum periodicity 
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest number of particles to the largest 
//                      number of particles and the smallest Ky to the largest Ky
// nbrParticlesSectors = reference on the array on the particle number sectors that have been generated through the symmetrization procedure
// kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
// architecture = pointer to the architecture
// return value = number of states that have been generated through the symmetrization procedure

int ParticleOnTorus::SymmetrizeSingleStatePeriodicSubsetOrbitals (RealVector& inputVector, int firstOrbitalIndex, int periodicity, 
								    RealVector*& symmetrizedVectors, int*& nbrParticlesSectors, int*& kySectors, 
								    AbstractArchitecture* architecture)
{
  double NormError = MACHINE_PRECISION;
  int TargetSpaceNbrOrbitals = this->GetNbrOrbitals() / periodicity;
  ComplexVector TmpInputVector(inputVector, true);
  ComplexVector** TmpVectors = new ComplexVector*[this->GetNbrParticles() + 1];
  for (int i = 0; i <= this->GetNbrParticles(); ++i)
    {
      TmpVectors[i] = new ComplexVector[TargetSpaceNbrOrbitals];
    }
  this->SymmetrizeSingleStatePeriodicSubsetOrbitalCore(TmpInputVector, TmpVectors, firstOrbitalIndex, periodicity, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i <= this->GetNbrParticles(); ++i)
    {
      for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
	{
	  if (TmpVectors[i][j].GetVectorDimension() != 0)	
	    {
	      ++NbrGeneratedSectors;
	    }     
	}
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new RealVector[NbrGeneratedSectors];
  kySectors = new int[NbrGeneratedSectors];
  nbrParticlesSectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i <= this->GetNbrParticles(); ++i)
    {
      for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
	{
	  if (TmpVectors[i][j].GetVectorDimension() != 0)	
	    {
	      double TmpNorm = TmpVectors[i][j].Norm();
	      if (TmpNorm > NormError)
		{
		  TmpVectors[i][j] /= TmpNorm; 
		  symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i][j];
		  kySectors[NbrGeneratedSectors] = j;
		  nbrParticlesSectors[NbrGeneratedSectors] = i;
		  ++NbrGeneratedSectors;	  
		}
	    }     
	}
    }  
  return NbrGeneratedSectors;
}

// symmetrize a vector by keeping only a subset of equally separated orbitals
//
// inputVector = reference on the vector to symmetrize
// firstOrbitalIndex = index of the first orbital to keep
// periodicity = momentum periodicity 
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest number of particles to the largest 
//                      number of particles and the smallest Ky to the largest Ky
// nbrParticlesSectors = reference on the array on the particle number sectors that have been generated through the symmetrization procedure
// kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
// architecture = pointer to the architecture
// return value = number of states that have been generated through the symmetrization procedure

int ParticleOnTorus::SymmetrizeSingleStatePeriodicSubsetOrbitals (ComplexVector& inputVector, int firstOrbitalIndex, int periodicity, 
								    ComplexVector*& symmetrizedVectors, int*& nbrParticlesSectors, int*& kySectors, 
								    AbstractArchitecture* architecture)
{
  double NormError = MACHINE_PRECISION;
  int TargetSpaceNbrOrbitals = this->GetNbrOrbitals() / periodicity;
  ComplexVector** TmpVectors = new ComplexVector*[this->GetNbrParticles() + 1];
  for (int i = 0; i <= this->GetNbrParticles(); ++i)
    {
      TmpVectors[i] = new ComplexVector[TargetSpaceNbrOrbitals];
    }
  this->SymmetrizeSingleStatePeriodicSubsetOrbitalCore(inputVector, TmpVectors, firstOrbitalIndex, periodicity, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i <= this->GetNbrParticles(); ++i)
    {
      for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
	{
	  if (TmpVectors[i][j].GetVectorDimension() != 0)	
	    {
	      ++NbrGeneratedSectors;
	    }     
	}
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new ComplexVector[NbrGeneratedSectors];
  kySectors = new int[NbrGeneratedSectors];
  nbrParticlesSectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i <= this->GetNbrParticles(); ++i)
    {
      for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
	{
	  if (TmpVectors[i][j].GetVectorDimension() != 0)	
	    {
	      double TmpNorm = TmpVectors[i][j].Norm();
	      if (TmpNorm > NormError)
		{
		  TmpVectors[i][j] /= TmpNorm; 
		  symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i][j];
		  kySectors[NbrGeneratedSectors] = j;
		  nbrParticlesSectors[NbrGeneratedSectors] = i;
		  ++NbrGeneratedSectors;	  
		}     
	    }
	}
    }  
  return NbrGeneratedSectors;
}

// symmetrize a vector by keeping only a subset of equally separated orbitals
//
// inputVector = reference on the vector to symmetrize
// firstOrbitalIndex = index of the first orbital to keep
// periodicity = momentum periodicity 
// phase = an optional phase (in pi units) that can be added for each kept and occupied orbital
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest number of particles to the largest 
//                      number of particles and the smallest Ky to the largest Ky
// nbrParticlesSectors = reference on the array on the particle number sectors that have been generated through the symmetrization procedure
// kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
// architecture = pointer to the architecture
// return value = number of states that have been generated through the symmetrization procedure

int ParticleOnTorus::SymmetrizeSingleStatePeriodicSubsetOrbitals (ComplexVector& inputVector, int firstOrbitalIndex, int periodicity, double phase, 
								  ComplexVector*& symmetrizedVectors, int*& nbrParticlesSectors, int*& kySectors, 
								  AbstractArchitecture* architecture)
{
  double NormError = MACHINE_PRECISION;
  int TargetSpaceNbrOrbitals = this->GetNbrOrbitals() / periodicity;
  ComplexVector** TmpVectors = new ComplexVector*[this->GetNbrParticles() + 1];
  for (int i = 0; i <= this->GetNbrParticles(); ++i)
    {
      TmpVectors[i] = new ComplexVector[TargetSpaceNbrOrbitals];
    }
  this->SymmetrizeSingleStatePeriodicSubsetOrbitalCore(inputVector, TmpVectors, firstOrbitalIndex, periodicity, phase, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i <= this->GetNbrParticles(); ++i)
    {
      for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
	{
	  if (TmpVectors[i][j].GetVectorDimension() != 0)	
	    {
	      ++NbrGeneratedSectors;
	    }     
	}
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new ComplexVector[NbrGeneratedSectors];
  kySectors = new int[NbrGeneratedSectors];
  nbrParticlesSectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i <= this->GetNbrParticles(); ++i)
    {
      for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
	{
	  if (TmpVectors[i][j].GetVectorDimension() != 0)	
	    {
	      double TmpNorm = TmpVectors[i][j].Norm();
	      if (TmpNorm > NormError)
		{
		  TmpVectors[i][j] /= TmpNorm; 
		  symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i][j];
		  kySectors[NbrGeneratedSectors] = j;
		  nbrParticlesSectors[NbrGeneratedSectors] = i;
		  ++NbrGeneratedSectors;	  
		}     
	    }
	}
    }  
  return NbrGeneratedSectors;
}

// symmetrized a product of two uncoupled states 
//
// outputVector = reference on the vector which will contain the symmetrozed state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

void ParticleOnTorus::SymmetrizeU1U1StateCore (RealVector& symmetrizedVector, RealVector& leftVector, RealVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
}

// symmetrized a product of two uncoupled states 
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

void ParticleOnTorus::SymmetrizeU1U1StateCore (ComplexVector& symmetrizedVector, ComplexVector& leftVector, ComplexVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
}

// symmetrize a product of several uncoupled states 
//
// outputState = reference on the output state
// inputStates = states which will be symmetrized
// inputSpaces = Hilbert spaces attached to each states
// nbrStates = number of states to symmetrize
// firstComponent = first component to symmetrize within the first Hilbert space of inputSpaces
// nbrComponents = number of components to symmetrize within the first Hilbert space of inputSpaces

void ParticleOnTorus::SymmetrizeU1U1StateCore (RealVector& outputState, RealVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, unsigned long firstComponent, unsigned long nbrComponents)
{
}
  
// symmetrized a product of several uncoupled states 
//
// outputState = reference on the output state
// inputStates = states which will be symmetrized
// inputSpaces = Hilbert spaces attached to each states
// nbrStates = number of states to symmetrize
// firstComponent = first component to symmetrize within the first Hilbert space of inputSpaces
// nbrComponents = number of components to symmetrize within the first Hilbert space of inputSpaces

void ParticleOnTorus::SymmetrizeU1U1StateCore (ComplexVector& outputState, ComplexVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, unsigned long firstComponent, unsigned long nbrComponents)
{
}
  

// symmetrize a vector by grouping neighbouring orbitals, core part
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
// first component = index of the first vector component 
// last component = index of the last component

void ParticleOnTorus::SymmetrizeSingleStateGroupingNeighbouringOrbitalsCore (ComplexVector& inputVector, ComplexVector* symmetrizedVectors, int nbrOrbitals, 
									     unsigned long firstComponent, unsigned long nbrComponents)
{
}

// symmetrize a vector by grouping distant and equally separated orbitals, core part
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
// first component = index of the first vector component 
// last component = index of the last component

void ParticleOnTorus::SymmetrizeSingleStateGroupingDistantOrbitalsCore (ComplexVector& inputVector, ComplexVector* symmetrizedVectors, 
									int nbrOrbitals, unsigned long firstComponent, unsigned long nbrComponents, bool twistedTorus)
{
}

// symmetrize a vector by keeping only a subset of equally separated orbitals
//
// inputVector = reference on the vector to symmetrize
// firstOrbitalIndex = index of the first orbital to keep
// symmetrizedVectors = array on the symmetrize states ranging from the smallest Ky to the largest Ky
// periodicity = momentum periodicity (should be a multiple of the number of orbitals)
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = symmetrized state

void ParticleOnTorus::SymmetrizeSingleStatePeriodicSubsetOrbitalCore (ComplexVector& inputVector, ComplexVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, 
								      unsigned long firstComponent, unsigned long nbrComponents)
{
}

// symmetrize a vector by keeping only a subset of equally separated orbitals
//
// inputVector = reference on the vector to symmetrize
// firstOrbitalIndex = index of the first orbital to keep
// symmetrizedVectors = array on the symmetrize states ranging from the smallest Ky to the largest Ky
// periodicity = momentum periodicity (should be a multiple of the number of orbitals)
// phase = an optional phase (in pi units) that can be added for each kept and occupied orbital
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = symmetrized state

void ParticleOnTorus::SymmetrizeSingleStatePeriodicSubsetOrbitalCore (ComplexVector& inputVector, ComplexVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, double phase, 
								      unsigned long firstComponent, unsigned long nbrComponents)
{
}

// create a state from its MPS description
//
// bMatrices = array that gives the B matrices 
// twistMatrix = reference on the twist matrix to insert in the trace
// state = reference to vector that will contain the state description
// mPSSumIndices = diagonal indices that have to be kept in the trace
// nbrMPSSumIndices = number of diagonal indices that have to be kept in the trace
// memory = amount of memory that can be use to precompute matrix multiplications  
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void ParticleOnTorus::CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, SparseRealMatrix& twistMatrix, RealVector& state, 
						     int* mPSSumIndices, int nbrMPSSumIndices,
						     long memory, long initialIndex, long nbrComponents)
{
}

// create a state from its MPS description
//
// bMatrices = array that gives the B matrices 
// twistMatrix = reference on the twist matrix to insert in the trace
// state = reference to vector that will contain the state description
// mPSSumIndices = diagonal indices that have to be kept in the trace
// nbrMPSSumIndices = number of diagonal indices that have to be kept in the trace
// memory = amount of memory that can be use to precompute matrix multiplications  
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void ParticleOnTorus::CreateStateFromMPSDescription (SparseComplexMatrix* bMatrices, SparseRealMatrix& twistMatrix, ComplexVector& state, 
						     int* mPSSumIndices, int nbrMPSSumIndices,
						     long memory, long initialIndex, long nbrComponents)
{
}

// transform a state expressed on a torus with a given angle to a state expressed on the same trous but a different angle
//
// inputVector = reference on the input vector
// inputAngle = angle (in radian) between the two vectors that span the torus on which the input state is defined
// inputAspectRatio = length ratio of the two vectors that span the torus on which the input state is defined
// outputAngle = angle (in radian) between the two vectors that span the torus on which the output state is defined
// outputAspectRatio = length ratio of the two vectors that span the torus on which the output state is defined
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = transformed state

ComplexVector ParticleOnTorus::ChangeTorusAngle (ComplexVector& inputVector, double inputAngle, double inputAspectRatio, double outputAngle, double outputAspectRatio,
						 unsigned long firstComponent, unsigned long nbrComponents)
{
  return ComplexVector();
}

