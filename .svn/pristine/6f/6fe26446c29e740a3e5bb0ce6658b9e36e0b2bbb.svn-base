////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of hardcore boson on lattice in real space        //
//                                                                            //
//                        last modification : 10/09/2014                      //
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
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

BosonOnLatticeGutzwillerProjectionRealSpace::BosonOnLatticeGutzwillerProjectionRealSpace ()
{
  this->NbrBosons = 0;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrSite = 0;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateLzMax = 0;  
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// memory = amount of memory granted for precalculations

BosonOnLatticeGutzwillerProjectionRealSpace::BosonOnLatticeGutzwillerProjectionRealSpace (int nbrBosons, int nbrSite, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSite - 1, 0l);
      this->StateLzMax = new int [this->LargeHilbertSpaceDimension];  
      int CurrentLzMax = this->NbrLzValue;
      while (((this->StateDescription[0] >> CurrentLzMax) & 0x1ul) == 0x0ul)
	--CurrentLzMax;
      this->StateLzMax[0] = CurrentLzMax;
      for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  while (((this->StateDescription[i] >> CurrentLzMax) & 0x1ul) == 0x0ul)
	    --CurrentLzMax;
	  this->StateLzMax[i] = CurrentLzMax;
 	}
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space, " << TmpLargeHilbertSpaceDimension << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
	}
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

BosonOnLatticeGutzwillerProjectionRealSpace::BosonOnLatticeGutzwillerProjectionRealSpace(const BosonOnLatticeGutzwillerProjectionRealSpace & bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->NbrSite = bosons.NbrSite;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  this->SignLookUpTable = bosons.SignLookUpTable;
  this->SignLookUpTableMask = bosons.SignLookUpTableMask;
  this->MaximumSignLookUp = bosons.MaximumSignLookUp;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnLatticeGutzwillerProjectionRealSpace::~BosonOnLatticeGutzwillerProjectionRealSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeGutzwillerProjectionRealSpace& BosonOnLatticeGutzwillerProjectionRealSpace::operator = (const BosonOnLatticeGutzwillerProjectionRealSpace& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
    }
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrSite = bosons.NbrSite;
  this->NbrLzValue = bosons.NbrLzValue;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeGutzwillerProjectionRealSpace::Clone()
{
  return new BosonOnLatticeGutzwillerProjectionRealSpace(*this);
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnLatticeGutzwillerProjectionRealSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  if (nbrParticleSector == this->NbrBosons)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  BosonOnLatticeGutzwillerProjectionRealSpace SubsytemSpace (nbrParticleSector, this->NbrSite);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnLatticeGutzwillerProjectionRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;

  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// apply a_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator
// return value =  multiplicative factor 

double BosonOnLatticeGutzwillerProjectionRealSpace::A (int index, int n)
{
  this->ProdATemporaryState = this->StateDescription[index];
  int StateLzMax = this->StateLzMax[index];

  if ((n > StateLzMax) || ((this->ProdATemporaryState & (0x1ul << n)) == 0x0ul))
    return 0.0;

  this->ProdALzMax = this->StateLzMax[index];

  double Coefficient = 1.0;
  
  this->ProdATemporaryState &= ~(0x1ul << n);

  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdALzMax = 0;
      return Coefficient;      
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0x0ul)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a^+_n1  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
//
// index = index of the state on which the operator has to be applied
// m = index for annihilation operator
// return value =  multiplicative factor 

double BosonOnLatticeGutzwillerProjectionRealSpace::Ad (int index, int m)
{
  this->ProdATemporaryState = this->StateDescription[index];
  int StateLzMax = this->StateLzMax[index];

  if ((this->ProdATemporaryState & (0x1ul << m)) != 0x0ul)
    return 0.0;

  double Coefficient = 1.0;
   
  this->ProdATemporaryState |= (0x1ul << m);
  this->ProdALzMax = this->StateLzMax[index];
  if (m > this->ProdALzMax)
    this->ProdALzMax = m;
  return Coefficient;
}

// apply a^+_m operator to the state produced using the A or Ad method (without destroying it)
//
// m = index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnLatticeGutzwillerProjectionRealSpace::Ad (int m, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
    
  if ((TmpState & (0x1ul << m))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = this->ProdALzMax;
  if (m > NewLzMax)
    NewLzMax = m;

  TmpState |= (0x1ul << m);
  return this->FindStateIndex(TmpState, NewLzMax);
}


// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnLatticeGutzwillerProjectionRealSpace::AdA (int index, int m, int n, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  unsigned long State = this->StateDescription[index];
  if ((n > StateLzMax) || ((State & (0x1ul << n)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  unsigned long TmpState = State;
  coefficient = 1.0;
  TmpState &= ~(0x1ul << n);
  if ((TmpState != 0x0ul))
    {
      while ((TmpState >> NewLzMax) == 0)
	--NewLzMax;
    }
  else
    NewLzMax = 0;
  if ((TmpState & (0x1ul << m))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLzMax)
    {
      NewLzMax = m;
    }

  TmpState |= (0x1ul << m);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

long BosonOnLatticeGutzwillerProjectionRealSpace::AdA (long index, int m, int n, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  unsigned long State = this->StateDescription[index];
  if ((n > StateLzMax) || ((State & (0x1ul << n)) == 0x0ul))
    {
      coefficient = 0.0;
      return this->LargeHilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  unsigned long TmpState = State;
  coefficient = 1.0;
  TmpState &= ~(0x1ul << n);
  if (TmpState != 0x0ul)
    {
      while ((TmpState >> NewLzMax) == 0)
	--NewLzMax;
    }
  else
    NewLzMax = 0;
  if ((TmpState & (0x1ul << m))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->LargeHilbertSpaceDimension;
    }
  if (m > NewLzMax)
    {
      NewLzMax = m;
    }

  TmpState |= (0x1ul << m);
  return this->FindStateIndex(TmpState, NewLzMax);
}


// evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index 
// nbrKeptOrbitals = array of orbitals that have to be kept
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem

ComplexMatrix BosonOnLatticeGutzwillerProjectionRealSpace::EvaluatePartialEntanglementMatrix (int nbrParticleSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  if ((nbrParticleSector >  nbrKeptOrbitals) || 
      (ComplementaryNbrParticles > (this->NbrSite - nbrKeptOrbitals)))
    {
      ComplexMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
  if (nbrKeptOrbitals == 0)
    {
      if (nbrParticleSector == 0)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, this->HilbertSpaceDimension, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      TmpEntanglementMatrix[i][0] = groundState[i];
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  if (nbrKeptOrbitals == this->NbrSite)
    {
      if (nbrParticleSector == this->NbrBosons)
	{
	  ComplexMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      TmpEntanglementMatrix[0][i] = groundState[i];
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  this->KeptOrbitals = new int [nbrKeptOrbitals];
  for (int i = 0 ; i < nbrKeptOrbitals; ++i) 
    this->KeptOrbitals[i] = keptOrbitals[i];
  BosonOnLatticeGutzwillerProjectionRealSpace SubsytemSpace (nbrParticleSector, nbrKeptOrbitals);
  BosonOnLatticeGutzwillerProjectionRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite - nbrKeptOrbitals);
  ComplexMatrix TmpEntanglementMatrix (SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySpace.HilbertSpaceDimension, true);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;

  long TmpEntanglementMatrixZero = this->EvaluatePartialEntanglementMatrixCore(0, ComplementarySpace.HilbertSpaceDimension, &ComplementarySpace, &SubsytemSpace, groundState, &TmpEntanglementMatrix);
//   FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpEntanglementMatrix);
//   Operation.ApplyOperation(architecture);
//   if (Operation.GetNbrNonZeroMatrixElements() > 0)	
  if (TmpEntanglementMatrixZero > 0)
     return TmpEntanglementMatrix;
   else
     {
       ComplexMatrix TmpEntanglementMatrixZero;
       return TmpEntanglementMatrixZero;
     }    
}
  


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnLatticeGutzwillerProjectionRealSpace::PrintState (ostream& Str, int state)
{
  unsigned long * TemporaryState = new unsigned long[this->NbrBosons];
  unsigned long TmpState = this->StateDescription[state];
 this->ConvertToMonomial(TmpState,TemporaryState);  
  for (int i = 0; i < this->NbrBosons; ++i)
    cout <<TemporaryState[i]<<" ";
  cout <<endl;

  for (int i = 0; i < this->NbrSite; ++i)
    Str << ((TmpState >> i) & 0x1ul) << " ";
  cout <<TmpState<<endl;
delete TemporaryState;
  return Str;
}



void BosonOnLatticeGutzwillerProjectionRealSpace::GetCompositeFermionWavefunction(ComplexVector & trialState, ComplexMatrix & jastrowEigenVecs,ComplexMatrix & cFEigenVecs)
{
#ifdef __LAPACK__
  ComplexLapackDeterminant SlaterCF(NbrBosons);
  ComplexLapackDeterminant SlaterJastrow(NbrBosons);
#else
  ComplexMatrix SlaterCF(NbrBosons, NbrBosons);
  ComplexMatrix SlaterJastrow(NbrBosons, NbrBosons);
#endif
  
  unsigned long * TemporaryState = new unsigned long [this->NbrBosons];
  
 for(int i = 0; i < this->HilbertSpaceDimension ; i++)
 {
   this->ConvertToMonomial(this->StateDescription[i],TemporaryState);
   
   for (int p = 0; p < NbrBosons; ++p)
     {
       for (int q = 0; q < NbrBosons; ++q)
	 {
	   SlaterCF.SetMatrixElement(p,q,cFEigenVecs[p][(int) TemporaryState[q]]);
	   SlaterJastrow.SetMatrixElement(p,q,jastrowEigenVecs[p][(int) TemporaryState[q]]);
	   
	   //   SlaterCF.SetMatrixElement(p,q,Conj(cFEigenVecs[p][(int) TemporaryState[q]]));
	   //   SlaterJastrow.SetMatrixElement(p,q,Conj(jastrowEigenVecs[p][(int) TemporaryState[q]]));
	 }	      
     }
   trialState[i] +=  SlaterCF.Determinant() * SlaterJastrow.Determinant();
 }
 delete [] TemporaryState;
}


// apply a_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator
// return value =  multiplicative factor 

double BosonOnLatticeGutzwillerProjectionRealSpace::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  int StateLzMax = this->StateLzMax[index];

  if ((n1 > StateLzMax) || (n2 > StateLzMax) || ((this->ProdATemporaryState & (0x1ul << n1)) == 0x0ul) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0x0ul))
    return 0.0;

  this->ProdALzMax = this->StateLzMax[index];

  double Coefficient = 1.0;
  
  this->ProdATemporaryState &= ~(0x1ul << n2);
  this->ProdATemporaryState &= ~(0x1ul << n1);

  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdALzMax = 0;
      return Coefficient;      
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0x0ul)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a^+_m operator to the state produced using the A or Ad method (without destroying it)
//
// m1 = index for creation operator
// m2 = index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnLatticeGutzwillerProjectionRealSpace::AdAd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
    
  if ( ((TmpState & (0x1ul << m1))!= 0x0ul) || ((TmpState & (0x1ul << m2))!= 0x0ul) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m1 > NewLzMax)
    NewLzMax = m1;
  
  if (m2 > NewLzMax)
    NewLzMax = m2;

  
  TmpState |= (0x1ul << m2);
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnLatticeGutzwillerProjectionRealSpace::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  
												     ParticleOnSphere* destinationHilbertSpace,
												     ComplexVector& groundState,  HermitianMatrix* densityMatrix)
{
  BosonOnLatticeGutzwillerProjectionRealSpace* TmpHilbertSpace =  (BosonOnLatticeGutzwillerProjectionRealSpace*) complementaryHilbertSpace;
  BosonOnLatticeGutzwillerProjectionRealSpace* TmpDestinationHilbertSpace =  (BosonOnLatticeGutzwillerProjectionRealSpace*) destinationHilbertSpace;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  BinomialCoefficients TmpBinomial (this->NbrBosons);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrBosons, TmpDestinationHilbertSpace->NbrBosons));
  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
 	      int TmpLzMax = this->LzMax;
	      unsigned long TmpState3 = TmpState | TmpState2;
	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  double Coefficient = TmpInvBinomial;
		  unsigned long Sign = 0x0ul;
		  int Pos2 = TmpDestinationHilbertSpace->LzMax;
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient[Pos] = Coefficient;
		  ++Pos;
		}
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
}

