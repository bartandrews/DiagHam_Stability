////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//          class of hard-core boson with a single quantum number             //
//                                                                            //
//                        last modification : 11/02/2008                      //
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
#include "HilbertSpace/FermionOnLatticeGenericMomentumSpace.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/StringTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "QuantumNumber/NumberParticleQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"


#include <bitset>
#include <iostream>
#include <cstdlib>

using std::bitset;
using std::cout;
using std::endl;

// default constructor
FermionOnLatticeGenericMomentumSpace::FermionOnLatticeGenericMomentumSpace()
{
  this->NbrFermions=0;
  this->HilbertSpaceDimension=0;
}


// basic constructor -> yields a square lattice in Landau gauge
// 
// nbrFermions = number of bosons
// nx = number of momenta in simulation cell in x-direction
// ny = number of momenta in simulation cell in y-direction
// kx = total momentum in x-direction
// ky = total momentum in y-direction
// nbrBands = number of single particle bands
// memory = memory that can be allocated for precalculations
// verbose = flag indicating if any output is wanted
FermionOnLatticeGenericMomentumSpace::FermionOnLatticeGenericMomentumSpace (int nbrFermions, int nx, int ny, int kx, int ky, int nbrBands, unsigned long memory, bool verbose)
{
  this->NbrFermions = nbrFermions;
  this->Nx = nx;
  this->Ny = ny;
  this->Kx=kx;
  this->Ky=ky;
  this->NbrBands = nbrBands;
  this->NbrKPoints = Nx*Ny;
  this->NbrStates = Nx*Ny*NbrBands;

#ifdef __64_BITS__  
  if (this->NbrStates>64)  
#else
  if (this->NbrStates>32)    
#endif
    {
      cout<<"FermionOnLatticeGenericMomentumSpace: Cannot represent the "<<NbrStates<<" states requested in a single word"<<endl;
      exit(1);
    }

  this->Flag.Initialize();

  //this->LargeHilbertSpaceDimension = EvaluateHilbertSpaceDimension(this->NbrFermions, this->Nx-1, this->Ny-1, 0, 0);
  this->LargeHilbertSpaceDimension = EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrStates-1, 0, 0);
  cout << "Dimension count="<<LargeHilbertSpaceDimension<<endl;
  this->StateDescription=new unsigned long[this->LargeHilbertSpaceDimension];
  this->StateHighestBit=new int[this->LargeHilbertSpaceDimension];  
  long TmpDim = this->GenerateStates(this->NbrFermions, this->NbrStates-1, 0, 0, 0);
  if (TmpDim!=this->LargeHilbertSpaceDimension)
    {
      cout << "Dimension mismatch in FermionOnLatticeGenericMomentumSpace!"<<endl;
      exit(1);
    }
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->TargetSpace=this;
  
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);

//   for (int i=0; i<HilbertSpaceDimension; ++i)
//     {
//       this->PrintState(cout,i);
//       cout << endl;
//     }

  this->Flag.Initialize();
#ifdef __DEBUG__
  // for (int i=0; i<HilbertSpaceDimension; ++i)
//     {
//       this->PrintState(cout,i);
//       cout << endl;
//     }
  unsigned long UsedMemory = 0;  
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  // review here!
  UsedMemory += this->NbrStates * sizeof(int);
  UsedMemory += this->NbrStates * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
  UsedMemory +=  (1 << this->MaximumSignLookUp) * sizeof(double);
  if (verbose)
    {
      cout << "memory requested for Hilbert space = ";
      PrintMemorySize(cout,UsedMemory)<<endl;
    }
#endif
  if (verbose)
    cout << "Hilbert space dimension = " << this->HilbertSpaceDimension<<endl;
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
FermionOnLatticeGenericMomentumSpace::FermionOnLatticeGenericMomentumSpace(const FermionOnLatticeGenericMomentumSpace& fermions)
{
  this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->Nx = fermions.Nx;
  this->Ny = fermions.Ny;
  this->Kx = fermions.Kx;
  this->Ky = fermions.Ky;
  this->NbrBands = fermions.NbrBands;
  this->NbrStates = fermions.NbrStates;
  this->NbrKPoints = fermions.NbrKPoints;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}


// virtual destructor
//

FermionOnLatticeGenericMomentumSpace::~FermionOnLatticeGenericMomentumSpace ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrStates; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// assignment (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space
FermionOnLatticeGenericMomentumSpace& FermionOnLatticeGenericMomentumSpace::operator = (const FermionOnLatticeGenericMomentumSpace& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrStates; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
  if (this->TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->Nx = fermions.Nx;
  this->Ny = fermions.Ny;
  this->Kx = fermions.Kx;
  this->Ky = fermions.Ky;
  this->NbrBands = fermions.NbrBands;
  this->NbrStates = fermions.NbrStates;
  this->NbrKPoints = fermions.NbrKPoints;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  return *this;
}


// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space
AbstractHilbertSpace* FermionOnLatticeGenericMomentumSpace::Clone()
{
  return new FermionOnLatticeGenericMomentumSpace(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnLatticeGenericMomentumSpace::SetTargetSpace(ParticleOnLattice* targetSpace)
{
  this->TargetSpace=(FermionOnLatticeGenericMomentumSpace*)targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int FermionOnLatticeGenericMomentumSpace::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->GetHilbertSpaceDimension();
}

// get the particle statistic 
//
// return value = particle statistic
int FermionOnLatticeGenericMomentumSpace::GetParticleStatistic()
{
  return AbstractQHEParticle::FermionicStatistic;
}

// get the quantization axis 
//
// return value = particle statistic
char FermionOnLatticeGenericMomentumSpace::GetLandauGaugeAxis()
{
  return '\0';
}


// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id

int FermionOnLatticeGenericMomentumSpace::GetHilbertSpaceAdditionalSymmetry()
{
  return ParticleOnLattice::XTranslations|ParticleOnLattice::YTranslations;
}


// check whether HilbertSpace implements ordering of operators
//
bool FermionOnLatticeGenericMomentumSpace::HaveOrder ()
{
  return true;
}


// check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
int FermionOnLatticeGenericMomentumSpace::CheckOrder (int* m, int* n, int nbrIndices)
{
  unsigned long CreationValue=0x0ul;
  unsigned long AnnihilationValue=0x0ul;
  for (int i=0; i<nbrIndices; ++i)
    {
      CreationValue |= 0x1ul << m[i];
      AnnihilationValue |= 0x1ul << n[i];
    }
  if (CreationValue > AnnihilationValue)
    return 1;
  else if (CreationValue < AnnihilationValue)
    return -1;
  else return 0;
}


// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number
List<AbstractQuantumNumber*> FermionOnLatticeGenericMomentumSpace::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new NumberParticleQuantumNumber(this->NbrFermions);
  TmpList += new PeriodicMomentumQuantumNumber (this->Kx, this->Nx);
  TmpList += new PeriodicMomentumQuantumNumber (this->Ky, this->Ny);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number
AbstractQuantumNumber* FermionOnLatticeGenericMomentumSpace::GetQuantumNumber (int index)
{
  if (index==0)
    return new NumberParticleQuantumNumber(this->NbrFermions);
  else if (index==1)
    return new PeriodicMomentumQuantumNumber (this->Kx, this->Nx);
  else if (index==2)
    return new PeriodicMomentumQuantumNumber (this->Ky, this->Ny);
  return NULL;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace
AbstractHilbertSpace* FermionOnLatticeGenericMomentumSpace::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}


// get the number of sites
//
// return value = number of sites
int FermionOnLatticeGenericMomentumSpace::GetNbrSites()
{
  return this->Nx*this->Ny;
}

// it is possible to change the flux through the simulation cell
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
void FermionOnLatticeGenericMomentumSpace::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  cout << "NbrFluxQuanta not defined in FermionOnLatticeGenericMomentumSpace"<<endl;
}

// change flux through cell and periodic boundary conditions
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
void FermionOnLatticeGenericMomentumSpace::SetNbrFluxQuanta(int nbrFluxQuanta, double solenoidX, double solenoidY)
{
  cout << "NbrFluxQuanta not defined in FermionOnLatticeGenericMomentumSpace"<<endl;
}

// request solenoid fluxes
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
//
void FermionOnLatticeGenericMomentumSpace::GetSolenoidFluxes(double &solenoidX, double &solenoidY)
{
  solenoidX=0.0;
  solenoidY=0.0;
}


// obtain the current setting of the flux piercing this lattice
int FermionOnLatticeGenericMomentumSpace::GetNbrFluxQuanta()
{
  return 0;
}


// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
unsigned long FermionOnLatticeGenericMomentumSpace::Ad (unsigned long state, int q, double& coefficient)
{  
  if ((state & (((unsigned long) (0x1)) << q))!= 0)
    {
      return 0x0l;
    }
  coefficient=1.0;
  state |= (((unsigned long) (0x1)) << q);
  return state;
}


// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnLatticeGenericMomentumSpace::AdAdAA (int index, int m1, int m2, int n1, int n2, double &coefficient)
{
  coefficient = 1.0;
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((State & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  // perform annihilation operators
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif  
  TmpState &= ~(((unsigned long) (0x1)) << n2);
  if (NewHighestBit == n2)
    while ((TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n1);
  if (NewHighestBit == n1)
    while ((TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  // create particle at m2
  if ((TmpState & (((unsigned long) (0x1)) << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewHighestBit)
    {
      NewHighestBit = m2;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m2);
  // create particle at m1
  if ((TmpState & (((unsigned long) (0x1)) << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    {
      NewHighestBit = m1;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m1);
  return this->TargetSpace->FindStateIndex(TmpState, NewHighestBit);
}


// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double FermionOnLatticeGenericMomentumSpace::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];

  if (((ProdATemporaryState & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((ProdATemporaryState & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2))
    return 0.0;

  this->ProdAHighestBit = this->StateHighestBit[index];

  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n2);

  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n1);

  while ((this->ProdATemporaryState >> this->ProdAHighestBit) == 0)
    --this->ProdAHighestBit;
  return Coefficient;
}


// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on a multiplicative, trivially one here (unreferenced)
// return value = index of the destination state 

int FermionOnLatticeGenericMomentumSpace::AdAd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if ((TmpState & (((unsigned long) (0x1)) << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = this->ProdAHighestBit;
  if (m2 > NewHighestBit)
    NewHighestBit = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m2);
  if ((TmpState & (((unsigned long) (0x1)) << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    NewHighestBit = m1;
  else
    {      
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m1);
  return this->TargetSpace->FindStateIndex(TmpState, NewHighestBit);
}


// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored (always 1.0)
// return value = index of the destination state 

int FermionOnLatticeGenericMomentumSpace::AdA (int index, int m, int n, double &coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];  
  if (m!=n)
    {
      if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      int NewHighestBit = StateHighestBit;
      unsigned long TmpState = State;
      // perform annihilation operators
      TmpState &= ~(((unsigned long) (0x1)) << n);
      coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
      coefficient *= this->SignLookUpTable[(TmpState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
      if (NewHighestBit == n)
	while ((NewHighestBit > 0)&&((TmpState >> NewHighestBit) == 0))
	  --NewHighestBit;
      // create particle at m
      if ((TmpState & (((unsigned long) (0x1)) << m))!= 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (m > NewHighestBit)
	{
	  NewHighestBit = m;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
	}
      TmpState |= (((unsigned long) (0x1)) << m);
      return this->TargetSpace->FindStateIndex(TmpState, NewHighestBit);
    }
  else
    {
      if ((m > StateHighestBit) || ((State & (((unsigned long) (0x1)) << m)) == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      else
	{
	  coefficient = 1.0;
	  return index;
	}
    }
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double FermionOnLatticeGenericMomentumSpace::AdA (int index, int m)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  if ((m > StateHighestBit) || ((State & (((unsigned long) (0x1)) << m)) == 0))
    return 0.0;
  else
    return 1.0;
}

// apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
// index = index of the state on which the operator has to be applied
// nbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
// qValues = array of quantum numbers where an interaction is present
// interactionPerQ = coefficient U_q of the interaction
//
// no diagonal interaction present, derelict function from inheritance
//
double FermionOnLatticeGenericMomentumSpace::AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  return 0.0;
}



// obtain a list of quantum numbers in state
// index = index of many-body state to be considered
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
// normalization = indicating the multiplicity of the state for bosonic spaces
void FermionOnLatticeGenericMomentumSpace::ListQuantumNumbers(int index, int *quantumNumbers, double &normalization)
{
  normalization=1.0;
  this->ListQuantumNumbers(index, quantumNumbers);
}

// obtain a list of quantum numbers in state
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
void FermionOnLatticeGenericMomentumSpace::ListQuantumNumbers(int index, int *quantumNumbers)
{
  unsigned long State = this->StateDescription[index];
  int HighestBit = this->StateHighestBit[index];
  int NbrQ=0;
  for (int q=0; q<=HighestBit; ++q)
    if (State&(0x1u<<q))
      quantumNumbers[NbrQ++]=q;
}

// translate a state by a multiple of the lattice vectors
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// translationPhase = returns phase inccurred by translation
// return value = index of translated state
int FermionOnLatticeGenericMomentumSpace::TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase)
{
  translationPhase = 0.0;
  return index;
}

// find whether there is a translation vector from state i to state f
// i = index of initial state
// f = index of final state
// shiftX = length of translation in x-direction
// shiftY = length of translation in y-direction
// return value = final state can be reached by translation
bool FermionOnLatticeGenericMomentumSpace::IsTranslation(int i, int f, int &shiftX, int &shiftY)
{
  shiftX=0;
  shiftY=0;
  if (i==f)
    return true;
  else
    return false;
}


  // evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex FermionOnLatticeGenericMomentumSpace::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// return value = wave function evaluated at the given location

Complex FermionOnLatticeGenericMomentumSpace::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, 
						     this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location
Complex FermionOnLatticeGenericMomentumSpace::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
						int firstComponent, int nbrComponent)
{
  return Complex(0.0, 0.0);
}


// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnLatticeGenericMomentumSpace::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, 
								 int nextCoordinates, int firstComponent, 
								 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnLatticeGenericMomentumSpace::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnLatticeGenericMomentumSpace::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrStates; ++i)
    Str << ((TmpState >> i) & ((unsigned long) 0x1)) << " ";
  Str << " position = " << this->FindStateIndex(TmpState, this->StateHighestBit[state]) ;// << ", hb="<<this->StateHighestBit[state];
  if (state !=  this->FindStateIndex(TmpState, this->StateHighestBit[state]))
    Str << " error! ";
  return Str;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnLatticeGenericMomentumSpace::FindStateIndex(unsigned long stateDescription, int highestBit)
{
  // bitset <32> b = stateDescription;
//   cout << "in FindStateIndex: desc=" << b<<", highest bit=" << highestBit << endl;
  long PosMax = stateDescription >> this->LookUpTableShift[highestBit];
  long PosMin = this->LookUpTable[highestBit][PosMax];
  PosMax = this->LookUpTable[highestBit][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMin;
}

// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found
int FermionOnLatticeGenericMomentumSpace::CarefulFindStateIndex(unsigned long stateDescription, int highestBit)
{
//   bitset<20> b=stateDescription;
//   cout << "Searching " << b <<", bitcount="<<bitcount(stateDescription);
  if (bitcount(stateDescription)!=this->NbrFermions)
    {
      return this->HilbertSpaceDimension;
    }
  if (highestBit<0)
    {
      highestBit = getHighestBit(stateDescription)-1;
    }
//   cout << ", highestBit="<<highestBit;
  if (highestBit >= this->NbrStates)
    {
      return this->HilbertSpaceDimension;
    }
//   cout <<", Found index="<<this->FindStateIndex(stateDescription, highestBit)<<endl;
//   this->PrintState(cout, this->FindStateIndex(stateDescription, highestBit));
  return this->FindStateIndex(stateDescription, highestBit);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentQ = current largest combined quantum number
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnLatticeGenericMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentQ, int currentTotalKx, int currentTotalKy)
{
  if ((nbrFermions == 0) || (currentQ < 0))
    return 0l;
  int CurrentKx, CurrentKy, CurrentBand;
  DecodeQuantumNumber(currentQ, CurrentKx, CurrentKy, CurrentBand);
  long Count = 0;

  if (nbrFermions == 1)
    {
      int ky = this->Ky - (currentTotalKy % this->Ny);
      if (ky < 0)
	ky += this->Ny;
      int kx = this->Kx - (currentTotalKx % this->Nx);
      if (kx < 0)
	kx += this->Nx;
      int Q;
      for (int b=0; b <= CurrentBand; ++b)
	{
	  Q = FastEncodeQuantumNumber(kx, ky, b);
	  if (Q<=currentQ)
	    ++Count;
	}
      return Count;
    }
  
  int ReducedCurrentQ = currentQ - 1;
  Count = this->EvaluateHilbertSpaceDimension(nbrFermions-1, ReducedCurrentQ, currentTotalKx+CurrentKx, currentTotalKy+CurrentKy);
  return Count+EvaluateHilbertSpaceDimension(nbrFermions, ReducedCurrentQ, currentTotalKx, currentTotalKy);
}


// generate states with kx and ky symmetries
//
// nbrFermions = number of fermions
// currentQ = current consolidated quantum number
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// pos = position where to start filling states
// return value = Hilbert space dimension

long FermionOnLatticeGenericMomentumSpace::GenerateStates(int nbrFermions, int currentQ, int currentTotalKx, int currentTotalKy, long pos)
{
  if ((nbrFermions == 0) || (currentQ < 0))
    return pos;
  int CurrentKx, CurrentKy, CurrentBand;
  DecodeQuantumNumber(currentQ, CurrentKx, CurrentKy, CurrentBand);
  
  if (nbrFermions == 1)
    {
      int ky = this->Ky - (currentTotalKy % this->Ny);
      if (ky < 0)
	ky += this->Ny;
      int kx = this->Kx - (currentTotalKx % this->Nx);
      if (kx < 0)
	kx += this->Nx;
      int Q;
      for (int b=CurrentBand; b >= 0; --b)
	{
	  Q = FastEncodeQuantumNumber(kx, ky, b);
	  if (Q<=currentQ)
	    {
	      this->StateDescription[pos] = 0x1ul << Q;
	      this->StateHighestBit[pos] = Q;
	      ++pos;
	    }
	}
      return pos;
    }
  
  int ReducedCurrentQ = currentQ - 1;
  int TmpPos = this->GenerateStates(nbrFermions-1, ReducedCurrentQ, currentTotalKx+CurrentKx, currentTotalKy+CurrentKy,pos);
  unsigned long Mask = ((unsigned long) 1) << currentQ;
  for (int i = pos; i < TmpPos; i++)
    {
      this->StateDescription[i] |= Mask;
      this->StateHighestBit[i] = currentQ;
    }
  return GenerateStates(nbrFermions, ReducedCurrentQ, currentTotalKx, currentTotalKy, TmpPos);
}


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnLatticeGenericMomentumSpace::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrStates);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrStates)
    this->MaximumLookUpShift = this->NbrStates;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrStates];
  this->LookUpTableShift = new int [this->NbrStates];
  for (int i = 0; i < this->NbrStates; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentHighestBit = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
  if (CurrentHighestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentHighestBit] = 0;
  else
    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentHighestBit];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentHighestBit != this->StateHighestBit[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  /*	  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
	    cout << TmpLookUpTable[j] << " ";
	    cout << endl << "-------------------------------------------" << endl;*/
 	  CurrentHighestBit = this->StateHighestBit[i];
	  TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
	  if (CurrentHighestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentHighestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentHighestBit];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;
  /*  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
    cout << TmpLookUpTable[j] << " ";
    cout << endl << "-------------------------------------------" << endl;*/
  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	this->SignLookUpTable[j] = -1.0;
      else
	this->SignLookUpTable[j] = 1.0;
    }
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif
}

