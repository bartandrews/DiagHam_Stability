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
#include "HilbertSpace/BosonOnLatticeGenericMomentumSpace.h"
#include "HilbertSpace/HardCoreBoson.h"
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
BosonOnLatticeGenericMomentumSpace::BosonOnLatticeGenericMomentumSpace()
{
  this->NbrBosons=0;
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
BosonOnLatticeGenericMomentumSpace::BosonOnLatticeGenericMomentumSpace (int nbrBosons, int nx, int ny, int kx, int ky, int nbrBands, unsigned long memory, bool verbose)
{
  this->NbrBosons = nbrBosons;
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
      cout<<"BosonOnLatticeGenericMomentumSpace: Cannot represent the "<<NbrStates<<" states requested in a single word"<<endl;
      exit(1);
    }

  this->Flag.Initialize();

  //this->LargeHilbertSpaceDimension = EvaluateHilbertSpaceDimension(this->NbrBosons, this->Nx-1, this->Ny-1, 0, 0);
  this->LargeHilbertSpaceDimension = EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrStates-1, 0, 0);
  cout << "Dimension count="<<LargeHilbertSpaceDimension<<endl;
  this->StateDescription=new unsigned long[this->LargeHilbertSpaceDimension];
  this->StateHighestBit=new int[this->LargeHilbertSpaceDimension];  
  long TmpDim = this->GenerateStates();
  if (TmpDim!=this->LargeHilbertSpaceDimension)
    {
      cout << "Dimension mismatch in BosonOnLatticeGenericMomentumSpace!"<<endl;
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
BosonOnLatticeGenericMomentumSpace::BosonOnLatticeGenericMomentumSpace(const BosonOnLatticeGenericMomentumSpace& fermions)
{
  this->TargetSpace = this;
  this->NbrBosons = fermions.NbrBosons;
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

BosonOnLatticeGenericMomentumSpace::~BosonOnLatticeGenericMomentumSpace ()
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
BosonOnLatticeGenericMomentumSpace& BosonOnLatticeGenericMomentumSpace::operator = (const BosonOnLatticeGenericMomentumSpace& fermions)
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
  this->NbrBosons = fermions.NbrBosons;
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
AbstractHilbertSpace* BosonOnLatticeGenericMomentumSpace::Clone()
{
  return new BosonOnLatticeGenericMomentumSpace(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnLatticeGenericMomentumSpace::SetTargetSpace(ParticleOnLattice* targetSpace)
{
  this->TargetSpace=(BosonOnLatticeGenericMomentumSpace*)targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int BosonOnLatticeGenericMomentumSpace::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->GetHilbertSpaceDimension();
}

// get the particle statistic 
//
// return value = particle statistic
int BosonOnLatticeGenericMomentumSpace::GetParticleStatistic()
{
  return AbstractQHEParticle::BosonicStatistic;
}

// get the quantization axis 
//
// return value = particle statistic
char BosonOnLatticeGenericMomentumSpace::GetLandauGaugeAxis()
{
  return '\0';
}


// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id

int BosonOnLatticeGenericMomentumSpace::GetHilbertSpaceAdditionalSymmetry()
{
  return ParticleOnLattice::XTranslations|ParticleOnLattice::YTranslations;
}


// check whether HilbertSpace implements ordering of operators
//
bool BosonOnLatticeGenericMomentumSpace::HaveOrder ()
{
  return true;
}


// check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
int BosonOnLatticeGenericMomentumSpace::CheckOrder (int* m, int* n, int nbrIndices)
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
List<AbstractQuantumNumber*> BosonOnLatticeGenericMomentumSpace::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new NumberParticleQuantumNumber(this->NbrBosons);
  TmpList += new PeriodicMomentumQuantumNumber (this->Kx, this->Nx);
  TmpList += new PeriodicMomentumQuantumNumber (this->Ky, this->Ny);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number
AbstractQuantumNumber* BosonOnLatticeGenericMomentumSpace::GetQuantumNumber (int index)
{
  if (index==0)
    return new NumberParticleQuantumNumber(this->NbrBosons);
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
AbstractHilbertSpace* BosonOnLatticeGenericMomentumSpace::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}


// get the number of sites
//
// return value = number of sites
int BosonOnLatticeGenericMomentumSpace::GetNbrSites()
{
  return this->Nx*this->Ny;
}

// it is possible to change the flux through the simulation cell
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
void BosonOnLatticeGenericMomentumSpace::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  cout << "NbrFluxQuanta not defined in BosonOnLatticeGenericMomentumSpace"<<endl;
}

// change flux through cell and periodic boundary conditions
// Attention: this does require the Hamiltonian to be recalculated!!
// nbrFluxQuanta = number of quanta of flux piercing the simulation cell
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
void BosonOnLatticeGenericMomentumSpace::SetNbrFluxQuanta(int nbrFluxQuanta, double solenoidX, double solenoidY)
{
  cout << "NbrFluxQuanta not defined in BosonOnLatticeGenericMomentumSpace"<<endl;
}

// request solenoid fluxes
// solenoidX = new solenoid flux through torus in x-direction
// solenoidY = new solenoid flux through torus in y-direction
//
void BosonOnLatticeGenericMomentumSpace::GetSolenoidFluxes(double &solenoidX, double &solenoidY)
{
  solenoidX=0.0;
  solenoidY=0.0;
}


// obtain the current setting of the flux piercing this lattice
int BosonOnLatticeGenericMomentumSpace::GetNbrFluxQuanta()
{
  return 0;
}

// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
unsigned long BosonOnLatticeGenericMomentumSpace::Ad (unsigned long state, int q, double& coefficient)
{
  int StateHighestBit = getHighestBit(state)-1;
  this->FermionToBoson(state, StateHighestBit, this->TemporaryState, this->TemporaryStateHighestBit);
  if (q>this->TemporaryStateHighestBit)
    {
      for (int i = this->TemporaryStateHighestBit + 1; i < q; ++i)
	this->TemporaryState[i] = 0x0ul;
      this->TemporaryState[q] = 0x1ul;
      this->TemporaryStateHighestBit = q;
      coefficient = 1.0;
    }
  else
    {
      ++this->TemporaryState[q];
      coefficient = sqrt((double)this->TemporaryState[q]);
    }
  return this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit);
}


// apply a^+_q1 a^+_q2 a_r1 a_r2 operator to a given state (with q1+q2=n1+r2)
//
// index = index of the state on which the operator has to be applied
// q1 = first index for creation operator
// q2 = second index for creation operator
// r1 = first index for annihilation operator
// r2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnLatticeGenericMomentumSpace::AdAdAA (int index, int q1, int q2, int r1, int r2, double& coefficient)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  if ((r1 > this->TemporaryStateHighestBit) || (r2 > this->TemporaryStateHighestBit) || (this->TemporaryState[r1] == 0) || (this->TemporaryState[r2] == 0) || ((r1 == r2) && (this->TemporaryState[r1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = this->TemporaryState[r2];
  --this->TemporaryState[r2];
  if ((r2==TemporaryStateHighestBit)&&(this->TemporaryState[r2]==0)&&(TemporaryStateHighestBit>0))
    {
      --TemporaryStateHighestBit;
      while ((TemporaryStateHighestBit>0)&&(TemporaryState[TemporaryStateHighestBit] == 0))
	--TemporaryStateHighestBit;
    }
  coefficient *= this->TemporaryState[r1];
  --this->TemporaryState[r1];
  if ((r1==TemporaryStateHighestBit)&&(this->TemporaryState[r1]==0)&&(TemporaryStateHighestBit>0))
    {
      --TemporaryStateHighestBit;
      while ((TemporaryStateHighestBit>0)&&(TemporaryState[TemporaryStateHighestBit] == 0))
	--TemporaryStateHighestBit;
    }
  if (q2 > TemporaryStateHighestBit)
    {
      for (int i = this->TemporaryStateHighestBit + 1; i <= q2; ++i)
	this->TemporaryState[i] = 0ul;
      TemporaryStateHighestBit = q2;
    }
  ++this->TemporaryState[q2];
  coefficient *= this->TemporaryState[q2];
  if (q1 > TemporaryStateHighestBit)
    {
      for (int i = this->TemporaryStateHighestBit + 1; i <= q1; ++i)
	this->TemporaryState[i] = 0ul;
      TemporaryStateHighestBit = q1;
    }
  ++this->TemporaryState[q1];
  coefficient *= this->TemporaryState[q1];  
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
}

// apply a_r1 a_r2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// r1 = first index for annihilation operator
// r2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnLatticeGenericMomentumSpace::AA (int index, int r1, int r2)
{
#ifdef DEBUG_OUTPUT
  cout << "AA ("<<index <<" " <<r1 <<" " <<r2<<") ";
  this->PrintState(cout,index);
#endif
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->ProdATemporaryState, this->ProdATemporaryStateHighestBit);
  if ((r1 > this->ProdATemporaryStateHighestBit) || (r2 > this->ProdATemporaryStateHighestBit) || 
      (this->ProdATemporaryState[r1] == 0) || (this->ProdATemporaryState[r2] == 0) || ((r1 == r2) && (this->ProdATemporaryState[r1] == 1)))
    {
#ifdef DEBUG_OUTPUT
      if ((index==4)||(index==2))
	cout << "Vanishing on "<<index<<": "<<r1<<" "<<r2<<endl;
#endif
      return 0.0;
    }
  for (int i = this->ProdATemporaryStateHighestBit + 1; i < this->NbrStates; ++i)
    this->ProdATemporaryState[i] = 0ul;
  double Coefficient = this->ProdATemporaryState[r2];
  --this->ProdATemporaryState[r2];
  if ((r2==ProdATemporaryStateHighestBit)&&(ProdATemporaryStateHighestBit>0)&&(this->ProdATemporaryState[r2]==0))
    {
      --ProdATemporaryStateHighestBit;
      while ((ProdATemporaryStateHighestBit>0)&&(ProdATemporaryState[ProdATemporaryStateHighestBit] == 0))
	--ProdATemporaryStateHighestBit;
    }
  Coefficient *= this->ProdATemporaryState[r1];
  --this->ProdATemporaryState[r1];
  if ((r1==ProdATemporaryStateHighestBit)&&(ProdATemporaryStateHighestBit>0)&&(this->ProdATemporaryState[r1]==0))
    {
      --ProdATemporaryStateHighestBit;
      while ((ProdATemporaryStateHighestBit>0)&&(ProdATemporaryState[ProdATemporaryStateHighestBit] == 0))
	--ProdATemporaryStateHighestBit;
    }
#ifdef DEBUG_OUTPUT
  if ((index==4)||(index==2))
    cout << "Acting on "<<index<<": "<<r1<<" "<<r2<<endl;
  cout << "ProdATemporaryState="<<this->ProdATemporaryState[0];
  for (int i=1; i<NbrStates; ++i) cout << " " <<this->ProdATemporaryState[i];
  cout << " (h.b. "<<ProdATemporaryStateHighestBit<<")"<<endl;
#endif
  return sqrt(Coefficient);
}


// apply a^+_q1 a^+_q2 operator to the state produced using AA method (without destroying it)
//
// q1 = first index for creation operator
// q2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnLatticeGenericMomentumSpace::AdAd (int q1, int q2, double& coefficient)
{
  for (int i = 0; i < this->NbrStates; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  TemporaryStateHighestBit=ProdATemporaryStateHighestBit;
  ++this->TemporaryState[q2];
  coefficient = this->TemporaryState[q2];
  if (q2 > TemporaryStateHighestBit)
    TemporaryStateHighestBit = q2;
  ++this->TemporaryState[q1];
  coefficient *= this->TemporaryState[q1];
  if (q1 > TemporaryStateHighestBit)
    TemporaryStateHighestBit = q1;
  coefficient = sqrt(coefficient);
#ifdef DEBUG_OUTPUT
  cout << "Temporary State="<<this->TemporaryState[0];
  for (int i=1; i<NbrStates; ++i) cout << " " <<this->TemporaryState[i];
  cout << " (h.b. "<<TemporaryStateHighestBit<<")"<<endl;
  unsigned long TargetState=this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit);
  int TargetIndex=this->FindStateIndex(TargetState, this->TemporaryStateHighestBit + this->NbrBosons - 1);
  bitset<32> b=TargetState;
  bitset<32> b2=HardCoreBasis->StateDescription[TargetIndex];
  cout << q1 << " " << q2 <<" Binary Target: "<<b<<" index "<<TargetIndex<<" vs "<<b2<<endl;
  return TargetIndex;
#else
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
#endif
}


// apply a^+_q a_q operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_q a_q

double BosonOnLatticeGenericMomentumSpace::AdA (int index, int q)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  if (q > this->NbrStates)  
    return 0.0;
  return (double) (this->TemporaryState[q]);  
}

// apply a^+_q a_r operator to a given state 
//
// index = index of the state on which the operator has to be applied
// q = index of the creation operator
// r = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnLatticeGenericMomentumSpace::AdA (int index, int q, int r, double& coefficient)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  if ((r > this->TemporaryStateHighestBit) || (this->TemporaryState[r] == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = this->TemporaryState[r];
  --this->TemporaryState[r];
  if ((TemporaryStateHighestBit==r)&&(this->TemporaryState[r]==0)&&(TemporaryStateHighestBit>0))
    {
      --TemporaryStateHighestBit;
      while ((TemporaryState[TemporaryStateHighestBit] == 0)&&(TemporaryStateHighestBit>0))
	--TemporaryStateHighestBit;
    }
  if (q > TemporaryStateHighestBit)
    {
      for (int i = this->TemporaryStateHighestBit + 1; i <= q; ++i)
	this->TemporaryState[i] = 0ul;
      TemporaryStateHighestBit = q;
    }
  ++this->TemporaryState[q];
  coefficient *= this->TemporaryState[q];
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateHighestBit), this->TemporaryStateHighestBit + this->NbrBosons - 1);
}

// apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
// index = index of the state on which the operator has to be applied
// NbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
// QValues = array of quantum numbers where an interaction is present
// InteractionPerQ = coefficient U_q of the interaction
//
double BosonOnLatticeGenericMomentumSpace::AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  this->FermionToBoson(this->StateDescription[index], this->StateHighestBit[index], this->TemporaryState, this->TemporaryStateHighestBit);
  int BosonsLeft=this->NbrBosons;
  double result=0.0;
  if (nbrInteraction==NbrStates)
    {
      int Q=TemporaryStateHighestBit;
      while ((Q>-1) && (BosonsLeft>0))
	{
	  if (this->TemporaryState[Q]!=0)
	    {
	      result+=interactionPerQ[Q]*TemporaryState[Q]*(TemporaryState[Q]-1);
	      BosonsLeft-=TemporaryState[Q];
	    }
	  --Q;
	}
    }
  else // cannot assume ordered set of qValues
    {
      for (int i=0; (i<nbrInteraction)&&(BosonsLeft>0); ++i)
	{
	  int Q=qValues[i];
	  result+=interactionPerQ[Q]*TemporaryState[Q]*(TemporaryState[Q]-1);
	  BosonsLeft-=TemporaryState[Q];	  
	}
    }
  return result;
}


// obtain a list of quantum numbers in state
// index = index of many-body state to be considered
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
// normalization = indicating the multiplicity of the state for bosonic spaces
void BosonOnLatticeGenericMomentumSpace::ListQuantumNumbers(int index, int *quantumNumbers, double &normalization)
{
  normalization=1.0;
  this->ListQuantumNumbers(index, quantumNumbers);
}

// obtain a list of quantum numbers in state
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
void BosonOnLatticeGenericMomentumSpace::ListQuantumNumbers(int index, int *quantumNumbers)
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
int BosonOnLatticeGenericMomentumSpace::TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase)
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
bool BosonOnLatticeGenericMomentumSpace::IsTranslation(int i, int f, int &shiftX, int &shiftY)
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

Complex BosonOnLatticeGenericMomentumSpace::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
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

Complex BosonOnLatticeGenericMomentumSpace::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
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
Complex BosonOnLatticeGenericMomentumSpace::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
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

Complex BosonOnLatticeGenericMomentumSpace::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, 
								 int nextCoordinates, int firstComponent, 
								 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void BosonOnLatticeGenericMomentumSpace::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnLatticeGenericMomentumSpace::PrintState (ostream& Str, int state)
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

int BosonOnLatticeGenericMomentumSpace::FindStateIndex(unsigned long stateDescription, int highestBit)
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
int BosonOnLatticeGenericMomentumSpace::CarefulFindStateIndex(unsigned long stateDescription, int highestBit)
{
//   bitset<20> b=stateDescription;
//   cout << "Searching " << b <<", bitcount="<<bitcount(stateDescription);
  if (bitcount(stateDescription)!=this->NbrBosons)
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
// nbrBosons = number of fermions
// currentQ = current largest combined quantum number
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonOnLatticeGenericMomentumSpace::EvaluateHilbertSpaceDimension(int nbrBosons, int currentQ, int currentTotalKx, int currentTotalKy)
{
  if ((nbrBosons == 0) || (currentQ < 0))
    return 0l;
  int CurrentKx, CurrentKy, CurrentBand;
  DecodeQuantumNumber(currentQ, CurrentKx, CurrentKy, CurrentBand);
  long Count = 0;

  if (nbrBosons == 1)
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
  Count = 0; 
  while (nbrBosons > 0)
    {
      Count += this->EvaluateHilbertSpaceDimension(nbrBosons, ReducedCurrentQ, currentTotalKx, currentTotalKy);
      --nbrBosons;
      currentTotalKx+=CurrentKx;
      currentTotalKy+=CurrentKy;
    }
  return Count;
}

long BosonOnLatticeGenericMomentumSpace::GenerateStates()
{
  HardCoreBoson *FullBasis = new HardCoreBoson(this->NbrBosons, this->NbrStates + this->NbrBosons - 1, 0);
  long FullDim = FullBasis->GetLargeHilbertSpaceDimension();
  long Pos=0;
  int *MomentaKx = new int[this->NbrStates];
  int *MomentaKy = new int[this->NbrStates];
  int TmpBand;
  int StateKx, StateKy;
  for (int q=0; q<NbrStates; ++q)
    DecodeQuantumNumber(q, MomentaKx[q], MomentaKy[q], TmpBand);

  for (long i=0; i<FullDim; ++i)
    {
      this->FermionToBoson(FullBasis->StateDescription[i], FullBasis->StateHighestBit[i], this->TemporaryState, this->TemporaryStateHighestBit);
      StateKx=0;
      StateKy=0;
      for (int q=0; q<this->TemporaryStateHighestBit; ++q)
	{
	  if (this->TemporaryState[q]!=0)
	    {
	      StateKx+=this->TemporaryState[q]*MomentaKx[q];
	      StateKy+=this->TemporaryState[q]*MomentaKy[q];
	    }
	}
      if ( ((StateKx % this->Nx)==this->Kx) && ((StateKy % this->Ny)==this->Ky))
	{
	  this->StateDescription[Pos] = FullBasis->StateDescription[i];
	  this->StateHighestBit[Pos] = FullBasis->StateHighestBit[i];
	  ++Pos;
	}
    }
  delete [] MomentaKx;
  delete [] MomentaKy;
  delete FullBasis;
  return Pos;
}
  



// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnLatticeGenericMomentumSpace::GenerateLookUpTable(unsigned long memory)
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

