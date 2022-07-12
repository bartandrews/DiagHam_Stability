////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
//                                                                            //
//                                                                            //
//                       class of Hamiltonian H=S*S on lattice                //
//                                                                            //
//                        last modification : 31/07/2009                      //
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
#include "SUNSpinOnLatticeQuadraticHamiltonian.h"
#include "HilbertSpace/GenericSUNSpinCollection.h"
#include "Hamiltonian/AbstractSUNSpinOnLatticeHamiltonian.h"
#include "Tools/LatticeConnections.h"
#include "GeneralTools/StringTools.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/SUNSpinPrecalculationOperation.h"



// constructor
// space = Hilbert space for problem
// lattice = class providing size and geometry / connections on lattice
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// cyclic = prefactor of optional cyclic permutations around plaquettes
SUNSpinOnLatticeQuadraticHamiltonian::SUNSpinOnLatticeQuadraticHamiltonian(GenericSUNSpinCollection *space, LatticeConnections *lattice, AbstractArchitecture* architecture, long memory, char* precalculationFileName, double cyclic)
{ 
  this->Spins = space;
  this->Lattice = lattice;
  this->NbrSpins = Lattice->GetNbrSites();
  this->FastMultiplicationFlag = false;
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  this->CyclicPrefactor = cyclic;
  this->EvaluateInteractionTerms();
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  this->DiskStorageFlag = false;
  
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout  << "fast = ";
	  PrintMemorySize(cout, TmpMemory)<<endl;
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);

}

// destructor
//
SUNSpinOnLatticeQuadraticHamiltonian::~SUNSpinOnLatticeQuadraticHamiltonian()
{
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian
//SUNSpinOnLatticeQuadraticHamiltonian::SUNSpinOnLatticeQuadraticHamiltonian* Clone ()
//{
//}


// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SUNSpinOnLatticeQuadraticHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Spins->GetHilbertSpaceDimension();
}


// evaluate all interaction factors
//   
void SUNSpinOnLatticeQuadraticHamiltonian::EvaluateInteractionTerms()
{
  this->HaveComplexInteractions = false;
  int *Partners;
  int NbrPartners;
  int TmpNbrPermutationTerms=0;
  for (int s=0; s<NbrSpins; ++s)
    {      
      this->Lattice->GetPartners(s, Partners, NbrPartners);
      TmpNbrPermutationTerms+=NbrPartners;
    }
  this->NbrPermutationTerms=TmpNbrPermutationTerms;
  this->PermutationPrefactors = new double[this->NbrPermutationTerms];
  this->PermutationI = new int[this->NbrPermutationTerms];
  this->PermutationJ = new int[this->NbrPermutationTerms];
  int Pos=0;
  for (int s=0; s<NbrSpins; ++s)
    {      
      this->Lattice->GetPartners(s, Partners, NbrPartners);
      for (int p=0; p<NbrPartners; ++p)
	{
	  this->PermutationI[Pos]=s;
	  this->PermutationJ[Pos]=Partners[p];
	  this->PermutationPrefactors[Pos]=1.0;
	  ++Pos;
	}
    }
  if (fabs((double)this->CyclicPrefactor)>0.0)
    {
      this->NbrCyclicPermutations = 2*this->Lattice->GetNbrPlaquettes();
      if (this->NbrCyclicPermutations==0)
	{
	  cout << "Attention, no plaquettes defined for this lattice!"<<endl;
	  return;
	}
      this->CyclicPermutationPrefactors = new Complex[NbrCyclicPermutations];
      this->CyclicPermutationLength = new int[NbrCyclicPermutations];
      this->CyclicPermutationIndices = new int*[NbrCyclicPermutations];
      for (int i=0; i<NbrCyclicPermutations>>1; ++i)
	{
	  this->Lattice->GetPlaquetteSpins(i, CyclicPermutationIndices[2*i], CyclicPermutationLength[2*i]);
	  this->CyclicPermutationLength[2*i+1]=this->CyclicPermutationLength[2*i];
	  int l=this->CyclicPermutationLength[2*i+1];
	  this->CyclicPermutationIndices[2*i+1]=new int[l];
	  for (int k=0; k<l; ++k)
	    this->CyclicPermutationIndices[2*i+1][k]=this->CyclicPermutationIndices[2*i][l-1-k];
	  this->CyclicPermutationPrefactors[2*i] = Complex(0.0, this->CyclicPrefactor);
	  this->CyclicPermutationPrefactors[2*i+1] = Complex(0.0, -this->CyclicPrefactor);
	}
    }
}
