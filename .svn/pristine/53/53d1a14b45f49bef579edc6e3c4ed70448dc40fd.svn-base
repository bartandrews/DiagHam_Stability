////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                 class of quatum Hall hamiltonian associated                //
//   to particles with contact interactions on a lattice in magnetic field    //
//                                                                            //
//                      last modification : 13/02/2008                        //
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
#include "Hamiltonian/ParticleOnLatticeKapitMuellerHamiltonian.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/StringTools.h"
#include <iostream>
#include <algorithm>
#include <cassert>
using std::cout;
using std::endl;

using std::ostream;

// switch for debugging output:
// #define DEBUG_OUTPUT



// constructor for contact interactions on a square lattice
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lx = length of simulation cell in x-direction
// ly = length of simulation cell in y-direction
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// contactInteractionU = strength of on-site delta interaction
// reverseHopping = flag to indicate if sign of hopping terms should be reversed
// deltaPotential = strength of a delta potential at site (0,0)
// randomPotential = magnitude of random potential to add to all sites
// range = maximum length cut-off for single-particle hoppings
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// hermitianFlag = flag indicating whether to use hermitian symmetry
ParticleOnLatticeKapitMuellerHamiltonian::ParticleOnLatticeKapitMuellerHamiltonian(ParticleOnLattice* particles, int nbrParticles, int lx, int ly, int nbrFluxQuanta, double contactInteractionU, bool reverseHopping, double deltaPotential, double randomPotential, double range, AbstractArchitecture* architecture, int nbrBody, unsigned long memory, char* precalculationFileName, bool hermitianFlag)
{
  cout << "creating ParticleOnLatticeKapitMuellerHamiltonian:\n"
       << "nbrParticles="<<nbrParticles<<", lx="<<lx<<", ly="<<ly<<", nbrFluxQuanta="<<nbrFluxQuanta<<"\nU="<<contactInteractionU<<", reverseHopping="<<reverseHopping<<", deltaPotential="<<deltaPotential<<", randomPotential="<<randomPotential<<", range="<<range<<", nbrBody="<<nbrBody<<", hermitianFlag = "<<hermitianFlag<<endl;
  this->Particles=particles;
  this->NbrParticles=nbrParticles;
  this->Lx=lx;
  this->Ly=ly;
  this->NbrSublattices=1;
  this->HaveKySymmetry=false;
  this->KyMax=0;  
  this->NbrCells=lx*ly;
  this->NbrSites=NbrCells*NbrSublattices;
  this->NbrFluxQuanta=nbrFluxQuanta;
  this->HamiltonianShift=0.0;
  this->FluxDensity=((double)nbrFluxQuanta)/NbrCells;
  this->ContactInteractionU=contactInteractionU;
  this->ReverseHopping = reverseHopping;
  this->DeltaPotential = deltaPotential;
  this->RandomPotential = randomPotential;
  this->Range = range;
  if (range < 1.0)
    {
      std::cerr << "Error: range of Kapit-Mueller Hamiltonian needs to be >=1.0"<<endl;
      exit(1);
    }
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  if (hermitianFlag)
    this->HermitianSymmetrizeInteractionFactors();
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  PrintMemorySize(cout, TmpMemory);
	  cout << endl;
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
ParticleOnLatticeKapitMuellerHamiltonian::~ParticleOnLatticeKapitMuellerHamiltonian()
{
  if (NbrHoppingTerms>0)
    {
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
    }
  if (NbrInteractionFactors>0)
    {
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      delete [] this->Q3Value;
      delete [] this->Q4Value;
    }
  if (NbrQ12Indices>0)
    {
      for (int i=0; i<NbrQ12Indices; ++i)
	{
	  delete [] this->Q3PerQ12[i];
	  delete [] this->Q4PerQ12[i];
	}
      delete [] this->NbrQ34Values;
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      
    }
  if (NbrDiagonalInteractionFactors>0)
    {
      delete [] this->DiagonalInteractionFactors;
      delete [] this->DiagonalQValues;
    }
}


// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream
ostream& operator << (ostream& Str, ParticleOnLatticeKapitMuellerHamiltonian& H)
{
  Str << "Need to implement ostream& operator << for ParticleOnLatticeKapitMuellerHamiltonian!" << endl;
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream
MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeKapitMuellerHamiltonian& H)
{
  Str << "Need to implement MathematicaOutput& operator << for ParticleOnLatticeKapitMuellerHamiltonian!\n";
  return Str;
}


// evaluate all interaction factors
//   
void ParticleOnLatticeKapitMuellerHamiltonian::EvaluateInteractionFactors()
{  
  // hopping terms are present independent of statistics:
  this->NbrHoppingTerms=this->NbrSites*this->NbrSites;
  if (this->DeltaPotential != 0.0) ++this->NbrHoppingTerms;
  if (this->RandomPotential != 0.0) this->NbrHoppingTerms+=this->NbrSites;  
  this->HoppingTerms = new Complex[NbrHoppingTerms];
  this->KineticQi = new int[NbrHoppingTerms];
  this->KineticQf = new int[NbrHoppingTerms];

  int TmpNumberTerms=0;
  double HoppingSign = (this->ReverseHopping ? 1.0 : -1.0);
  Complex TranslationPhase;

  // manual parameters, for now
  int images = 5; // number of images of the simulation cell
  int maxRange = std::sqrt(30./(1.0-this->FluxDensity)); // larger exponents will yields numerical zero
  if (this->Range > maxRange)
    {
      // cout << "Testing: range not reduced to new maxRange="<<maxRange<<endl;
      this->Range = maxRange;
    }
  images = this->Range / (Lx < Ly? Lx:Ly) + 2;
  Complex amplitude;
  //loop over initial sites (i,j)

#define FINITE_RANGE
#ifdef FINITE_RANGE
  for (int i=0; i<Lx; ++i)
    {
      for (int j=0; j<Ly; ++j)
	{
	  int qi = Particles->EncodeQuantumNumber(i, j, 0, TranslationPhase);
	      	      
	  // have long-range hopping, so sum over all possible final sites (k,l)
	  for (int k=0; k<Lx; ++k)
	    {
	      for (int l=0; l<Ly; ++l)
		{	
		  KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(k, l, 0, TranslationPhase);
		  KineticQi[TmpNumberTerms] = qi;

		  HoppingTerms[TmpNumberTerms] = 0.0;
		  for (int dX = -images; dX <= images; ++dX)
		    for (int dY = -images; dY <= images; ++dY)
		      {
			if ( this->GetDistance(k+dX*Lx - i, l+dY*Ly - j) <= this->Range // allow hard cut-off
			     // && ( ( dX!=0 || dY!=0 ) || ( i!=k || j!=l ) ) // optionally, exclude onsite terms (needed to get zero energy for lowest band)
			     )
			  {
			    int qf = Particles->EncodeQuantumNumber(k+dX*Lx, l+dY*Ly, 0, TranslationPhase);
			    assert(qf==KineticQf[TmpNumberTerms]);
			    amplitude = -HoppingSign*this->KapitMuellerHopping(k+dX*Lx, l+dY*Ly, i, j) * Conj(TranslationPhase);			       
#ifdef DEBUG_OUTPUT
			    //if (TranslationPhase!=1.0)
			    if (Norm(amplitude) > 1e-15)
			      cout << "image ("<<dX<<", "<<dY<<"): sites ("<<i<<", "<<j<<")->("<<k+dX*Lx<<", "<<l+dY*Ly<<") with dL=("<<dX<<"," <<dY<<") : "
				   <<"Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
				   << TranslationPhase << " amp="<< amplitude<<endl;
#endif
			    HoppingTerms[TmpNumberTerms] += amplitude;
			  }
		      }
#ifdef DEBUG_OUTPUT		     
		  if (Norm(HoppingTerms[TmpNumberTerms])>1e-15)
		    cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
		  if (Norm(HoppingTerms[TmpNumberTerms])>1e-15)
		    ++TmpNumberTerms; // only take into account terms with non-zero magnitude
		}
	    }
	}
    }
#else
  for (int i=0; i<Lx; ++i)
    {
      for (int j=0; j<Ly; ++j)
	{
	  int qi = Particles->EncodeQuantumNumber(i, j, 0, TranslationPhase);
	      	      
	  for (int k=0; k<Lx; ++k)
	    {
	      for (int l=0; l<Ly; ++l)
		{	
		  KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(k, l, 0, TranslationPhase);
		  KineticQi[TmpNumberTerms] = qi;

		  HoppingTerms[TmpNumberTerms] = this->SumImagesForHoppings(k,l,i,j,images);
#ifdef DEBUG_OUTPUT
		  if (Norm(HoppingTerms[TmpNumberTerms])>1e-15)
		    cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
		  if (Norm(HoppingTerms[TmpNumberTerms])>1e-15)
		    ++TmpNumberTerms; // only take into account terms with non-zero magnitude
		}
	    }
	}
    }
#endif


  if (this->DeltaPotential != 0.0)
    {
      KineticQi[TmpNumberTerms] = Particles->EncodeQuantumNumber(0, 0, 0, TranslationPhase);
      KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(0, 0, 0, TranslationPhase);
      HoppingTerms[TmpNumberTerms] = this->DeltaPotential;
      //cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
      ++TmpNumberTerms;
    }
  if (this->RandomPotential != 0.0)
    {
      NumRecRandomGenerator G;
      G.UseTimeSeed();
      for (int x=0; x<Lx; ++x)
	for (int y=0; y<Ly; ++y)
	  {
	    KineticQi[TmpNumberTerms] = Particles->EncodeQuantumNumber(x, y, 0, TranslationPhase);
	    KineticQf[TmpNumberTerms] = KineticQi[TmpNumberTerms];
	    HoppingTerms[TmpNumberTerms] = this->RandomPotential*(-0.5+G.GetRealRandomNumber());
	    //cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
	    ++TmpNumberTerms;
	  }
    }


  // reassign tables to minimize memory footprint
  if (TmpNumberTerms < this->NbrHoppingTerms)
    {
        Complex *NewHoppingTerms = new Complex[TmpNumberTerms];
	int *NewKineticQi = new int[TmpNumberTerms];
	int *NewKineticQf = new int[TmpNumberTerms];

	for (int i=0; i<TmpNumberTerms; ++i)	 
	  {
	    NewHoppingTerms[i] = this->HoppingTerms[i];
	    NewKineticQi[i] = this->KineticQi[i];
	    NewKineticQf[i] = this->KineticQf[i];
	  }
	delete [] this->HoppingTerms;
	delete [] this->KineticQi;
	delete [] this->KineticQf;
	this->NbrHoppingTerms = TmpNumberTerms;
	this->HoppingTerms = NewHoppingTerms;
	this->KineticQi = NewKineticQi;
	this->KineticQf = NewKineticQf;
    }

  // we have no general four-particle interactions:
  this->NbrInteractionFactors=0;
  this->NbrQ12Indices=0;
  
  // contact interactions come to play for bosons, only!
  if ((this->Particles->GetParticleStatistic() == ParticleOnLattice::BosonicStatistic) && (this->ContactInteractionU!=0.0))
    {
      cout << "adding interaction terms"<<endl;
      this->NbrDiagonalInteractionFactors = this->NbrSites;
      this->DiagonalInteractionFactors = new double[NbrDiagonalInteractionFactors];
      this->DiagonalQValues=new int[NbrDiagonalInteractionFactors];
      for (int i=0; i<NbrDiagonalInteractionFactors; ++i)
	{
	  this->DiagonalQValues[i]=i;
	  this->DiagonalInteractionFactors[i]=this->ContactInteractionU;
	}
    }
  else // no such interactions
    {
      NbrDiagonalInteractionFactors=0;
    }
}
