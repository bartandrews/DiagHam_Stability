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
#include "Hamiltonian/ParticleOnLatticeDeltaHamiltonian.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
using std::cout;
using std::endl;

using std::ostream;

// switch for debugging output:
//#define DEBUG_OUTPUT



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
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
ParticleOnLatticeDeltaHamiltonian::ParticleOnLatticeDeltaHamiltonian(ParticleOnLattice* particles, int nbrParticles, int lx, int ly, int nbrFluxQuanta, double contactInteractionU, bool reverseHopping, double deltaPotential, double randomPotential, AbstractArchitecture* architecture, int memory, char* precalculationFileName)
{
  this->Particles=particles;
  this->NbrParticles=nbrParticles;
  this->Lx=lx;
  this->Ly=ly;
  this->SubLattices=1;
  this->HaveKySymmetry=false;
  this->KyMax=0;  
  this->NbrCells=lx*ly;
  this->NbrSites=NbrCells*SubLattices;
  this->NbrFluxQuanta=nbrFluxQuanta;
  this->HamiltonianShift=0.0;
  this->FluxDensity=((double)nbrFluxQuanta)/NbrCells;
  this->ContactInteractionU=contactInteractionU;
  this->ReverseHopping = reverseHopping;
  this->DeltaPotential = deltaPotential;
  this->RandomPotential = randomPotential;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  int TmpMemory = this->FastMultiplicationMemory(memory);
	  if (TmpMemory < 1024)
	    cout  << "fast = " <<  TmpMemory << "b ";
	  else
	    if (TmpMemory < (1 << 20))
	      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	    else
	      if (TmpMemory < (1 << 30))
		cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	      else
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
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
ParticleOnLatticeDeltaHamiltonian::~ParticleOnLatticeDeltaHamiltonian()
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
ostream& operator << (ostream& Str, ParticleOnLatticeDeltaHamiltonian& H)
{
  Str << "Need to implement ostream& operator << for ParticleOnLatticeDeltaHamiltonian!" << endl;
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream
MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeDeltaHamiltonian& H)
{
  Str << "Need to implement MathematicaOutput& operator << for ParticleOnLatticeDeltaHamiltonian!\n";
  return Str;
}


// evaluate all interaction factors
//   
void ParticleOnLatticeDeltaHamiltonian::EvaluateInteractionFactors()
{  
  // hopping terms are present independent of statistics:
  this->NbrHoppingTerms=4*this->NbrSites;
  if (this->DeltaPotential != 0.0) ++this->NbrHoppingTerms;
  if (this->RandomPotential != 0.0) this->NbrHoppingTerms+=this->NbrSites;  
  this->HoppingTerms = new Complex[NbrHoppingTerms];
  this->KineticQi = new int[NbrHoppingTerms];
  this->KineticQf = new int[NbrHoppingTerms];

  int TmpNumberTerms=0;
  double HoppingSign = (this->ReverseHopping ? 1.0 : -1.0);
  Complex TranslationPhase;
  switch (this->Particles->GetLandauGaugeAxis())
    {
    case 'y': {
      for (int i=0; i<Lx; ++i)
	{
	  Complex Phase=Polar(1.0,2.0*M_PI*this->FluxDensity*(double)i);
	  for (int j=0; j<Ly; ++j)
	    {
	      KineticQi[TmpNumberTerms] = Particles->EncodeQuantumNumber(i, j, 0, TranslationPhase);
	      KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(i+1, j, 0, TranslationPhase);
	      HoppingTerms[TmpNumberTerms] = HoppingSign*TranslationPhase;
#ifdef DEBUG_OUTPUT
	      if (TranslationPhase!=1.0)
		cout << "(i="<<i<<"->"<<i+1<<") Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
		     <<TranslationPhase<<endl;
	      cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	      ++TmpNumberTerms;
	      KineticQi[TmpNumberTerms] = KineticQi[TmpNumberTerms-1];
	      KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(i-1, j, 0, TranslationPhase);
	      HoppingTerms[TmpNumberTerms] = HoppingSign*TranslationPhase;
#ifdef DEBUG_OUTPUT
	      if (TranslationPhase!=1.0)
		cout << "(i="<<i<<"->"<<i-1<<") Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
		     <<TranslationPhase<<endl;
	      cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	      ++TmpNumberTerms;
	      KineticQi[TmpNumberTerms] = KineticQi[TmpNumberTerms-1];
	      KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(i, j+1, 0, TranslationPhase);
	      HoppingTerms[TmpNumberTerms] = HoppingSign*Conj(Phase)*TranslationPhase;
#ifdef DEBUG_OUTPUT
	      if (TranslationPhase!=1.0)
		cout << "(j="<<j<<"->"<<j+1<<") Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
		     <<TranslationPhase<<endl;
	      cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	      ++TmpNumberTerms;
	      KineticQi[TmpNumberTerms] = KineticQi[TmpNumberTerms-1];
	      KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(i, j-1, 0, TranslationPhase);
	      HoppingTerms[TmpNumberTerms] = HoppingSign*Phase*TranslationPhase;
#ifdef DEBUG_OUTPUT
	      if (TranslationPhase!=1.0)
		cout << "(j="<<j<<"->"<<j-1<<") Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
		     <<TranslationPhase<<endl;
	      cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	      ++TmpNumberTerms;
	    }
	}
      break;
    }
    case 'x': {
      for (int j=0; j<Lx; ++j)
	{
	  Complex Phase=Polar(1.0,-2.0*M_PI*this->FluxDensity*(double)j);
	  for (int i=0; i<Ly; ++i)
	    {
	      KineticQi[TmpNumberTerms] = Particles->EncodeQuantumNumber(i, j, 0, TranslationPhase);
	      KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(i+1, j, 0, TranslationPhase);
	      HoppingTerms[TmpNumberTerms] = HoppingSign*Conj(Phase)*TranslationPhase;
#ifdef DEBUG_OUTPUT
	      if (TranslationPhase!=1.0)
		cout << "(i="<<i<<"->"<<i+1<<") Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
		     <<TranslationPhase<<endl;
	      cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	      ++TmpNumberTerms;
	      KineticQi[TmpNumberTerms] = KineticQi[TmpNumberTerms-1];
	      KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(i-1, j, 0, TranslationPhase);
	      HoppingTerms[TmpNumberTerms] = HoppingSign*Phase*TranslationPhase;
#ifdef DEBUG_OUTPUT
	      if (TranslationPhase!=1.0)
		cout << "(i="<<i<<"->"<<i-1<<") Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
		     <<TranslationPhase<<endl;
	      cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	      ++TmpNumberTerms;
	      KineticQi[TmpNumberTerms] = KineticQi[TmpNumberTerms-1];
	      KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(i, j+1, 0, TranslationPhase);
	      HoppingTerms[TmpNumberTerms] = HoppingSign*TranslationPhase;
#ifdef DEBUG_OUTPUT
	      if (TranslationPhase!=1.0)
		cout << "(j="<<j<<"->"<<j+1<<") Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
		     <<TranslationPhase<<endl;
	      cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	      ++TmpNumberTerms;
	      KineticQi[TmpNumberTerms] = KineticQi[TmpNumberTerms-1];
	      KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(i, j-1, 0, TranslationPhase);
	      HoppingTerms[TmpNumberTerms] = HoppingSign*TranslationPhase;
#ifdef DEBUG_OUTPUT
	      if (TranslationPhase!=1.0)
		cout << "(j="<<j<<"->"<<j-1<<") Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
		     <<TranslationPhase<<endl;
	      cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	      ++TmpNumberTerms;
	    }
	}
      break;
    }
    default:
      cout << "Invalid Landau quantization axis encountered in ParticleOnLatticeDeltaHamiltonian."<<endl;
      exit(1);
      break;
    }

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
	    HoppingTerms[TmpNumberTerms] = this->DeltaPotential*(-0.5+G.GetRealRandomNumber());
	    //cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
	    ++TmpNumberTerms;
	  }
    }

  // we have no general four-particle interactions:
  this->NbrInteractionFactors=0;
  this->NbrQ12Indices=0;
  
  // contact interactions come to play for bosons, only!
  if ((this->Particles->GetParticleStatistic() == ParticleOnLattice::BosonicStatistic) && (this->ContactInteractionU!=0.0))
    {
      cout << "adding interaction terms"<<endl;
      this->NbrDiagonalInteractionFactors=this->NbrSites;
      this->DiagonalInteractionFactors=new double[NbrDiagonalInteractionFactors];
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
