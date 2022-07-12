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
#include "GeneralTools/StringTools.h"

#include <iostream>
using std::cout;
using std::endl;

using std::ostream;

// switch for debugging output:
#define DEBUG_OUTPUT



// constructor for contact interactions on a square lattice
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lx = length of simulation cell in x-direction
// ly = length of simulation cell in y-direction
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// contactInteractionU = strength of on-site delta interaction for bosons / NN interaction for fermions
// reverseHopping = flag to indicate if sign of hopping terms should be reversed
// deltaPotential = strength of a delta potential at site (0,0)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// hermitianFlag = flag indicating whether to use hermitian symmetry
// cylinder_geometry = flag indicating whether to omit periodic boundary condition in the x-direction
ParticleOnLatticeDeltaHamiltonian::ParticleOnLatticeDeltaHamiltonian(ParticleOnLattice* particles, int nbrParticles, int lx, int ly, int nbrFluxQuanta, double contactInteractionU, bool reverseHopping, double deltaPotential, double randomPotential, AbstractArchitecture* architecture, int nbrBody, unsigned long memory, char* precalculationFileName, bool hermitianFlag, bool cylinder_geometry)
{
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
  this->CylinderGeometry = cylinder_geometry;
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
	  PrintMemorySize(cout, TmpMemory) << endl;
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
	      if (!this->CylinderGeometry || i < Lx-1)
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
		}
	      if (!this->CylinderGeometry || i > 0)
		{
		  KineticQi[TmpNumberTerms] = Particles->EncodeQuantumNumber(i, j, 0, TranslationPhase);
		  KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(i-1, j, 0, TranslationPhase);
		  HoppingTerms[TmpNumberTerms] = HoppingSign*TranslationPhase;
#ifdef DEBUG_OUTPUT
		  if (TranslationPhase!=1.0)
		    cout << "(i="<<i<<"->"<<i-1<<") Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
			 <<TranslationPhase<<endl;
		  cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
		  ++TmpNumberTerms;
		}
	      KineticQi[TmpNumberTerms] = Particles->EncodeQuantumNumber(i, j, 0, TranslationPhase);
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
      for (int j=0; j<Ly; ++j)
	{
	  Complex Phase=Polar(1.0,-2.0*M_PI*this->FluxDensity*(double)j);
	  for (int i=0; i<Lx; ++i)
	    {
	      if (!this->CylinderGeometry || i < Lx-1)
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
		}
	      if (!this->CylinderGeometry || i > 0)
		{		  
		  KineticQi[TmpNumberTerms] = Particles->EncodeQuantumNumber(i, j, 0, TranslationPhase);
		  KineticQf[TmpNumberTerms] = Particles->EncodeQuantumNumber(i-1, j, 0, TranslationPhase);
		  HoppingTerms[TmpNumberTerms] = HoppingSign*Phase*TranslationPhase;
#ifdef DEBUG_OUTPUT
		  if (TranslationPhase!=1.0)
		    cout << "(i="<<i<<"->"<<i-1<<") Translation ["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="
			 <<TranslationPhase<<endl;
		  cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<" tP="<<TranslationPhase<<endl;
#endif
		  ++TmpNumberTerms;
		}
	      KineticQi[TmpNumberTerms] = Particles->EncodeQuantumNumber(i, j, 0, TranslationPhase);
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
	    HoppingTerms[TmpNumberTerms] = this->RandomPotential*(-0.5+G.GetRealRandomNumber());
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
      this->NbrDiagonalInteractionFactors = this->NbrSites;
      this->NbrRhoRhoInteractionFactors = 0;
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
      if ((this->Particles->GetParticleStatistic() == ParticleOnLattice::FermionicStatistic) && (this->Particles->GetMinNbrParticles()>1))
	{
	  this->NbrDiagonalInteractionFactors=0;
	  this->NbrRhoRhoInteractionFactors = 2*this->NbrSites - this->Ly * (this->CylinderGeometry==true);
	  this->RhoRhoInteractionFactors = new double[this->NbrRhoRhoInteractionFactors];
	  this->RhoRhoQ12Values = new int[2*this->NbrRhoRhoInteractionFactors];
	  int Qi, posQ=0, posU=0;
	  Complex TranslationPhase;
	  for (int i=0; i<Lx; ++i)
	    for (int j=0; j<Ly; ++j)
	      {
		Qi = Particles->EncodeQuantumNumber(i, j, 0, TranslationPhase);
		// horizontal interactions
		if (!this->CylinderGeometry || i < Lx-1)
		  {		    
		    this->RhoRhoQ12Values[posQ++]=Qi;
		    this->RhoRhoQ12Values[posQ++]=Particles->EncodeQuantumNumber(i+1, j, 0, TranslationPhase);
		    this->RhoRhoInteractionFactors[posU++]=this->ContactInteractionU;
		  }
		// vertical interactions
		this->RhoRhoQ12Values[posQ++]=Qi;
		this->RhoRhoQ12Values[posQ++]=Particles->EncodeQuantumNumber(i, j+1, 0, TranslationPhase);
		this->RhoRhoInteractionFactors[posU++]=this->ContactInteractionU;
	      }
	}
      else
	{
	  this->NbrRhoRhoInteractionFactors=0;
	  this->NbrDiagonalInteractionFactors=0;
	}
    }
}
