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
#include "Hamiltonian/ParticleOnLatticeGenericHamiltonian.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
using std::cout;
using std::endl;

using std::ostream;

// switch for debugging output:
#define DEBUG_GENERIC_LATTICE_HAMILTONIAN



// constructor for contact interactions on a square lattice
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// latticeGeometry = geometry of lattice system is living on
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// contactInteractionU = strength of on-site delta interaction
// reverseHopping = flag to indicate if sign of hopping terms should be reversed
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// hermitianFlag = flag indicating whether to use hermitian symmetry
// hoppingOnly = evaluate only energy of hopping terms, excluding local potentials
ParticleOnLatticeGenericHamiltonian::ParticleOnLatticeGenericHamiltonian(ParticleOnLattice* particles, int nbrParticles, LatticePhases *latticeGeometry, int nbrFluxQuanta, double contactInteractionU, bool reverseHopping, AbstractArchitecture* architecture, unsigned long memory, char* precalculationFileName, double overrideFluxDensity, bool hermitianFlag, bool hoppingOnly)
{
  this->Particles=particles;
  this->NbrParticles=nbrParticles;
  this->LatticeGeometry = latticeGeometry;
  this->LatticeDimension = LatticeGeometry->GetLatticeDimension();
  this->NbrSublattices = LatticeGeometry->GetNbrSubLattices();
  this->NbrCells=LatticeGeometry->GetNbrCells();
  this->Length = new int[LatticeDimension];
  for (int i=0; i<LatticeDimension; ++i)
    this->Length[i] = LatticeGeometry->GetLatticeLength(i);
  this->HaveKySymmetry=false;
  this->KyMax=0;  
  this->NbrSites = this->LatticeGeometry->GetNbrSites();
  this->NbrFluxQuanta=nbrFluxQuanta;
  this->HamiltonianShift=0.0;
  if (LatticeGeometry->HavePredefinedFlux())
    this->FluxDensity=(double)(LatticeGeometry->GetPredefinedFlux())/(double)(LatticeGeometry->GetNbrSites());
  else
    this->FluxDensity=((double)nbrFluxQuanta)/NbrSites;
  if (overrideFluxDensity!=0.0)
    {
      if (LatticeGeometry->AllowContinuousPhases())
	{
	  this->FluxDensity=overrideFluxDensity;
	}
      else
	{
	  cout << "This lattice does not support non-quantized fluxes"<<endl;
	  exit(1);
	}
    }
  this->ContactInteractionU=contactInteractionU;
  this->ReverseHopping = reverseHopping;
  this->Architecture = architecture;
  this->HoppingOnly = hoppingOnly;
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
	  PrintMemorySize(cout, TmpMemory)<< endl;
	  if (memory > 0)
	    this->EnableFastMultiplication();
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//
ParticleOnLatticeGenericHamiltonian::~ParticleOnLatticeGenericHamiltonian()
{
  delete [] this->Length;
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
ostream& operator << (ostream& Str, ParticleOnLatticeGenericHamiltonian& H)
{
  Str << "Need to implement ostream& operator << for ParticleOnLatticeGenericHamiltonian!" << endl;
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream
MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeGenericHamiltonian& H)
{
  Str << "Need to implement MathematicaOutput& operator << for ParticleOnLatticeGenericHamiltonian!\n";
  return Str;
}


// evaluate all interaction factors
//   
void ParticleOnLatticeGenericHamiltonian::EvaluateInteractionFactors()
{  
  // hopping terms are present independent of statistics:
  this->NbrHoppingTerms=LatticeGeometry->GetNbrHoppingTerms();
  if (!this->HoppingOnly)
    this->NbrHoppingTerms+=LatticeGeometry->GetNbrLocalPotentials();
  this->HoppingTerms = new Complex[NbrHoppingTerms];
  this->KineticQi = new int[NbrHoppingTerms];
  this->KineticQf = new int[NbrHoppingTerms];

  int TmpNumberTerms=0;
  int *Neighbors;
  double *Phases;
  int **Translations;
  double *Amplitudes;
  int NbrNeighbors;
  double HoppingSign = (this->ReverseHopping ? 1.0 : -1.0);
  double SolenoidX, SolenoidY;
  Particles->GetSolenoidFluxes(SolenoidX, SolenoidY);
  
  for (int s=0; s<NbrSites; ++s)
    {
      this->LatticeGeometry->GetNeighbors(s, NbrNeighbors, Neighbors, Phases, Translations, Amplitudes);
      for (int n=0; n<NbrNeighbors; ++n)
	{
	  KineticQi[TmpNumberTerms] = s;
	  KineticQf[TmpNumberTerms] = Neighbors[n];
#ifdef DEBUG_GENERIC_LATTICE_HAMILTONIAN
	  if (LatticeGeometry->HavePredefinedFlux())
	    cout << "Using hard-coded flux-density "<<this->FluxDensity<<" and Phase="<<Phases[n]<<endl;
	  else
	    cout << "Using flux density="<<this->FluxDensity<<" and Phase="<<Phases[n]<<endl;
#endif
	  
	  HoppingTerms[TmpNumberTerms] = Polar(1.0,Translations[n][0]*SolenoidX)*Polar(1.0,Translations[n][1]*SolenoidY)*HoppingSign*
	    Polar(1.0,-2.0*M_PI*this->FluxDensity*Phases[n]);
	  if (Amplitudes!=NULL)
 	    HoppingTerms[TmpNumberTerms] *= Amplitudes[n];
#ifdef DEBUG_GENERIC_LATTICE_HAMILTONIAN
	  cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms];
	  cout << " including solenoid flux "<< Polar(1.0,Translations[n][0]*SolenoidX)<<"_x, ";
	  cout << Polar(1.0,Translations[n][1]*SolenoidY)<<"_y ("<<Translations[n][0]<<","<<Translations[n][1]<<")"<<endl;
#endif
	  ++TmpNumberTerms;
	}
    }

  if ((!this->HoppingOnly)&&(LatticeGeometry->HaveOneParticlePotentials()))
    {
      cout << "Adding one particle potentials in Hamiltonian"<<endl;
      int NbrPotentials;
      int *PotentialPositions;
      double *Potentials = LatticeGeometry->GetOneParticlePotentials(NbrPotentials, PotentialPositions);
      for (int n=0; n<NbrPotentials; ++n)
	{
	  KineticQi[TmpNumberTerms] = PotentialPositions[n];
	  KineticQf[TmpNumberTerms] = KineticQi[TmpNumberTerms];
	  HoppingTerms[TmpNumberTerms] = Potentials[n];
#ifdef DEBUG_GENERIC_LATTICE_HAMILTONIAN
	  cout << "H["<<KineticQi[TmpNumberTerms]<<"->"<<KineticQf[TmpNumberTerms]<<"]="<<HoppingTerms[TmpNumberTerms]<<endl;
#endif
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
