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
#include "Hamiltonian/ParticleOnLatticeWithKyDeltaHamiltonian.h"
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
// kyMax = maximum value of momentum in y-direction
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// contactInteractionU = strength of on-site delta interaction
// reverseHopping = flag to indicate if sign of hopping terms should be reversed
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
ParticleOnLatticeWithKyDeltaHamiltonian::ParticleOnLatticeWithKyDeltaHamiltonian(ParticleOnLattice* particles, int nbrParticles, int lx, int ly, int kyMax, int nbrFluxQuanta, double contactInteractionU, bool reverseHopping, double randomPotential, AbstractArchitecture* architecture, int memory, char* precalculationFileName)
{
  this->Particles=particles;
  this->NbrParticles=nbrParticles;
  this->Lx=lx;
  this->Ly=ly;
  this->SubLattices=1;
  this->HaveKySymmetry=true;
  this->KyMax=kyMax;  
  this->NbrCells=lx*ly;
  this->NbrSites=NbrCells*SubLattices;
  this->NbrFluxQuanta=nbrFluxQuanta;
  this->HamiltonianShift=0.0;
  this->FluxDensity=((double)nbrFluxQuanta)/NbrCells;
  this->ContactInteractionU=contactInteractionU;
  this->ReverseHopping = reverseHopping;
  this->DeltaPotential = 0.0;
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
ParticleOnLatticeWithKyDeltaHamiltonian::~ParticleOnLatticeWithKyDeltaHamiltonian()
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
ostream& operator << (ostream& Str, ParticleOnLatticeWithKyDeltaHamiltonian& H)
{
  Str << "Need to implement ostream& operator << for ParticleOnLatticeWithKyDeltaHamiltonian!" << endl;
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream
MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeWithKyDeltaHamiltonian& H)
{
  Str << "Need to implement MathematicaOutput& operator << for ParticleOnLatticeWithKyDeltaHamiltonian!\n";
  return Str;
}


// evaluate all interaction factors
//
void ParticleOnLatticeWithKyDeltaHamiltonian::EvaluateInteractionFactors()
{  
  // hopping terms are present independent of statistics:
  int MaxNumberTerms=4*this->NbrSites;

  if (this->DeltaPotential != 0.0) ++MaxNumberTerms;
  if (this->RandomPotential != 0.0) MaxNumberTerms+=this->NbrSites;  
  this->HoppingTerms = new Complex[MaxNumberTerms];
  this->KineticQi = new int[MaxNumberTerms];
  this->KineticQf = new int[MaxNumberTerms];

  this->NbrHoppingTerms=0;
  
  double HoppingSign = (this->ReverseHopping ? 1.0 : -1.0);
  Complex TranslationPhase;
  int p = Ly/KyMax;
  int compositeKy, compositeKy2;
  for (int i=0; i<Lx; ++i) 
    {
      Complex Phase=Polar(1.0,2.0*M_PI*this->FluxDensity*(double)i);
      for (int k=0; k<KyMax; ++k)
	for (int s=0; s<p; ++s)
	  {
	    compositeKy = k*p+s;
	    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy, 0, TranslationPhase);
	    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i+1, compositeKy, 0, TranslationPhase);
	    HoppingTerms[this->NbrHoppingTerms] = HoppingSign*TranslationPhase;	    

#ifdef DEBUG_OUTPUT
	    if (TranslationPhase!=1.0)
	      cout << "(i="<<i<<"->"<<i+1<<") Translation ["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<TranslationPhase<<endl;
	    cout << "k="<<k<<", "<<"s="<<s<<" x_0="<<i<<endl;
	    cout << "x: H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	    ++this->NbrHoppingTerms;
	    
	    KineticQi[this->NbrHoppingTerms] = KineticQi[this->NbrHoppingTerms-1];
	    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i-1, compositeKy, 0, TranslationPhase);
	    HoppingTerms[this->NbrHoppingTerms] = HoppingSign*TranslationPhase;

#ifdef DEBUG_OUTPUT
	    if (TranslationPhase!=1.0)
	    cout << "(i="<<i<<"->"<<i-1<<") Translation ["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="
	    	 <<TranslationPhase<<endl;	    
	    cout << "x: H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<" tP="<<TranslationPhase<<endl;
#endif
	    ++this->NbrHoppingTerms;

	    // coupling in y-direction: distinguish p=0 and p>0!
	    if (p==1)
	      {
		KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, k, 0, TranslationPhase);
		KineticQf[this->NbrHoppingTerms] = KineticQi[this->NbrHoppingTerms];
		HoppingTerms[this->NbrHoppingTerms] = HoppingSign*2.0*cos(2.0*M_PI*((double)k/Ly-this->FluxDensity*(double)i));
#ifdef DEBUG_OUTPUT
		cout << "y - p=1: H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<endl;
#endif
		++this->NbrHoppingTerms;
	      }
	    else
	      {
		if (s<p-1)
		  {
		    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy, 0, TranslationPhase);
		    compositeKy2 = k*p+s+1;
		    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy2, 0, TranslationPhase);
		    HoppingTerms[this->NbrHoppingTerms] = HoppingSign*Conj(Phase)*TranslationPhase;
		  }
		else		  
		  {
		    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy, 0, TranslationPhase);
		    compositeKy2 = k*p;
		    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy2, 0, TranslationPhase);
		    HoppingTerms[this->NbrHoppingTerms] = HoppingSign*Phase*TranslationPhase
		      *Polar(1.0,-2.0*M_PI*((double)k/KyMax-this->FluxDensity*(double)i));
		  }		
#ifdef DEBUG_OUTPUT
		if (TranslationPhase!=1.0)
		  cout << "(sk)="<<compositeKy<<"->"<<compositeKy2<<") Translation ["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="
				 <<TranslationPhase<<endl;
		cout << "y - p>1: H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<" tP="<<TranslationPhase<<endl;
#endif
		++this->NbrHoppingTerms;
		if (s==0)
		  {
		    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy, 0, TranslationPhase);
		    compositeKy2 = k*p+p-1;
		    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy2, 0, TranslationPhase);
		    HoppingTerms[this->NbrHoppingTerms] = HoppingSign*Conj(Phase)*TranslationPhase
		      *Polar(1.0,2.0*M_PI*((double)k/KyMax-this->FluxDensity*(double)i));
		  }
		else
		  {
		    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy, 0, TranslationPhase);
		    compositeKy2 = k*p+s-1;
		    KineticQf[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(i, compositeKy2, 0, TranslationPhase);
		    HoppingTerms[this->NbrHoppingTerms] = HoppingSign*Phase*TranslationPhase;
		  }
#ifdef DEBUG_OUTPUT
		if (TranslationPhase!=1.0)
		  cout << "(sk)="<<compositeKy<<"->"<<compositeKy2<<") Translation ["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="
				 <<TranslationPhase<<endl;
		cout << "y - p>1: H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<" tP="<<TranslationPhase<<endl;
#endif
		++this->NbrHoppingTerms;
	      }
	  }
    }
  if (this->RandomPotential != 0.0)
    {
//       NumRecRandomGenerator G;
//       G.UseTimeSeed();
//       for (int x=0; x<Lx; ++x)
// 	for (int y=0; y<Ly; ++y)
// 	  {
// 	    KineticQi[this->NbrHoppingTerms] = Particles->EncodeQuantumNumber(x, y, 0, TranslationPhase);
// 	    KineticQf[this->NbrHoppingTerms] = KineticQi[this->NbrHoppingTerms];
// 	    HoppingTerms[this->NbrHoppingTerms] = this->DeltaPotential*(-0.5+G.GetRealRandomNumber());
//#ifdef DEBUG_OUTPUT
// 	    //cout << "H["<<KineticQi[this->NbrHoppingTerms]<<"->"<<KineticQf[this->NbrHoppingTerms]<<"]="<<HoppingTerms[this->NbrHoppingTerms]<<" tP="<<TranslationPhase<<endl;
//#endif
// 	    ++this->NbrHoppingTerms;
// 	  }
    }

  // interactions are not fully diagonal in k-space:
  this->NbrDiagonalInteractionFactors=0;
  this->DiagonalInteractionFactors=0;
  this->DiagonalQValues=0;
  
  // contact interactions come to play for bosons, only!
  if ((this->Particles->GetParticleStatistic() == ParticleOnLattice::BosonicStatistic) && (this->ContactInteractionU!=0.0))
    {
      // general four-particle interactions:
      int Pos=0;
      int k4;
      double* TmpCoefficient = new double [Lx * p * KyMax * KyMax * KyMax];
      double Strength = this->ContactInteractionU / KyMax;  // Check Amplitude!!
      double MaxCoefficient = 0.0;
      for (int i=0; i<Lx; ++i)
	for (int s=0; s<p; ++s)
	  for (int k1 = 0; k1 < this->KyMax; ++k1)
	    for (int k2 = 0; k2 <= k1; ++k2)
	      for (int k3 = 0; k3 < this->KyMax; ++k3)
		{
		  k4 = k1 + k2 - k3;
		  if (k4 < 0)
		    k4 += this->KyMax;
		  else
		    if (k4 >= this->KyMax)
		      k4 -= this->KyMax;
		  if (k3 > k4)
		    {
		      if (k1 != k2)
			{
			  TmpCoefficient[Pos] = 4.0*Strength;
			}
		      else
			TmpCoefficient[Pos] = 2.0*Strength;
		      if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
			MaxCoefficient = fabs(TmpCoefficient[Pos]);
		      ++Pos;
		    }
		  else
		    if (k3 == k4)
		      {
			if (k1 != k2)
			  TmpCoefficient[Pos] = 2.0*Strength;
			else
			  TmpCoefficient[Pos] = Strength;
			if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
			  MaxCoefficient = fabs(TmpCoefficient[Pos]);
			++Pos;
		      }
		}
      this->NbrInteractionFactors = 0;
      this->Q1Value = new int [Pos];
      this->Q2Value = new int [Pos];
      this->Q3Value = new int [Pos];
      this->Q4Value = new int [Pos];
      this->InteractionFactors = new Complex [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i=0; i<Lx; ++i)
	for (int s=0; s<p; ++s)
	  for (int k1 = 0; k1 < this->KyMax; ++k1)
	    for (int k2 = 0; k2 <= k1; ++k2)
	      for (int k3 = 0; k3 < this->KyMax; ++k3)
		{
		  k4 = k1 + k2 - k3;
		  if (k4 < 0)
		    k4 += this->KyMax;
		  else
		    if (k4 >= this->KyMax)
		      k4 -= this->KyMax;
		  if (k3 >= k4)
		    {
		      if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
			{
			  this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
			  this->Q1Value[this->NbrInteractionFactors] = Particles->EncodeQuantumNumber(i, k1*p+s, 0, TranslationPhase);
			  this->Q2Value[this->NbrInteractionFactors] = Particles->EncodeQuantumNumber(i, k2*p+s, 0, TranslationPhase);
			  this->Q3Value[this->NbrInteractionFactors] = Particles->EncodeQuantumNumber(i, k3*p+s, 0, TranslationPhase);
			  this->Q4Value[this->NbrInteractionFactors] = Particles->EncodeQuantumNumber(i, k4*p+s, 0, TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "4-body term: "<<this->InteractionFactors[this->NbrInteractionFactors]<< " a^\\dag_"<< this->Q4Value[this->NbrInteractionFactors]<<"a^\\dag_"<<this->Q3Value[this->NbrInteractionFactors]<<"a_"<<this->Q2Value[this->NbrInteractionFactors]<<"a_"<<this->Q1Value[this->NbrInteractionFactors]<<endl;
#endif
			  ++this->NbrInteractionFactors;
			}
		      ++Pos;
		    }
		}
      delete [] TmpCoefficient;
    }
  else // no such interactions
    {
      this->NbrInteractionFactors = 0;
    }

  
    
  
}
