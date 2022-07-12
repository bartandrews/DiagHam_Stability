////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2016 Gunnar Moeller                    //
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
#include "ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian.h"
#include "ParticleOnLatticeKapitMuellerHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeKapitMuellerMultiLayerHamiltonian.h"
#include "HilbertSpace/SingleParticleOnLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/StringTools.h"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
using std::cout;
using std::endl;

using std::ostream;

// switch for debugging output:
//#define DEBUG_OUTPUT


// class for storing information about branch cuts
class KMBranchCut
{
public:
  // initial point
  RealVector A;
  // final point
  RealVector B;
  // unit vector in direction AB
  RealVector E1; 
  // unit vector in direction e_z ^ AB 
  RealVector E2;

  // default constructor
  KMBranchCut():
    A(2, true),
    B(2, true),
    E1(2,true),
    E2(2,true)
  {}

  // constructor
  /// @param r coordinates in order of (x_A, y_A, x_B, y_B)
  /// @param shift modulo added to layer index when crossing branch cut
  KMBranchCut(double *r, int shift=1):
    A(2),
    B(2),
    E2(2,true)
  {
    A[0]=r[0];    A[1]=r[1];
    B[0]=r[2];    B[1]=r[3];
    E1 = B - A;
    E1.Normalize();
    E2[0]=-E1[1];
    E2[1]=E1[0];
  }
  
  // copy constructor
  KMBranchCut(const KMBranchCut& cut):
    A(cut.A, true),
    B(cut.B, true),
    E1(cut.E1, true),
    E2(cut.E2, true)
  {}

  // destructor
  ~KMBranchCut()
  {
  }

  // Output Stream overload
  //
  // str = reference on output stream
  // cut = branch cut to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, KMBranchCut& c)
  {
    str << "("<<c.A[0]<<", "<<c.A[1]<<")->("<<c.B[0]<<", "<<c.B[1]<<")";
    return str;
  }

};




// constructor for contact interactions on a square lattice
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// numStates = number of single particle states to consider (projection to lowest numState SP states)
// flatBand = apply the flat-band projection for all kept single-particle states
// lx = length of simulation cell in x-direction
// ly = length of simulation cell in y-direction
// nbrLayers = number of layers to simulate
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// contactInteractionU = strength of on-site delta interaction
// contactInteractionW = strength of inter-layer on-site delta interaction
// reverseHopping = flag to indicate if sign of hopping terms should be reversed
// deltaPotential = strength of a delta potential at site (0,0)
// randomPotential = magnitude of random potential to add to all sites
// range = maximum length cut-off for single-particle hoppings
// nbrBranchCuts = number of branch cuts
// branchCoordinates coordinates of branch cuts in successive blocks of order (x_A, y_A, x_B, y_B)_i, i=1..nbrBranchCuts
// branchShift shift of the layer index at the i-th branch cut
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// hermitianFlag = flag indicating whether to use hermitian symmetry
ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian::ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian(ParticleOnLattice* particles, int nbrParticles, int numStates, bool flatBand, int lx, int ly, int nbrLayers, int nbrFluxQuanta, double contactInteractionU, double contactInteractionW, bool reverseHopping, double deltaPotential, double randomPotential, double range, int nbrBranchCuts, double *branchCoordinates, int *branchShift, AbstractArchitecture* architecture, int nbrBody, unsigned long memory, char* precalculationFileName, bool hermitianFlag)
{
  this->Particles=particles;
  this->NbrParticles=nbrParticles;
  this->NbrProjectorStates = numStates;
  this->Lx=lx;
  this->Ly=ly;
  this->NbrLayers=nbrLayers;
  this->HaveKySymmetry=false;
  this->KyMax=0;
  this->NbrCells=lx*ly;
  this->NbrSites=NbrCells*NbrLayers;
  if (this->NbrProjectorStates>this->NbrSites)
    this->NbrProjectorStates = this->NbrSites;
  this->NbrFluxQuanta=nbrFluxQuanta;
  this->HamiltonianShift=0.0;
  this->FluxDensity=((double)nbrFluxQuanta)/NbrCells;
  this->ContactInteractionU=contactInteractionU;
  this->ContactInteractionW=contactInteractionW;
  this->DeltaPotential = deltaPotential;
  this->RandomPotential = randomPotential;
  this->Range = range;
  this->NbrBranchCuts = nbrBranchCuts;
  this->BranchCoordinates = branchCoordinates;
  this->BranchShift = branchShift;
  this->ReverseHopping = reverseHopping;
  this->FlatBand = flatBand;
  this->EnergyLevels=NULL;

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
	  cout << "dimension="<<this->GetHilbertSpaceDimension()<<endl;
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
ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian::~ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian()
{
  if (NbrHoppingTerms>0)
    {
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
    }
  if (this->EnergyLevels!=NULL)
    delete [] this->EnergyLevels;
  if (NbrInteractionFactors>0 && NbrQ12Indices==0)
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
      delete [] this->Q3PerQ12;
      delete [] this->Q4PerQ12;
      delete [] this->NbrQ34Values;
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      
    }
  if (this->NbrDiagonalInteractionFactors>0)
    {
      delete [] this->DiagonalInteractionFactors;
      delete [] this->DiagonalQValues;
    }
  if (this->NbrRhoRhoInteractionFactors>0)
    {
      delete [] this->RhoRhoQ12Values;
      delete [] this->RhoRhoInteractionFactors;
    }
}


// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream
ostream& operator << (ostream& Str, ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian& H)
{
  Str << "Need to implement ostream& operator << for ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian!" << endl;
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream
MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian& H)
{
  Str << "Need to implement MathematicaOutput& operator << for ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian!\n";
  return Str;
}


// evaluate all interaction factors
//   
void ParticleOnLatticeProjectedKapitMuellerMultiLayerHamiltonian::EvaluateInteractionFactors()
{  
  // create single-particle space and Hamiltonian
  double solenoidX, solenoidY;
  this->Particles->GetSolenoidFluxes(solenoidX, solenoidY);
  SingleParticleOnLattice SingleParticleSpace(this->NbrParticles, this->Lx, this->Ly, this->NbrFluxQuanta, 10000000, solenoidX, solenoidY, this->Particles->GetLandauGaugeAxis(), this->NbrLayers);

  AbstractHamiltonian *SingleParticleHamiltonian;
    
  if (NbrLayers>1)
     SingleParticleHamiltonian = new ParticleOnLatticeKapitMuellerMultiLayerHamiltonian(&SingleParticleSpace, 1, this->Lx, this->Ly, this->NbrFluxQuanta, 0.0, 0.0, this->ReverseHopping, this->DeltaPotential, this->RandomPotential, this->Range, this->NbrBranchCuts, this->BranchCoordinates, this->BranchShift, this->Architecture, 2, /* memory = */ 0, /* precalculationFileName = */ 0, /* hermitianFlag = */ false);
  else
    SingleParticleHamiltonian = new ParticleOnLatticeKapitMuellerHamiltonian(&SingleParticleSpace, 1, this->Lx, this->Ly, this->NbrFluxQuanta, 0.0, this->ReverseHopping, this->DeltaPotential, this->RandomPotential, this->Range, this->Architecture, 2, /* memory = */ 0, /* precalculationFileName = */ 0, /* hermitianFlag = */ false);
  
  HermitianMatrix TmpOneBodyHamiltonian(this->NbrSites, true);
  SingleParticleHamiltonian->GetHamiltonian(TmpOneBodyHamiltonian);
  
  // output for testing
  if (SingleParticleHamiltonian->GetHilbertSpaceDimension() < 10)
    {
      cout << "Single-Particle Hamiltonian:\n"<<TmpOneBodyHamiltonian;
    }
  
  this->OneBodyBasis.Resize(this->NbrSites, this->NbrSites);
  this->OneBodyBasis.SetToIdentity();
  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, this->OneBodyBasis);
#else
  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, this->OneBodyBasis);
#endif
  this->EnergyLevels = new double[this->NbrSites];
  for (int i = 0; i < this->NbrSites; ++i)
    {
      this->EnergyLevels[i] = TmpDiag(i, i);
      cout << "E["<<i<<"]="<<this->EnergyLevels[i]<<endl;
    }

  // output for testing
  if (SingleParticleHamiltonian->GetHilbertSpaceDimension() < 10)
    for (int i = 0; i < this->NbrSites; ++i)
      {
	this->EnergyLevels[i] = TmpDiag(i, i);
	cout << "Phi["<<i<<"]="<<this->OneBodyBasis[i]<<endl;
      }

  delete SingleParticleHamiltonian;


  // hopping terms are present independent of statistics:
  this->NbrHoppingTerms=0;
  if (this->FlatBand == false)
    this->NbrHoppingTerms+=this->NbrProjectorStates;
  if ((this->DeltaPotential != 0.0) || (this->RandomPotential != 0.0)) 
    this->NbrHoppingTerms+=this->NbrProjectorStates * this->NbrProjectorStates;  
  this->HoppingTerms = new Complex[NbrHoppingTerms];
  this->KineticQi = new int[NbrHoppingTerms];
  this->KineticQf = new int[NbrHoppingTerms];

  Complex TranslationPhase;
  int TmpNumberTerms=0;
  if (this->FlatBand == false)
    for (int n=0; n<this->NbrProjectorStates; ++n)
      {
	KineticQi[TmpNumberTerms] = n;
	KineticQf[TmpNumberTerms] = n;
	HoppingTerms[TmpNumberTerms++] = EnergyLevels[n];
      }

  if (this->DeltaPotential != 0.0 || this->RandomPotential != 0.0)
    {
      double LocalPotentials[this->NbrSites];
      for (int i=0; i<this->NbrSites; ++i) LocalPotentials[i] = 0.0;
      for (int l=0; l<this->NbrLayers; ++l)
	{
	  int q = SingleParticleSpace.EncodeQuantumNumber(0, 0, l, TranslationPhase);
	  LocalPotentials[q] += this->DeltaPotential;
	}
      if (this->RandomPotential != 0.0)
	{
	  NumRecRandomGenerator G;
	  G.UseTimeSeed();
	  for (int x=0; x<this->Lx; ++x)
	    for (int y=0; y<this->Ly; ++y)
	      for (int l=0; l<this->NbrLayers; ++l)
		{
		  int q = SingleParticleSpace.EncodeQuantumNumber(x, y, l, TranslationPhase);
		  LocalPotentials[q] += this->RandomPotential*(-0.5+G.GetRealRandomNumber());
		}
	}
      Complex TmpHopping;
      for (int n1=0; n1<this->NbrProjectorStates; ++n1)
	for (int n2=0; n2<this->NbrProjectorStates; ++n2)
	  {
	    KineticQi[TmpNumberTerms] = n1;
	    KineticQf[TmpNumberTerms] = n2;
	    TmpHopping=0.0;
	    for (int q=0; q<this->NbrSites; ++q)
	      TmpHopping += LocalPotentials[q] * Conj(this->OneBodyBasis[n1][q]) * this->OneBodyBasis[n2][q];
	    HoppingTerms[TmpNumberTerms++] = TmpHopping;
	  }
    }

  bool symmetries = true;
  if (symmetries == false)
    {
      // general two-particle interactions:
      // first attempt: no optimization for symmetries under pairs of creation / annihilation operators

      if (this->Particles->GetParticleStatistic() == ParticleOnLattice::FermionicStatistic)
	this->NbrQ12Indices=this->NbrProjectorStates * (this->NbrProjectorStates-1);
      else
	this->NbrQ12Indices=this->NbrProjectorStates * this->NbrProjectorStates;
  
      this->Q1Value = new int[this->NbrQ12Indices];
      this->Q2Value = new int[this->NbrQ12Indices];
      this->NbrQ34Values = new int[this->NbrQ12Indices];
      this->Q3PerQ12 = new int*[this->NbrQ12Indices];
      this->Q4PerQ12 = new int*[this->NbrQ12Indices];
      for (int i=0; i<this->NbrQ12Indices; ++i)
	{
	  this->NbrQ34Values[i] = this->NbrQ12Indices;   
	  this->Q3PerQ12[i] = new int[this->NbrQ12Indices];
	  this->Q4PerQ12[i] = new int[this->NbrQ12Indices];
	}
      this->NbrInteractionFactors=this->NbrQ12Indices*this->NbrQ12Indices;
      this->InteractionFactors = new Complex[this->NbrInteractionFactors];
  
      int tmpNbrInteractionFactors = 0;
      int tmpQ12Index = 0;
      bool fermions = this->Particles->GetParticleStatistic() == ParticleOnLattice::FermionicStatistic;
      Complex tmpSum;
      cout << "adding interaction terms"<<endl;
		  
      for (int n1=0; n1<this->NbrProjectorStates; ++n1)
	for (int n2=0; n2<this->NbrProjectorStates; ++n2)
	  {
	    if (fermions && n1==n2) {
	      if (n2 == this->NbrProjectorStates-1) 
		break;
	      else 
		++n2;
	    }
	    this->Q1Value[tmpQ12Index] = n1;
	    this->Q2Value[tmpQ12Index] = n2;
	    int tmpQ34Index = 0;
	    int *TmpQ3Values = this->Q3PerQ12[tmpQ12Index];
	    int *TmpQ4Values = this->Q4PerQ12[tmpQ12Index];
	    for (int n3=0; n3<this->NbrProjectorStates; ++n3)
	      for (int n4=0; n4<this->NbrProjectorStates; ++n4)
		{
		  if (fermions && n3==n4) {
		    if (n4 == this->NbrProjectorStates-1) 
		      break;
		    else 
		      ++n4;
		  }
		  TmpQ3Values[tmpQ34Index] = n3;
		  TmpQ4Values[tmpQ34Index] = n4;
	      
		  tmpSum = 0.0;

		  // intra-layer interactions
		  if (this->ContactInteractionU!=0.0)
		    {
		      for (int i=0; i<this->Lx; ++i)
			for (int j=0; j<this->Ly; ++j)
			  for (int s=0; s<this->NbrLayers; ++s)
			    {
			      int q = SingleParticleSpace.EncodeQuantumNumber(i, j, s, TranslationPhase); 
			      // if (q>=this->NbrSites)
			      //   cout << "Error: exceed number of single body states at i="<<i<<", j"<<j<<", s="<<s<<" -> q="<<q<<" (NbrSites="<<NbrSites<<")"<<endl;
			      tmpSum += this->ContactInteractionU * Conj(this->OneBodyBasis[n3][q] * this->OneBodyBasis[n4][q])
				* this->OneBodyBasis[n1][q] * this->OneBodyBasis[n2][q];
			      // if (n1==0 && n2==0 && n3==0 && n4==0)
			      //   cout << "Contribution for n=(0,0,0,0) at q="<<q<<" for (i,j,s)="<<i<<", "<<j<<", "<<s<<") is = "<< Conj(this->OneBodyBasis[n3][q] * this->OneBodyBasis[n4][q])
			      //     * this->OneBodyBasis[n1][q] * this->OneBodyBasis[n2][q] << endl;
			    }
		    }

		  // inter-layer interactions
		  if (this->NbrLayers>1 && this->ContactInteractionW!=0.0)
		    {
		      for (int i=0; i<this->Lx; ++i)
			for (int j=0; j<this->Ly; ++j)
			  for (int s=0; s<this->NbrLayers; ++s)
			    {
			      int q1 = SingleParticleSpace.EncodeQuantumNumber(i, j, s, TranslationPhase);
			      for (int s2=s+1; s2<this->NbrLayers; ++s2)
				{
				  int q2 = SingleParticleSpace.EncodeQuantumNumber(i, j, s2, TranslationPhase);
				  tmpSum += this->ContactInteractionW * Conj(this->OneBodyBasis[n4][q1] * this->OneBodyBasis[n3][q2])
				    * this->OneBodyBasis[n2][q2] * this->OneBodyBasis[n1][q1];
				}
			    }
		    }
		  // cout << this->Q1Value[tmpQ12Index] << " " << this->Q2Value[tmpQ12Index] << " " << TmpQ3Values[tmpQ34Index] << " " << TmpQ4Values[tmpQ34Index] << " " << tmpSum << endl;
	      
		  if (Norm(tmpSum) > 10. * MACHINE_PRECISION)  // check for numerically non-zero coefficients
		    this->InteractionFactors[tmpNbrInteractionFactors++] = tmpSum;
		  else
		    this->InteractionFactors[tmpNbrInteractionFactors++] = 0.0; // waste some time by keeping these coefficients, for now.
		  ++tmpQ34Index;

		}
	    ++tmpQ12Index;
	  }

      if (tmpQ12Index < NbrQ12Indices)
	{
	  cout << "Inconsistent count of NbrQ12Indices!"<<endl;
	  exit(1);
	}
        if (this->NbrInteractionFactors > tmpNbrInteractionFactors)
	  {
	    Complex *newIFs = new Complex[tmpNbrInteractionFactors];
	    for (int i=0; i<tmpNbrInteractionFactors; ++i)
	      newIFs[i] = this->InteractionFactors[i];
	    delete [] this->InteractionFactors;
	    this->InteractionFactors = newIFs;
	  }
	this->NbrInteractionFactors = tmpNbrInteractionFactors;
    }
  else
    { // matrix elements with symmetries considered:

      if (this->Particles->GetParticleStatistic() == ParticleOnLattice::FermionicStatistic)
	this->NbrQ12Indices=this->NbrProjectorStates * (this->NbrProjectorStates-1) / 2;
      else
	this->NbrQ12Indices=this->NbrProjectorStates * (this->NbrProjectorStates +1) / 2;
  
      this->Q1Value = new int[this->NbrQ12Indices];
      this->Q2Value = new int[this->NbrQ12Indices];
      this->NbrQ34Values = new int[this->NbrQ12Indices];
      this->Q3PerQ12 = new int*[this->NbrQ12Indices];
      this->Q4PerQ12 = new int*[this->NbrQ12Indices];
      for (int i=0; i<this->NbrQ12Indices; ++i)
	{
	  this->NbrQ34Values[i] = this->NbrQ12Indices;   
	  this->Q3PerQ12[i] = new int[this->NbrQ12Indices];
	  this->Q4PerQ12[i] = new int[this->NbrQ12Indices];
	}
      this->NbrInteractionFactors=this->NbrQ12Indices*this->NbrQ12Indices;
      this->InteractionFactors = new Complex[this->NbrInteractionFactors];
  
      int tmpNbrInteractionFactors = 0;
      int tmpQ12Index = 0;
      bool fermions = this->Particles->GetParticleStatistic() == ParticleOnLattice::FermionicStatistic;
      Complex tmpSum;
      cout << "adding interaction terms"<<endl;

      int fermionOffset = (fermions ? 1 : 0);
      int fermionSign = 1-2*fermions;
      for (int n1=0; n1<this->NbrProjectorStates; ++n1)
	for (int n2=n1+fermionOffset; n2<this->NbrProjectorStates; ++n2)
	  {
	    this->Q1Value[tmpQ12Index] = n1;
	    this->Q2Value[tmpQ12Index] = n2;
	    int tmpQ34Index = 0;
	    int *TmpQ3Values = this->Q3PerQ12[tmpQ12Index];
	    int *TmpQ4Values = this->Q4PerQ12[tmpQ12Index];
	    for (int n3=0; n3<this->NbrProjectorStates; ++n3)
	      for (int n4=n3+fermionOffset; n4<this->NbrProjectorStates; ++n4)
		{
		  TmpQ3Values[tmpQ34Index] = n3;
		  TmpQ4Values[tmpQ34Index] = n4;
		  tmpSum = 0.0;

		  // intra-layer interactions
		  if (!fermions && this->ContactInteractionU!=0.0)
		    {
		      for (int i=0; i<this->Lx; ++i)
			for (int j=0; j<this->Ly; ++j)
			  for (int s=0; s<this->NbrLayers; ++s)
			    {
			      int q = SingleParticleSpace.EncodeQuantumNumber(i, j, s, TranslationPhase); 
			      // if (q>=this->NbrSites)
			      //   cout << "Error: exceed number of single body states at i="<<i<<", j"<<j<<", s="<<s<<" -> q="<<q<<" (NbrSites="<<NbrSites<<")"<<endl;
			      // for bosons, this term is repeated four times. for fermions it is zero by Pauli exclusion!
			      tmpSum += 4.0 * this->ContactInteractionU * Conj(this->OneBodyBasis[n3][q] * this->OneBodyBasis[n4][q]) * this->OneBodyBasis[n1][q] * this->OneBodyBasis[n2][q];
			      // if (n1==0 && n2==0 && n3==0 && n4==0)
			      //   cout << "Contribution for n=(0,0,0,0) at q="<<q<<" for (i,j,s)="<<i<<", "<<j<<", "<<s<<") is = "<< Conj(this->OneBodyBasis[n3][q] * this->OneBodyBasis[n4][q])
			      //     * this->OneBodyBasis[n1][q] * this->OneBodyBasis[n2][q] << endl;
			    }
		    }

		  // inter-layer interactions
		  if (this->NbrLayers>1 && this->ContactInteractionW!=0.0)
		    {
		      for (int i=0; i<this->Lx; ++i)
			for (int j=0; j<this->Ly; ++j)
			  for (int s=0; s<this->NbrLayers; ++s)
			    {
			      int q1 = SingleParticleSpace.EncodeQuantumNumber(i, j, s, TranslationPhase);
			      for (int s2=s+1; s2<this->NbrLayers; ++s2)
				{
				  int q2 = SingleParticleSpace.EncodeQuantumNumber(i, j, s2, TranslationPhase);
				  tmpSum += this->ContactInteractionW * Conj(this->OneBodyBasis[n4][q1] * this->OneBodyBasis[n3][q2]) * this->OneBodyBasis[n2][q2] * this->OneBodyBasis[n1][q1];
				  tmpSum += fermionSign * this->ContactInteractionW * Conj(this->OneBodyBasis[n4][q1] * this->OneBodyBasis[n3][q2]) * this->OneBodyBasis[n1][q2] * this->OneBodyBasis[n2][q1];
				  tmpSum += fermionSign * this->ContactInteractionW * Conj(this->OneBodyBasis[n3][q1] * this->OneBodyBasis[n4][q2]) * this->OneBodyBasis[n2][q2] * this->OneBodyBasis[n1][q1];
				  tmpSum += this->ContactInteractionW * Conj(this->OneBodyBasis[n3][q1] * this->OneBodyBasis[n4][q2]) * this->OneBodyBasis[n1][q2] * this->OneBodyBasis[n2][q1];
				}
			    }
		    }

		  if (n1==n2) tmpSum *= 0.5;
		  if (n3==n4) tmpSum *= 0.5;

		  // cout << this->Q1Value[tmpQ12Index] << " " << this->Q2Value[tmpQ12Index] << " " << TmpQ3Values[tmpQ34Index] << " " << TmpQ4Values[tmpQ34Index] << " " << tmpSum << endl;
	      
		  if (Norm(tmpSum) > 10. * MACHINE_PRECISION)  // check for numerically non-zero coefficients
		    this->InteractionFactors[tmpNbrInteractionFactors++] = tmpSum;
		  else
		    this->InteractionFactors[tmpNbrInteractionFactors++] = 0.0; // waste some time by keeping these coefficients, for now.
		  ++tmpQ34Index;

		}
	    ++tmpQ12Index;
	  }
      if (tmpQ12Index < NbrQ12Indices)
	{
	  cout << "Inconsistent count of NbrQ12Indices!"<<endl;
	  exit(1);
	}
      if (this->NbrInteractionFactors > tmpNbrInteractionFactors)
	{
	  Complex *newIFs = new Complex[tmpNbrInteractionFactors];
	  for (int i=0; i<tmpNbrInteractionFactors; ++i)
	    newIFs[i] = this->InteractionFactors[i];
	  delete [] this->InteractionFactors;
	  this->InteractionFactors = newIFs;
	}
      this->NbrInteractionFactors = tmpNbrInteractionFactors;
    }




  // no explicit diagonal terms, here.
  this->NbrDiagonalInteractionFactors = 0;

  this->NbrRhoRhoInteractionFactors = 0;
}

