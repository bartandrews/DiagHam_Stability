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
#include "Hamiltonian/ParticleOnLatticeKapitMuellerMultiLayerHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
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
// lx = length of simulation cell in x-direction
// ly = length of simulation cell in y-direction
// nbrFluxQuanta = number of flux quanta piercing the simulation cell
// contactInteractionU = strength of on-site delta interaction
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
ParticleOnLatticeKapitMuellerMultiLayerHamiltonian::ParticleOnLatticeKapitMuellerMultiLayerHamiltonian(ParticleOnLattice* particles, int nbrParticles, int lx, int ly, int nbrFluxQuanta, double contactInteractionU, double contactInteractionW, bool reverseHopping, double deltaPotential, double randomPotential, double range, int nbrBranchCuts, double *branchCoordinates, int *branchShift, AbstractArchitecture* architecture, int nbrBody, unsigned long memory, char* precalculationFileName, bool hermitianFlag)
{
  this->Particles=particles;
  this->NbrParticles=nbrParticles;
  this->Lx=lx;
  this->Ly=ly;
  this->NbrLayers=this->Particles->GetNbrSublattices();
  this->NbrSublattices=1;
  this->HaveKySymmetry=false;
  this->KyMax=0;  
  this->NbrCells=lx*ly;
  this->NbrSites=NbrCells*NbrLayers;
  this->NbrFluxQuanta=nbrFluxQuanta;
  this->HamiltonianShift=0.0;
  this->FluxDensity=((double)nbrFluxQuanta)/NbrCells;
  this->ContactInteractionU=contactInteractionU;
  this->ContactInteractionW=contactInteractionW;
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


  // initialize branch cuts
  int NbrImagesOfCutsEachWay = 10;
  this->NbrBranchCuts = (2*NbrImagesOfCutsEachWay+1)*(2*NbrImagesOfCutsEachWay+1)*nbrBranchCuts;

  // array with information about branch cuts
  int index=0;
  if (this->NbrBranchCuts>0)
    {
      this->BranchCuts = new KMBranchCut[this->NbrBranchCuts];
      this->BranchShift = new int[this->NbrBranchCuts];
      
      double coords[4];
      for (int i=0; i<nbrBranchCuts; ++i) // loop over number of inputs
	{
	  // create images of the branch cuts (explicitly) in the nearby simulation cells only
	  for (int dX = -NbrImagesOfCutsEachWay; dX <=NbrImagesOfCutsEachWay; ++dX)
	    for (int dY = -NbrImagesOfCutsEachWay; dY <=NbrImagesOfCutsEachWay; ++dY)
	      {
		coords[0] = branchCoordinates[4*i] + dX*this->Lx;
		coords[1] = branchCoordinates[4*i+1] + dY*this->Ly;
		coords[2] = branchCoordinates[4*i+2] + dX*this->Lx;
		coords[3] = branchCoordinates[4*i+3] + dY*this->Ly;

		this->BranchCuts[index]=KMBranchCut(coords);
		if (branchShift!=NULL)
		  this->BranchShift[index]=branchShift[i];
		else
		  this->BranchShift[index]=1;

		// cout << "cut "<<index<<": "<<this->BranchCuts[index]<<" with shift "<<this->BranchShift[index]<<endl;
		++index;
	      }
	}
    }
  assert(index == this->NbrBranchCuts);

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
ParticleOnLatticeKapitMuellerMultiLayerHamiltonian::~ParticleOnLatticeKapitMuellerMultiLayerHamiltonian()
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
  if (this->NbrBranchCuts>0)
    {
      delete [] this->BranchCuts;
      delete [] this->BranchShift;
    }
}


// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream
ostream& operator << (ostream& Str, ParticleOnLatticeKapitMuellerMultiLayerHamiltonian& H)
{
  Str << "Need to implement ostream& operator << for ParticleOnLatticeKapitMuellerMultiLayerHamiltonian!" << endl;
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream
MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeKapitMuellerMultiLayerHamiltonian& H)
{
  Str << "Need to implement MathematicaOutput& operator << for ParticleOnLatticeKapitMuellerMultiLayerHamiltonian!\n";
  return Str;
}


// evaluate all interaction factors
//   
void ParticleOnLatticeKapitMuellerMultiLayerHamiltonian::EvaluateInteractionFactors()
{  
  // hopping terms are present independent of statistics:
  this->NbrHoppingTerms=this->NbrSites*this->NbrSites;
  cout << "this->DeltaPotential="<<this->DeltaPotential<<", this->RandomPotential="<<this->RandomPotential<<endl;
  cout << "this->NbrSites"<<this->NbrSites<<endl;
  cout << "NbrHoppingTerms="<<NbrHoppingTerms<<endl;
  if (this->DeltaPotential != 0.0) ++this->NbrHoppingTerms;
  if (this->RandomPotential != 0.0) this->NbrHoppingTerms+=this->NbrSites;
  this->HoppingTerms = new Complex[NbrHoppingTerms];
  this->KineticQi = new int[NbrHoppingTerms];
  this->KineticQf = new int[NbrHoppingTerms];

  int TmpNumberTerms=0;
  double HoppingSign = (this->ReverseHopping ? 1.0 : -1.0);
  Complex TranslationPhase;

  // manual parameters, for now
  int images = 50; // number of images of the simulation cell
  int maxRange = std::sqrt(30./(1.0-this->FluxDensity)); // larger exponents will yield numerical zero
  if (this->Range > maxRange)
    {
      // cout << "Testing: range not reduced to new maxRange="<<maxRange<<endl;
      this->Range = maxRange;
    }
  images = this->Range / (Lx < Ly? Lx:Ly) + 2;
  Complex amplitude;
  //loop over initial sites (i,j)
  
  int *qi = new int[this->NbrLayers];
  int *qf = new int[this->NbrLayers];

  for (int i=0; i<Lx; ++i)
    {
      for (int j=0; j<Ly; ++j)
	{
	  for (int s=0; s<this->NbrLayers; ++s)
	    qi[s] = Particles->EncodeQuantumNumber(i, j, s, TranslationPhase);
	      	      
	  // have long-range hopping, so sum over all possible final sites (k,l)
	  for (int k=0; k<Lx; ++k)
	    {
	      for (int l=0; l<Ly; ++l)
		{		  
		  for (int s=0; s<this->NbrLayers; ++s)
		    qf[s] = Particles->EncodeQuantumNumber(k, l, s, TranslationPhase);
		  
		  
		  ComplexMatrix sumHopping(this->NbrLayers, this->NbrLayers, true);
		  for (int dX = -images; dX <= images; ++dX)
		    for (int dY = -images; dY <= images; ++dY)
		      {
			if ( this->GetDistance(k+dX*Lx - i, l+dY*Ly - j) <= this->Range // allow hard cut-off
			     // && ( ( dX!=0 || dY!=0 ) || ( i!=k || j!=l ) ) // optionally, exclude onsite terms (needed to get zero energy for lowest band)
			     )
			  {
			    Particles->EncodeQuantumNumber(k+dX*Lx, l+dY*Ly, 0, TranslationPhase); // just need to get new TranslationPhase, here
			    amplitude = -HoppingSign*this->KapitMuellerHopping(k+dX*Lx, l+dY*Ly, i, j) * Conj(TranslationPhase);

			    if (Norm(amplitude) > 1e-13)
			      {
				// check for any branch cuts being crossed
				int shiftModulo = this->EvaluateBranchCrossings(i, j, k+dX*Lx, l+dY*Ly);
				while (shiftModulo<0) shiftModulo+=this->NbrLayers;
				shiftModulo %= this->NbrLayers;		    

				for (int si=0; si<this->NbrLayers; ++si)
				  {
				    int sf = (si+shiftModulo);
				    while (sf<0) sf+=this->NbrLayers;
				    sf %= this->NbrLayers;		    
				    sumHopping.AddToMatrixElement(si, sf, amplitude);
				  }
#ifdef DEBUG_OUTPUT
				//if (TranslationPhase!=1.0)
				cout << "xx ("<<i<<", "<<j<<")->("<<k+dX*Lx<<", "<<l+dY*Ly<<") amplitude="<<amplitude<<", shift="<< shiftModulo <<", dL=("<<dX<<"," <<dY<<") : "
				     <<"Translation ["<<qi[0]<<"->"<<qf[0]<<"]="
				     << TranslationPhase << endl;
#endif
			      }
			  }
		      }
#ifdef DEBUG_OUTPUT		     
		  //if (Norm(sumHopping[0][0])>1e-15)
		  {
		    cout << "H[("<<i<<", "<<j<<")->("<<k<<", "<<l<<")]= ";
		    for (int si=0; si<this->NbrLayers; ++si)
		      for (int sf=0; sf<this->NbrLayers; ++sf)
			cout << sumHopping[si][sf]<<", ";
		    cout<<" tP="<<TranslationPhase<<endl;
		  }
#endif
		  for (int si=0; si<this->NbrLayers; ++si)
		    for (int sf=0; sf<this->NbrLayers; ++sf)
		      {
			// only take into account terms with non-zero magnitude
			sumHopping.GetMatrixElement(si, sf, amplitude);
			if (Norm(amplitude)>1e-15)
			  {			    
			    KineticQi[TmpNumberTerms] = qi[si];
			    KineticQf[TmpNumberTerms] = qf[sf];
			    HoppingTerms[TmpNumberTerms] = amplitude;
			    ++TmpNumberTerms;
			  }
		      }
		}
	    }
	}
    }
  delete [] qi;
  delete [] qf;

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
  
  // inter-layer interactions
  if (this->NbrLayers>0 && this->ContactInteractionW!=0.0)
    {
      Complex translationPhase;
      this->NbrRhoRhoInteractionFactors = this->NbrCells;
      this->RhoRhoQ12Values = new int[2*this->NbrRhoRhoInteractionFactors];
      this->RhoRhoInteractionFactors = new double[this->NbrRhoRhoInteractionFactors];
      cout << "adding interlayer interaction terms"<<endl;
      int index=0;
      for (int i=0; i<this->Lx; ++i)
	for (int j=0; j<this->Ly; ++j)
	  {
	    for (int s=0; s<this->NbrLayers; ++s)
	      {
		int q1 = this->Particles->EncodeQuantumNumber(i, j, s, translationPhase);
		for (int s2=s+1; s2<this->NbrLayers; ++s2)
		  {
		    int q2 = this->Particles->EncodeQuantumNumber(i, j, s2, translationPhase);
		    this->RhoRhoQ12Values[index<<1]=q1;
		    this->RhoRhoQ12Values[(index<<1)+1]=q2;
		    this->RhoRhoInteractionFactors[index++]=this->ContactInteractionW;
		  }
	      }
	  }
    }
  else // no such interactions
    {
      this->NbrRhoRhoInteractionFactors=0;
    }
}

/// evaluate branch crossings for a link from site i to j
/// @param xi, yi, xj, yj coordinates of the initial (i) and final site (j)
/// @return overall shift for crossing any of the known branch cuts
int ParticleOnLatticeKapitMuellerMultiLayerHamiltonian::EvaluateBranchCrossings(int xi, int yi, int xj, int yj)
{
  RealVector Ri(2), Rj(2);
  Ri[0]=xi;
  Ri[1]=yi;
  Rj[0]=xj;
  Rj[1]=yj;
  RealVector Rij(Rj, true);
  Rij-=Ri;
  double xiA, xiB, xiJ, xiS, etaA, etaB, etaJ;
  int shift=0;
  for (int i=0; i<NbrBranchCuts; ++i)
    {
      RealVector RA(this->BranchCuts[i].A, true);
      RA -= Ri;
      RealVector RB(this->BranchCuts[i].B, true);
      RB -= Ri;
      xiA=this->BranchCuts[i].E1 * RA;
      etaA=this->BranchCuts[i].E2 * RA;
      xiB=this->BranchCuts[i].E1 * RB;
      etaB=this->BranchCuts[i].E2 * RB;
      xiJ=this->BranchCuts[i].E1 * Rij;
      etaJ=this->BranchCuts[i].E2 * Rij;
      if (fabs(etaJ)>1e-13) // Rij not parallel to cut?
	{
	  xiS = xiJ * etaA / etaJ;
	  double etaS = etaA;
	  // crossing point in interval of branch cut?
	  if ( (xiS >= fmin(xiA,xiB)-1e-13) && (xiS <= fmax(xiA,xiB)+1e-13)
	       && (xiS <= fmax(0.0,xiJ)+1e-13) && (xiS >= fmin(0.0,xiJ)-1e-13)
	       && ( etaS <=  fmax(0.0,etaJ)+1e-13) && (etaS >= fmin(0.0,etaJ)-1e-13) )
	  // if ( ( (xiS >= xiA && xiS <= xiB) || (xiS <= xiA && xiS >= xiB) )
	  //      && ( (xiS >= 0 && xiS <= xiJ) || (xiS <= 0 && xiS >= xiJ) )
	  //      && ( (etaS >= 0 && etaS <= etaJ) || (etaS <= 0 && etaS >= etaJ) ) )
	    {
	      shift += this->BranchShift[i] * (1 - 2 * (etaJ < 0));
	      // RealVector S(Ri, true);
	      // S.AddLinearCombination(etaS/etaJ, Rij);
	      // cout << xi<<yi<<xj<<yj<<" ("<<xi<<", "<<yi<<")->("<<xj<<", "<<yj<<") crossing "<<this->BranchCuts[i]<<" at ("<<S[0]<<", "<<S[1]<<")"<<endl;
	      // cout << "RA="<<RA<<"RB="<<RB<<"E1="<<this->BranchCuts[i].E1<<"E2="<<this->BranchCuts[i].E2<<" xiA="<<xiA<<" etaA="<<etaA<<" xiB="<<xiB<<" etaB="<<etaB<<" xiJ="<<xiJ<<" etaJ="<<etaJ<<", xiS="<<xiS<<endl;
	    }
	  else
	    {
	      // cout << xi<<yi<<xj<<yj<<" ("<<xi<<", "<<yi<<")->("<<xj<<", "<<yj<<") not crossing"<<this->BranchCuts[i]<<" "<<(xiS >= fmin(xiA,xiB)) <<" "<< (xiS <= fmax(xiA,xiB)) <<" "<< (xiS <= fmax(0.0,xiJ)) <<" "<< (xiS >= fmin(0.0,xiJ)) <<" "<< ( etaS <=  fmax(0.0,etaJ)) <<" "<< (etaS >= fmin(0.0,etaJ)) << " xiS="<<xiS<<" xiJ="<<xiJ<<endl;
	      // cout << "RA="<<RA<<"RB="<<RB<<"Ri="<<Ri<<"E1="<<this->BranchCuts[i].E1<<"E2="<<this->BranchCuts[i].E2<<" xiA="<<xiA<<" etaA="<<etaA<<" xiB="<<xiB<<" etaB="<<etaB<<" xiJ="<<xiJ<<" etaJ="<<etaJ<<", xiS="<<xiS<<endl;
	    }
	}
    }
  return shift;
}
