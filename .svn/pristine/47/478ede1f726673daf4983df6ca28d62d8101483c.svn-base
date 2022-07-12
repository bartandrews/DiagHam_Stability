////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Gunnar Möller                         //
//                                                                            //
//  class of tight binding model for the square lattice with homogeneous flux //
//                                                                            //
//                        last modification : 08/05/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelHofstadterFiniteCylinder.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>
using std::cout;
using std::endl;

// #define DEBUG_OUTPUT


// default constructor
//
// nbrCellsX = number of unit cells in the x direction
// nbrCellsY = number of unit cella in the y direction
// unitCellX = number of sites in unit cell in x direction
// unitCellY = number of sites in unit cell in y direction
// nbrFlux = number of flux quanta per unit cell
// axis = direction of Landau gauge within cell ('x' or 'y')
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
TightBindingModelHofstadterFiniteCylinder::TightBindingModelHofstadterFiniteCylinder(int nbrSiteX, int nbrSiteY, int nbrFlux, double * tunnelElementX, double * tunnelElementY, char axis,double gammaX, double gammaY,  AbstractArchitecture* architecture, double  fluxInserted, bool storeOneBodyMatrices, bool torusGeometry)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY =  nbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->TunnelElementX = tunnelElementX;
  this->TunnelElementY = tunnelElementY;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands =  this->NbrSiteY;
  this->NbrStatePerBand = this->NbrSiteX;
  this->FluxDensity =  ((double) nbrFlux)/ (this->NbrSiteX * this->NbrSiteY) ;
  this->Architecture = architecture;
  this->LxTranslationPhase = Polar(1.0, 2.0*M_PI*FluxDensity*this->NbrSiteX);
  this->FluxInserted = fluxInserted;
  this->TorusGeometry = torusGeometry;
  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
    }
  else
    {
      this->OneBodyBasis = 0;
    }
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    }

  this->ComputeBandStructure();  
}

// destructor
//

TightBindingModelHofstadterFiniteCylinder::~TightBindingModelHofstadterFiniteCylinder()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelHofstadterFiniteCylinder::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double K1;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
	  int Index = kx;
 	  if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	    {
	      K1 = this->KxFactor*(((double) kx) + this->GammaX);

	      // construct magnetic unit cell:
	      
	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      
     Complex TranslationPhase;
     int   PosY = 0;

      for (int i=0; i<1; ++i)
      {
	  int InitialIndex = this->EncodeSublatticeIndex(i, PosY, K1,TranslationPhase); // TranlationPhase always one, so can be discarded
	  int FinalIndex = this->EncodeSublatticeIndex(i+1, PosY, K1,TranslationPhase);
	  if (InitialIndex>=FinalIndex)
	    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -this->TunnelElementX[PosY]*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<< -1.0*TranslationPhase<<endl;
#endif
			  FinalIndex = this->EncodeSublatticeIndex(i-1, PosY, K1, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -this->TunnelElementX[PosY]*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<< -1.0*Conj(Phase)*TranslationPhase<<endl;
#endif		  
			  FinalIndex = this->EncodeSublatticeIndex(i, PosY+1, K1, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -this->TunnelElementY[PosY]*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif
  }
		  for (PosY = 1; PosY < NbrSiteY - 1; ++PosY)
		    {
                      Complex PhaseX = Phase(2.0*M_PI*this->FluxDensity*(double)PosY);
		      for (int i=0; i<1; ++i)
			{

			  int InitialIndex = this->EncodeSublatticeIndex(i, PosY, K1,TranslationPhase); // TranlationPhase always one, so can be discarded
			  int FinalIndex = this->EncodeSublatticeIndex(i+1, PosY, K1,TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -this->TunnelElementX[PosY]*PhaseX*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<< -1.0*PhaseX*TranslationPhase<<endl;
#endif
			  FinalIndex = this->EncodeSublatticeIndex(i-1, PosY, K1, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -this->TunnelElementX[PosY]*Conj(PhaseX)*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<< -1.0*Conj(PhaseX)*TranslationPhase<<endl;
#endif		  
			  FinalIndex = this->EncodeSublatticeIndex(i, PosY+1, K1, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -this->TunnelElementY[PosY]*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif

			  FinalIndex = this ->EncodeSublatticeIndex(i, PosY-1, K1, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -this->TunnelElementY[PosY-1]*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif
			}
		    }


      PosY =  NbrSiteY - 1;
      Complex PhaseX = Phase(2.0*M_PI*this->FluxDensity*(double)PosY);
      for (int i=0; i<1; ++i)
      {
	  int InitialIndex = this->EncodeSublatticeIndex(i, PosY, K1,TranslationPhase); // TranlationPhase always one, so can be discarded
	  int FinalIndex = this->EncodeSublatticeIndex(i+1, PosY, K1,TranslationPhase);
	  if (InitialIndex>=FinalIndex)
	    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -this->TunnelElementX[PosY]*PhaseX*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<< -1.0*PhaseX*TranslationPhase<<endl;
#endif
			  FinalIndex = this->EncodeSublatticeIndex(i-1, PosY, K1, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -this->TunnelElementX[PosY]*Conj(PhaseX)*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<< -1.0*Conj(PhaseX)*TranslationPhase<<endl;
#endif		  
			  FinalIndex = this->EncodeSublatticeIndex(i, PosY-1, K1, TranslationPhase);
			  if (InitialIndex>=FinalIndex)
			    TmpOneBodyHamiltonian.AddToMatrixElement(InitialIndex, FinalIndex, -this->TunnelElementY[PosY-1]*TranslationPhase);
#ifdef DEBUG_OUTPUT
			  cout << "H["<<InitialIndex<<"->"<<FinalIndex<<"]="<<-TranslationPhase<<endl;
#endif
  }


#ifdef DEBUG_OUTPUT
	      cout << TmpOneBodyHamiltonian<< endl;
#endif
	    

	      if (this->OneBodyBasis != 0)
		{
		  ComplexMatrix TmpMatrix(this->NbrBands, this->NbrBands, true);
		  TmpMatrix.SetToIdentity();
		  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
		  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
		  this->OneBodyBasis[Index] = TmpMatrix;
		  for (int i = 0; i < this->NbrBands; ++i)
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		}
	      else
		{
		  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
		  TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
		  for (int i = 0; i < this->NbrBands; ++i)
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		}
	    }
    }
}

 
// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// KX = current momentum in x-direction
// translationPhase = phase factor associated with any crossings of unit cell boundary
//
int TightBindingModelHofstadterFiniteCylinder::EncodeSublatticeIndex(int posx, int posy, double KX, Complex &translationPhase)
{
  //cout << "Encoding " << posx<<", "<<posy<<": ";
  int numXTranslations=0;  
  while (posx<0)
    {
      posx+=this->NbrSiteX;
      ++numXTranslations;      
    }
  while (posx>=this->NbrSiteX)
    {
      posx-=this->NbrSiteX;
      --numXTranslations;
    }

  int rst = posx + this->NbrSiteX*posy;
  // determine phase for shifting site to the simulation cell:
  Complex tmpPhase(1.0,0.0);
  Complex tmpPhase2;
  translationPhase=tmpPhase;
  if (numXTranslations>0)
    tmpPhase2=LxTranslationPhase;
  else
    tmpPhase2=Conj(LxTranslationPhase);
  for (int i=0; i<abs(numXTranslations); ++i)
    tmpPhase*=tmpPhase2;
  //cout<<" tmpPhaseX="<<tmpPhase;
  for (int y=1;y<=posy; ++y)
    translationPhase*=tmpPhase;
  translationPhase*=Polar(1.0, KX*numXTranslations);

  //cout << "tX="<<numXTranslations<< ", tY="<<numYTranslations<<", translationPhase= " <<translationPhase<<endl;
  return rst;
}

// get the tight binding hamiltonian in real space 
// 
// return value = tight binding hamiltonian

HermitianMatrix TightBindingModelHofstadterFiniteCylinder::GetRealSpaceTightBindingHamiltonian()
{
  int* NbrConnectedOrbitals = new int [this->NbrBands];
  int** OrbitalIndices = new int* [this->NbrBands];
  int** SpatialIndices = new int* [this->NbrBands];
  Complex** HoppingAmplitudes = new Complex* [this->NbrBands];

  NbrConnectedOrbitals[0] = 3;
if (this->TorusGeometry == true)
   NbrConnectedOrbitals[0] = 4;

 for(int posY = 1 ; posY < this->NbrSiteY - 1; posY++)
      {
 	NbrConnectedOrbitals[posY] = 4;
      }
  NbrConnectedOrbitals[this->NbrSiteY - 1] = 3;
  if (this->TorusGeometry == true)
    NbrConnectedOrbitals[this->NbrSiteY - 1] = 4;


  for (int i = 0; i < this->NbrBands; ++i)
    {
      OrbitalIndices[i] = new int[NbrConnectedOrbitals[i]];
      SpatialIndices[i] = new int[NbrConnectedOrbitals[i]];
      HoppingAmplitudes[i] = new Complex[NbrConnectedOrbitals[i]];
    }

 int   NumXTranslations;
 Complex translationPhase=1.0;
 Complex tmpPhase, tmpPhase2;

  int PosY = 0;
  int TmpPosition =  this->EncodeSublatticeIndex(0, PosY,NumXTranslations,translationPhase);
  int TmpIndex = 0;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(1, PosY,NumXTranslations,translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -this->TunnelElementX[PosY];
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY+1,NumXTranslations,translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] = -this->TunnelElementY[PosY];
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(-1, PosY,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -this->TunnelElementX[PosY];
  ++TmpIndex;

if (this->TorusGeometry == true)
{
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY-1,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -this->TunnelElementY[this->NbrSiteY - 1];
}

for (PosY=1; PosY<this->NbrSiteY - 1; ++PosY)
{
  Complex PhaseX = Phase(2.0*M_PI*this->FluxDensity*(double)PosY);
  TmpPosition =  this->EncodeSublatticeIndex(0, PosY,NumXTranslations,translationPhase);
  int TmpIndex = 0;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(1, PosY,NumXTranslations , translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -this->TunnelElementX[PosY]*PhaseX;
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY+1,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] = -this->TunnelElementY[PosY];
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(-1, PosY,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -this->TunnelElementX[PosY]*Conj(PhaseX);
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY-1,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -this->TunnelElementY[PosY-1];
}


  PosY = this->NbrSiteY - 1;
  Complex PhaseX = Phase(2.0*M_PI*this->FluxDensity*(double)PosY);
   TmpPosition =  this->EncodeSublatticeIndex(0, PosY,NumXTranslations,translationPhase);
   TmpIndex = 0;
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(1, PosY,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -this->TunnelElementX[PosY]*PhaseX;
   ++TmpIndex;
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(-1, PosY,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -this->TunnelElementX[PosY]*Conj(PhaseX);
   ++TmpIndex;
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY-1,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -this->TunnelElementY[PosY-1];
   ++TmpIndex;

if (this->TorusGeometry == true)
{
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY+1,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -this->TunnelElementY[PosY];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  HermitianMatrix TmpMatrix = this->BuildTightBindingHamiltonianRealSpace(NbrConnectedOrbitals, OrbitalIndices, SpatialIndices, HoppingAmplitudes);
  for (int i = 0; i < this->NbrBands; ++i)
    {
      delete[] HoppingAmplitudes[i];
      delete[] SpatialIndices[i];
      delete[] OrbitalIndices[i];
    }
  delete[] HoppingAmplitudes;
  delete[] SpatialIndices;
  delete[] OrbitalIndices;
  delete[] NbrConnectedOrbitals;
  return TmpMatrix;
}


/*
// get the tight binding hamiltonian in real space 
// 
// return value = tight binding hamiltonian

HermitianMatrix TightBindingModelHofstadterFiniteCylinder::GetRealSpaceTightBindingHamiltonian()
{
  int* NbrConnectedOrbitals = new int [this->NbrBands];
  int** OrbitalIndices = new int* [this->NbrBands];
  int** SpatialIndices = new int* [this->NbrBands];
  Complex** HoppingAmplitudes = new Complex* [this->NbrBands];

  NbrConnectedOrbitals[0] = 3;
if (this->TorusGeometry == true)
   NbrConnectedOrbitals[0] = 4;

 for(int posY = 1 ; posY < this->NbrSiteY - 1; posY++)
      {
 	NbrConnectedOrbitals[posY] = 4;
      }
  NbrConnectedOrbitals[this->NbrSiteY - 1] = 3;
  if (this->TorusGeometry == true)
    NbrConnectedOrbitals[this->NbrSiteY - 1] = 4;


  for (int i = 0; i < this->NbrBands; ++i)
    {
      OrbitalIndices[i] = new int[NbrConnectedOrbitals[i]];
      SpatialIndices[i] = new int[NbrConnectedOrbitals[i]];
      HoppingAmplitudes[i] = new Complex[NbrConnectedOrbitals[i]];
    }

 int   NumXTranslations;
 Complex translationPhase=1.0;
 Complex tmpPhase, tmpPhase2;

  int PosY = 0;
  int TmpPosition =  this->EncodeSublatticeIndex(0, PosY,NumXTranslations,translationPhase);
  int TmpIndex = 0;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(1, PosY,NumXTranslations,translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY+1,NumXTranslations,translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] = -1.0;
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(-1, PosY,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
  ++TmpIndex;

if (this->TorusGeometry == true)
{
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY-1,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
}

for (PosY=1; PosY<this->NbrSiteY - 1; ++PosY)
{
  Complex PhaseX = Phase(2.0*M_PI*this->FluxDensity*(double)PosY);
  TmpPosition =  this->EncodeSublatticeIndex(0, PosY,NumXTranslations,translationPhase);
  int TmpIndex = 0;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(1, PosY,NumXTranslations , translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*PhaseX;
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY+1,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] = -1.0;
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(-1, PosY,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*Conj(PhaseX);
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY-1,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
}


  PosY = this->NbrSiteY - 1;
  Complex PhaseX = Phase(2.0*M_PI*this->FluxDensity*(double)PosY);
   TmpPosition =  this->EncodeSublatticeIndex(0, PosY,NumXTranslations,translationPhase);
   TmpIndex = 0;
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(1, PosY,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*PhaseX;
   ++TmpIndex;
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(-1, PosY,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*Conj(PhaseX);
   ++TmpIndex;
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY-1,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
   ++TmpIndex;

if (this->TorusGeometry == true)
{
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY+1,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  HermitianMatrix TmpMatrix = this->BuildTightBindingHamiltonianRealSpace(NbrConnectedOrbitals, OrbitalIndices, SpatialIndices, HoppingAmplitudes);
  for (int i = 0; i < this->NbrBands; ++i)
    {
      delete[] HoppingAmplitudes[i];
      delete[] SpatialIndices[i];
      delete[] OrbitalIndices[i];
    }
  delete[] HoppingAmplitudes;
  delete[] SpatialIndices;
  delete[] OrbitalIndices;
  delete[] NbrConnectedOrbitals;
  return TmpMatrix;
}

*/
  
// build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
//
// nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
// orbitalIndices = array that gives the orbital indices of the connected orbitals
// spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
// hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
// return value = tight binding hamiltonian in real space 

HermitianMatrix  TightBindingModelHofstadterFiniteCylinder::BuildTightBindingHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes)
{
  HermitianMatrix TmpHamiltonian(this->NbrSiteX * this->NbrSiteY, true);
  int   NumXTranslations;
  Complex TmpXPhase = Complex(1.0,0);
  Complex tmp;  
  for (int i = 0;i < this->NbrSiteX ; ++i)
    {
  for (int k = 0;k < this->NbrSiteY ; ++k)
    {
      int Index2 = this->GetRealSpaceTightBindingLinearizedIndex(i,k);
      for (int l = 0; l < nbrConnectedOrbitals[k]; ++l)
	{
	  int Index1 = this->GetRealSpaceTightBindingLinearizedIndexSafe(spatialIndices[k][l] + i, orbitalIndices[k][l],NumXTranslations);                 
          int TmpX = (spatialIndices[k][l] + i)%this->NbrSiteX;
          int TmpY = orbitalIndices[k][l];

          if (Index1 >= Index2)
		  {
	           if (spatialIndices[k][l] == 0)
 	           {
                      TmpHamiltonian.AddToMatrixElement(Index1, Index2, hoppingAmplitudes[k][l]);
                      tmp = hoppingAmplitudes[k][l];
		    }
		    else
		    {
                      TmpHamiltonian.AddToMatrixElement(Index1, Index2, hoppingAmplitudes[k][l]* Phase(2.0*M_PI* this->FluxInserted/((double) this->NbrSiteX)*spatialIndices[k][l]));
                      tmp = hoppingAmplitudes[k][l]* Phase(2.0*M_PI* this->FluxInserted/((double) this->NbrSiteX)*spatialIndices[k][l]);
                    }
                  }
   
#ifdef DEBUG_OUTPUT
         cout << "orbitalIndices[k][l] = "<<orbitalIndices[k][l]  <<" spatialIndices[k][l] = "<<spatialIndices[k][l]<< " hoppingAmplitudes[k][l] = "<<hoppingAmplitudes[k][l]<<endl;
	 cout <<"x = " <<i<< " y = " <<k << " going to X = " <<  TmpX  << " Y = "<<TmpY<<" index1 = "<<Index1 << "; Index2 =" << Index2 <<" Coefficient = " <<  tmp <<"NumTranslation X = "<<NumXTranslations  <<endl;
#endif
			}
		}

	    }
#ifdef DEBUG_OUTPUT
 cout <<TmpHamiltonian<<endl;
#endif
  return TmpHamiltonian;
}


/*
// build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
//
// nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
// orbitalIndices = array that gives the orbital indices of the connected orbitals
// spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
// hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
// return value = tight binding hamiltonian in real space 

HermitianMatrix  TightBindingModelHofstadterFiniteCylinder::BuildTightBindingHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes)
{
  cout <<"in the new function"<<endl;
  cout <<" NbrSites = "<<this->NbrSiteX * this->NbrSiteY<<endl;
  HermitianMatrix TmpHamiltonian(this->NbrSiteX * this->NbrSiteY, true);
  ComplexMatrix TmpHamiltonian2(this->NbrSiteX * this->NbrSiteY,this->NbrSiteX * this->NbrSiteY, true);
  int   NumXTranslations;
  Complex  tmpPhase, tmpPhase2;
  Complex TmpXPhase = Complex(1.0,0);
//  Complex PhaseY = Phase(2.0*M_PI*this->FluxDensity*(double)PosX);

//  Phase(2*M_PI* this->FluxInserted/ this->NbrSiteX);
  for (int i = 0; i < this->NbrSiteX; ++i)

  cout <<" this->LxTranslationPhase  = " << this->LxTranslationPhase<<endl;
  int Index1;

  for (int i = 0; i < this->NbrSiteX; ++i)
    {
  for (int k = 0;k < this->NbrSiteY ; ++k)
    {
      int Index2 = this->GetRealSpaceTightBindingLinearizedIndex(i,k);
      Index1 = this->GetRealSpaceTightBindingLinearizedIndexSafe(i+1,k,NumXTranslations);       

      if (Index1 >= Index2)
          TmpHamiltonian.AddToMatrixElement(Index1, Index2, -1 * Phase(2*M_PI* this->FluxInserted/((double) this->NbrSiteX)));
      Index1 = this->GetRealSpaceTightBindingLinearizedIndexSafe(i-1,k,NumXTranslations);
             
      if (Index1 >= Index2)
     {
        int TmpX = (i - 1)%this->NbrSiteX;
        if ( TmpX<0 )
         	TmpX += this->NbrSiteX;
         TmpHamiltonian.AddToMatrixElement(Index1, Index2, -1 * Phase(-2*M_PI* this->FluxInserted/((double) this->NbrSiteX)));
    }
    if ( k >0  )
    {
      Index1 = this->GetRealSpaceTightBindingLinearizedIndexSafe(i,k-1,NumXTranslations);       
      if (Index1 >= Index2)
          TmpHamiltonian.AddToMatrixElement(Index1, Index2, -1 * Phase(2.0*M_PI*this->FluxDensity*((double) i)));
   }
    if ( k < this->NbrSiteY - 1)
    {
      Index1 = this->GetRealSpaceTightBindingLinearizedIndexSafe(i,k+1,NumXTranslations);       
      if (Index1 >= Index2)
          TmpHamiltonian.AddToMatrixElement(Index1, Index2, -1 * Phase(-2.0*M_PI*this->FluxDensity*((double) i)));
    }
			}
		}

  //cout <<TmpHamiltonian<<endl;
  return TmpHamiltonian;
}

*/

/*
// get the tight binding hamiltonian in real space 
// 
// return value = tight binding hamiltonian

HermitianMatrix TightBindingModelHofstadterFiniteCylinder::GetRealSpaceTightBindingHamiltonian()
{
  int* NbrConnectedOrbitals = new int [this->NbrBands];
  int** OrbitalIndices = new int* [this->NbrBands];
  int** SpatialIndices = new int* [this->NbrBands];
  Complex** HoppingAmplitudes = new Complex* [this->NbrBands];

  NbrConnectedOrbitals[0] = 3;
if (this->TorusGeometry == true)
   NbrConnectedOrbitals[0] = 4;

 for(int posY = 1 ; posY < this->NbrSiteY - 1; posY++)
      {
 	NbrConnectedOrbitals[posY] = 4;
      }
  NbrConnectedOrbitals[this->NbrSiteY - 1] = 3;
  if (this->TorusGeometry == true)
    NbrConnectedOrbitals[this->NbrSiteY - 1] = 4;


  for (int i = 0; i < this->NbrBands; ++i)
    {
      OrbitalIndices[i] = new int[NbrConnectedOrbitals[i]];
      SpatialIndices[i] = new int[NbrConnectedOrbitals[i]];
      HoppingAmplitudes[i] = new Complex[NbrConnectedOrbitals[i]];
    }

 int   NumXTranslations;
 Complex translationPhase=1.0;
 Complex tmpPhase, tmpPhase2;

  int PosY = 0;
  int TmpPosition =  this->EncodeSublatticeIndex(0, PosY,NumXTranslations,translationPhase);
  int TmpIndex = 0;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(1, PosY,NumXTranslations,translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY+1,NumXTranslations,translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] = -1.0;
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(-1, PosY,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
  ++TmpIndex;

if (this->TorusGeometry == true)
{
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY-1,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
}

for (PosY=1; PosY<this->NbrSiteY - 1; ++PosY)
{
  Complex PhaseX = Phase(2.0*M_PI*this->FluxDensity*(double)PosY);
  TmpPosition =  this->EncodeSublatticeIndex(0, PosY,NumXTranslations,translationPhase);
  int TmpIndex = 0;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(1, PosY,NumXTranslations , translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*PhaseX;
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY+1,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] = -1.0;
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(-1, PosY,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*Conj(PhaseX);
  ++TmpIndex;
  OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY-1,NumXTranslations, translationPhase);
  SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
  HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
}


  PosY = this->NbrSiteY - 1;
  Complex PhaseX = Phase(2.0*M_PI*this->FluxDensity*(double)PosY);
   TmpPosition =  this->EncodeSublatticeIndex(0, PosY,NumXTranslations,translationPhase);
   TmpIndex = 0;
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(1, PosY,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*PhaseX;
   ++TmpIndex;
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(-1, PosY,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0*Conj(PhaseX);
   ++TmpIndex;
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY-1,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
   ++TmpIndex;

if (this->TorusGeometry == true)
{
   OrbitalIndices[TmpPosition][TmpIndex] =  this->EncodeSublatticeIndex(0, PosY+1,NumXTranslations, translationPhase);
   SpatialIndices[TmpPosition][TmpIndex] = NumXTranslations;
   HoppingAmplitudes[TmpPosition][TmpIndex] =  -1.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  HermitianMatrix TmpMatrix = this->BuildTightBindingHamiltonianRealSpace(NbrConnectedOrbitals, OrbitalIndices, SpatialIndices, HoppingAmplitudes);
  for (int i = 0; i < this->NbrBands; ++i)
    {
      delete[] HoppingAmplitudes[i];
      delete[] SpatialIndices[i];
      delete[] OrbitalIndices[i];
    }
  delete[] HoppingAmplitudes;
  delete[] SpatialIndices;
  delete[] OrbitalIndices;
  delete[] NbrConnectedOrbitals;
  return TmpMatrix;
}*/
