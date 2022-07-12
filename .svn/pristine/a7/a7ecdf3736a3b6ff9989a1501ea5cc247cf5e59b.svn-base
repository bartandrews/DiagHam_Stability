////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on lattice			              //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 12/09/2014                      //
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
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong.h"
#include "GeneralTools/ArrayTools.h"

#include <math.h>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong::BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// nbrSite = number of sites
// xMomentum = momentum sector in the x direction
// maxYMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction 
// memory = amount of memory granted for precalculations

BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong::BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong (int nbrBosons,  int lx, int ly, int xMomentum, int  maxXMomentum, int yMomentum, int maxYMomentum, unsigned long memory):
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong (nbrBosons,lx*ly,xMomentum,maxXMomentum,yMomentum,maxYMomentum,memory), Lx(lx),Ly(ly)
{  
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong::BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong(const BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong & bosons) :  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong(bosons)
{
  this->Lx = bosons.Lx;
  this->Ly = bosons.Ly;
}

// destructor
//

BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong::~BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong ()
{
}


// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong & BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong::operator = (const BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong& bosons)
{
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::operator = (bosons);
  this->Lx = bosons.Lx;
  this->Ly = bosons.Ly;
  return (*this);
}


// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong::Clone()
{
  return new BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong(*this);
}


/*
void BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::GetCompositeFermionWavefunction(ComplexVector & trialState, ComplexMatrix & jastrowEigenVecs,ComplexMatrix & cFEigenVecs, double phaseTranslationX)
{
  int Nx = this->Lx / this->MaxXMomentum;
  int Ny = this->Ly / this->MaxYMomentum;
  int PositionX;
  int PositionY;
  Complex ** ExponentialFactors = new Complex*[this->MaxXMomentum];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      ExponentialFactors[i] = new Complex[this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{ 
	  ExponentialFactors[i][j] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->MaxXMomentum))
							       + (this->YMomentum * ((double) j) / ((double) this->MaxYMomentum))));
	}
    }

#ifdef __LAPACK__
  ComplexLapackDeterminant SlaterCF(NbrBosons);
  ComplexLapackDeterminant SlaterJastrow(NbrBosons);
#else
  ComplexMatrix SlaterCF(NbrBosons, NbrBosons);
  ComplexMatrix SlaterJastrow(NbrBosons, NbrBosons);
#endif

 int NbrTranslation;
 unsigned long * TemporaryState = new unsigned long [this->NbrBosons];
 
 for(int i = 0; i < this->HilbertSpaceDimension ; i++)
   {
     unsigned long * TranslationOfRepresentant = new unsigned long [this->NbrStateInOrbit[i]];
     long * MomentumTable = new long [this->NbrStateInOrbit[i]];
     unsigned long NbrRepresentant = 0;
     unsigned long TmpStateDescription = this->StateDescription[i];
     
     for(int p = 0 ; p < this->MaxYMomentum ;p++)
       {
	 for(int q = 0 ; q < this->MaxXMomentum ;q++)
	   {
	     NbrRepresentant += SearchInArrayAndDefinedWeight(TmpStateDescription,TranslationOfRepresentant,MomentumTable, NbrRepresentant, q* this->MaxYMomentum+p);
	     this->ApplySingleXTranslation(TmpStateDescription);
	   }
	 this->ApplySingleYTranslation(TmpStateDescription);
       }
     
     
     if(NbrRepresentant != this->NbrStateInOrbit[i])
       {
	 cout <<"Wrong Number of Translation for state !"  << i <<endl;
	 cout <<NbrRepresentant<< " " << this->NbrStateInOrbit[i]<<endl;
       }
     
     for(int k =0 ;k <NbrRepresentant; k++)
       {
	 //   cout <<TranslationOfRepresentant[k]<<" "<<MomentumTable[k]<<endl;
	 this->ConvertToMonomial(TranslationOfRepresentant[k],TemporaryState);
	 this->GetPositionSum(TemporaryState, PositionX, PositionY);
	 for (int p = 0; p < NbrBosons; ++p)
	   {
	     for (int q = 0; q < NbrBosons; ++q)
	       {
		 SlaterCF.SetMatrixElement(p,q,cFEigenVecs[p][(int) TemporaryState[q]]);
		 SlaterJastrow.SetMatrixElement(p,q,jastrowEigenVecs[p][(int) TemporaryState[q]]);
	       }	      
	   }
	 //    trialState[i] +=  Conj(ExponentialFactors[MomentumTable[k]/this->MaxYMomentum][MomentumTable[k]%this->MaxYMomentum]) * SlaterCF.Determinant() * SlaterJastrow.Determinant() * Phase(-1.0*PositionY*(MomentumTable[k]/this->MaxYMomentum)*phaseTranslationX);
	 trialState[i] +=  ExponentialFactors[MomentumTable[k]/this->MaxYMomentum][MomentumTable[k]%this->MaxYMomentum] * SlaterCF.Determinant() * SlaterJastrow.Determinant() * Phase(1.0*PositionY*(MomentumTable[k]/this->MaxYMomentum)*phaseTranslationX);
	 
       }
     trialState[i] /= sqrt(this->NbrStateInOrbit[i]);
     delete [] TranslationOfRepresentant;
   }
 for (int i = 0; i < this->MaxXMomentum; ++i)
   { 
     delete [] ExponentialFactors[i];
   }
 delete []  ExponentialFactors;
}



*/
