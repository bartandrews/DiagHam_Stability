////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class authors: Benjamin Huddart & Gunnar Möller               //
//                                                                            //
// Class to compute Coulomb energies between any two sites within a 2D        //
// lattice using the Ewald summation method                                   //
//                                                                            //
// Potentials for all possible separation stored in vector U                  //
// Values accessed using class member function GetAmplitude                   //
// Energies in units of 1/(4*pi*epsilon_0*a).                                 //
// Lattice constant a=1                                                       //
//                                                                            //
//                        last modification : 05/08/2016                      //
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
#include "TightBindingInteractionCoulomb2DEwald.h"

#include "Vector/RealVector.h"
#include <cmath>
#include <cassert>

#ifdef HAVE_ERFC
using std::erfc;
#else
#define erfc(x)  (1.0-erf(x))
#endif

// Constructor
// tightBinding = pointer to the relevant tightBindingModel
// angle = obliquity angle of lattice (radians)
// N1 = number of lattice sites along e1 = (1,0) direction
// N2 = number of lattice sites along e2 = (sin(angle),cos(angle)) direction
// alpha = Ewald parameter
// ratio = aspect ratio of lengths |e2|/|e1|
// nbrImages = number of image charges to use
TightBindingInteractionCoulomb2DEwald::TightBindingInteractionCoulomb2DEwald(Abstract2DTightBindingModel *tightBinding,  double alpha, int nbrImages)
{
  this->TightBindingModel = tightBinding;
  this->NbrCells1 = tightBinding->GetNbrSiteX();
  this->NbrCells2 = tightBinding->GetNbrSiteY();
  this->Alpha=alpha;
  this->HalfInvAlpha = 0.5/this->Alpha;
  this->NbrImages=nbrImages;
  
  // get vectors spanning simulation cell, and their reciprocal lattice vectors
  tightBinding->GetLatticeVector(this->BasisVector1, 0);
  this->BasisVector1*=this->NbrCells1;
  tightBinding->GetLatticeVector(this->BasisVector2, 1);
  this->BasisVector2*=this->NbrCells2;

  this->SimulationCellArea = std::fabs(this->BasisVector2 * Cross(this->BasisVector1));
  this->SelfInteraction = -2.*this->Alpha/sqrt(M_PI);
  this->BackgroundCorrection = -2.*sqrt(M_PI) / (this->Alpha*this->SimulationCellArea);

  tightBinding->GetReciprocalLatticeVector(this->K1, 0);
  this->K1/=this->NbrCells1;
  tightBinding->GetReciprocalLatticeVector(this->K2, 1);
  this->K2/=this->NbrCells2;

  assert(fabs(this->BasisVector1*this->K1 - 2.*M_PI)<1e-13 && fabs(this->BasisVector1*this->K2)<1e-13
	 && fabs(this->BasisVector2*this->K2 - 2.*M_PI)<1e-13 && fabs(this->BasisVector2*this->K1)<1e-13 );

  // think about ratio
  this->Ratio=1.0; // ratio; /// @bug Ratio may not be well defined

  this->TmpVector.Resize(2);
  
  this->PrecalculateInteractions();
    
};


void TightBindingInteractionCoulomb2DEwald::GetLatticeVector(double i, double j, RealVector& Rst) {
  Rst.ClearVector();
  Rst.AddLinearCombination(i, BasisVector1, j, BasisVector2);
};

double TightBindingInteractionCoulomb2DEwald::CoulombErfc(double m, double n, RealVector &X) 
{
  this->TmpVector.Copy(X);
  TmpVector.AddLinearCombination(m, BasisVector1, n, BasisVector2);
  double dist=TmpVector.Norm();
  if (dist != 0)
    return erfc(this->Alpha*dist)/dist;
  else return 0;
};

double TightBindingInteractionCoulomb2DEwald::FourierTerm(double u, double v, RealVector &R) {
  this->TmpVector.ClearVector();
  TmpVector.AddLinearCombination(u, this->K1, v, this->K2);
  double kappa=TmpVector.Norm();
  double phase=TmpVector*R; //dot product
  if (kappa > 1e-15) {
    return 2.0*M_PI/kappa * (1.0+std::cos(phase)) * erfc(kappa*this->HalfInvAlpha);
  }
  else return 0.0;
};
    
// carry out Ewald summation
void TightBindingInteractionCoulomb2DEwald::PrecalculateInteractions() 
{
  this->U.ResizeAndClean(this->TightBindingModel->GetNbrBands(), this->TightBindingModel->GetNbrSites());

  int index2;
  double invA=1.0/this->SimulationCellArea;
  RealVector zero=RealVector(2,true);
  RealVector unitCellIndices(2);
  RealVector subl1(2);
  RealVector R(2);
        
  for (int s=0; s<this->TightBindingModel->GetNbrBands(); ++s)
    {
      this->TightBindingModel->GetSublatticeVector(subl1, s);
      for (int x=0; x<this->NbrCells1; ++x)
	for (int y=0; y<this->NbrCells2; ++y)
	  {
	    unitCellIndices[0]=x;
	    for (int s2=0; s2<this->TightBindingModel->GetNbrBands(); ++s2)
	      {
		unitCellIndices[1]=y;
		// get position of site 2
		this->TightBindingModel->GetSitePosition(R, unitCellIndices, s2);
		// calculate relative separation
		R-=subl1;
		index2 = this->TightBindingModel->GetRealSpaceTightBindingLinearizedIndex(x, y, s2);
		            
		// compute real space contribution
		double CoulombSum=0.0;
		// sum over periodic repetitions
		for(double m=-NbrImages;m<NbrImages+1; m++){
		  for(double n=-NbrImages; n<NbrImages+1; n++) {
		    CoulombSum += CoulombErfc(m, n, R); // interactions between different charges
		    CoulombSum += CoulombErfc(m, n, zero); // interaction between charge and its own periodic image
		  }
		}
      
		// compute Fourier space contribution
		double FourierSum=0.0;
		// compute Fourier term
		for (double u=-NbrImages; u<NbrImages+1; u++) {
		  for(double v=-NbrImages; v<NbrImages+1; v++) {
		    FourierSum += FourierTerm(u,v,R);
		  }
		}
      
		FourierSum *= invA;
		double V=CoulombSum+FourierSum;
     
		V += this->SelfInteraction;
		V += this->BackgroundCorrection;
		
		U[s][index2]=V;
	      }
	  }
    }
}

// accessor routine to get the magnitude of the interaction
// obtain Amplitude for a given separation
// s = sublattice of the initial site
// dR1 = displacement final site in unit cells along 1-direction
// dR2 = displacement final site in unit cells along 2-direction
// s2 = sublattice index of the final site
double TightBindingInteractionCoulomb2DEwald::GetAmplitude(int s, int dR1, int dR2, int s2)
{
  int idx2 = this->TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(dR1, dR2, s2);
  return U[s][idx2];
}
