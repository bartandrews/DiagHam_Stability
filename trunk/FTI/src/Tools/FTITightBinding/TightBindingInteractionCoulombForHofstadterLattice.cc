////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2017 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class authors: Benjamin Huddart & Gunnar Möller               //
//                                                                            //
// Class to compute Coulomb energies between any two sites within a 2D        //
// lattice using the Ewald summation method, direcly from Ben Huddart's source//
//                                                                            //
// Potentials for all possible separation stored in vector U                  //
// Values accessed using class member function GetAmplitude                   //
// Energies in units of 1/(4*pi*epsilon_0*a).                                 //
// Lattice constant a=1                                                       //
//                                                                            //
//                        last modification : 02/03/2017                      //
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

#include "TightBindingInteractionCoulombForHofstadterLattice.h"

#include "Vector/RealVector.h"
#include <cmath>

using std::sqrt;
#ifdef HAVE_ERFC
using std::erfc;
#else
#define erfc(x)  (1.0-erf(x))
#endif


// Constructor
// tightBinding = Tight Binding model, which supplies some of the parameters named, below
// angle = obliquity angle of lattice (radians)
// N1 = number of lattice sites along e1 = (1,0) direction
// N2 = number of lattice sites along e2 = (sin(angle),cos(angle)) direction
// alpha = Ewald parameter
// ratio = aspect ratio of lengths |e2|/|e1|
// images = number of image charges to use
TightBindingInteractionCoulombForHofstadterLattice::TightBindingInteractionCoulombForHofstadterLattice(Abstract2DTightBindingModel *tightBinding, double alpha, int images)
{
  this->TightBindingModel = tightBinding;
  int numX, numY;
  tightBinding->GetMUCDimensions(numX, numY);
  this->twist_angle=tightBinding->GetTwistAngle();
  this->N_sites1=tightBinding->GetNbrSiteX()*numX;
  this->N_sites2=tightBinding->GetNbrSiteY()*numY;
  this->alpha=alpha;
  this->HalfInvAlpha = 0.5/this->alpha;
  //this->r=ratio;
  this->r=1;  // correct choice for regular Hofstadter lattices.
  this->images=images;
    
  this->GetPolarVector(1,0,this->bv1);
  this->GetPolarVector(this->r,this->twist_angle,this->bv2);
    
  this->GetReciprocalLatticeVector((double)this->N_sites2*r, this->twist_angle, this->k1);
  this->GetReciprocalLatticeVector((double)this->N_sites1, 0, this->k2);
  this->k2 *= -1.0;

  this->TmpVector.Resize(2);
  
  this->DoSummation();    
};

double TightBindingInteractionCoulombForHofstadterLattice::GetUnitCellArea() {
  return r*std::sin(twist_angle);
};

int TightBindingInteractionCoulombForHofstadterLattice::GetTotalSites() {
  return N_sites1*N_sites2;
};

// define site indices for each position
// @param[in] index - integer site index
// @param[out] Coords - coordinates for given site index
void TightBindingInteractionCoulombForHofstadterLattice::ConvertSiteIndex(int index, RealVector &Coords){
  Coords[0]=index % N_sites1;
  Coords[1]=index / N_sites1;
};

void TightBindingInteractionCoulombForHofstadterLattice::GetLatticeVector(double i, double j, RealVector& Rst) {
  Rst.ClearVector();
  Rst.AddLinearCombination(i, bv1, j, bv2);
};

void TightBindingInteractionCoulombForHofstadterLattice::GetReciprocalLatticeVector(double M, double theta, RealVector& Rst)
{
  Rst.Resize(2);
  Rst[0]=2.*M_PI*M*std::sin(theta)/(GetTotalSites()*GetUnitCellArea());
  Rst[1]=-2.*M_PI*M*std::cos(theta)/(GetTotalSites()*GetUnitCellArea());
};

double TightBindingInteractionCoulombForHofstadterLattice::CoulombErfc(double m, double n, RealVector X) {
    this->TmpVector.ClearVector();
    TmpVector += X;
    TmpVector.AddLinearCombination(m*N_sites1, bv1, n*N_sites2, bv2);
  double dist=TmpVector.Norm();
  if (dist != 0) {
    return erfc(alpha*dist)/dist;
  }
  else return 0;
};

double TightBindingInteractionCoulombForHofstadterLattice::FourierTerm(double u, double v, RealVector r) {
  this->TmpVector.ClearVector();
  TmpVector.AddLinearCombination(u, k1, v, k2);
  double kappa=TmpVector.Norm();
  double phase=TmpVector*r; //dot product
  if (kappa > 1e-15) {
    return 2.0*M_PI/kappa * (1.0+std::cos(phase)) * erfc(kappa*this->HalfInvAlpha);
  }
  else return 0;
};

double TightBindingInteractionCoulombForHofstadterLattice::SelfInteraction() {
  return -2.*alpha/sqrt(M_PI);
};

double TightBindingInteractionCoulombForHofstadterLattice::BackgroundCorrection() {
  return -2.*sqrt(M_PI) /(alpha*GetTotalSites()*GetUnitCellArea());
};

// calculate 2D vector from polar form
void TightBindingInteractionCoulombForHofstadterLattice::GetPolarVector(double M, double theta, RealVector& Rst) {
  Rst.Resize(2);
  Rst[0]=M*std::cos(theta);
  Rst[1]=M*std::sin(theta);
};
    
// carry out Ewald summation
void TightBindingInteractionCoulombForHofstadterLattice::DoSummation() {
  int sites=GetTotalSites();
  double invA=1.0/(sites*GetUnitCellArea());
  this->U.Resize(sites);
  RealVector coordinateIndices(2);
  RealVector R(2);
        
  for (int i=0; i<sites; i++) {
    // find relative position for each pair of indices and assign correct energy
    ConvertSiteIndex(i, coordinateIndices);
            
    this->GetLatticeVector(coordinateIndices[0], coordinateIndices[1], R);
    RealVector zero=RealVector(2,true);
            
    // compute real space contribution
    double CoulombSum=0;
    // sum over periodic repetitions
    for(double m=-images;m<images+1; m++){
      for(double n=-images; n<images+1; n++) {
	CoulombSum += CoulombErfc(m, n, R); // interactions between different charges
	CoulombSum += CoulombErfc(m, n, zero); // interaction between charge and its own periodic image
      }
    }
      
    // compute Fourier space contribution
    double FourierSum=0.0;
    // compute Fourier term
    for (double u=-images; u<images+1; u++) {
      for(double v=-images; v<images+1; v++) {
	FourierSum += FourierTerm(u,v,R);
      }
    }
      
    FourierSum *= invA;
    double V=CoulombSum+FourierSum;
     
    V += SelfInteraction();
    V += BackgroundCorrection();
            
    U[i]=V;
  }
}

// obtain Amplitude for a given separation
double TightBindingInteractionCoulombForHofstadterLattice::GetAmplitude(int dx, int dy)
{   
  while (dx < 0) {dx += N_sites1;}
  while (dy < 0) {dy += N_sites2;}
  while( dx > N_sites1 ) { dx -= N_sites1; }
  while( dy > N_sites2 ) { dx -= N_sites2; }
  return U[this->GetSiteIndex(dx,dy)];
};

// accessor routine to get the magnitude of the interaction - overriding function provided in interface
// obtain Amplitude for a given separation
// s = sublattice of the initial site
// dR1 = displacement final site in unit cells along 1-direction
// dR2 = displacement final site in unit cells along 2-direction
// s2 = sublattice index of the final site
double TightBindingInteractionCoulombForHofstadterLattice::GetAmplitude(int s, int dR1, int dR2, int s2)
{

  int numX, numY;
  this->TightBindingModel->GetMUCDimensions(numX, numY);
  int xI, xF, yI, yF;
  this->TightBindingModel->DecodeSublatticeIndex(s, xI, yI);
  this->TightBindingModel->DecodeSublatticeIndex(s2, xF, yF);
  if ((s==s2) && (dR1==0) && (dR2==0))
    return 0.0;
  else
    {
      int dx = xF - dR1 * numX - xI;
      int dy = yF - dR2 * numY - yI;
      return this->GetAmplitude(dx, dy);
    }
}

