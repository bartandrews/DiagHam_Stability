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

#ifndef TIGHTBINDINGINTERACTIONCOULOMBFORHOFSTADTERLATTICE_H
#define TIGHTBINDINGINTERACTIONCOULOMBFORHOFSTADTERLATTICE_H

#include "Vector/RealVector.h"
#include "Abstract2DTightBindingModel.h"
#include "AbstractTightBindingInteraction.h"

class TightBindingInteractionCoulombForHofstadterLattice : public AbstractTightBindingInteraction
{
protected:

  Abstract2DTightBindingModel *TightBindingModel;    // pointer to a tight-binding model
  double twist_angle; 	  // angle between unit vectors (radians)
  int N_sites1;           // no. of sites in sim. cell in 1 direction
  int N_sites2;           // no. of sites in sim. cell in 2 direction
  double alpha;           // Ewald splitting parameter
  double HalfInvAlpha;    // 1/(2*alpha)
  double r;               // ratio of length of e2 to length of e1
  int images;             // number of images used for summation
    
  RealVector bv1;         // first basis vector
  RealVector bv2;         // second basis vector
    
  RealVector k1;          // first reciprocal lattice vector
  RealVector k2;          // second reciprocal lattice vector  
    
  RealVector U;           // vector storing the potential terms

    
public:
  // constructor
  TightBindingInteractionCoulombForHofstadterLattice(Abstract2DTightBindingModel *tightBinding, double alpha=1, int images=100);
    
  // total number of sites
  int GetTotalSites();
    
  // area of unit cell
  double GetUnitCellArea();
    
  // index each site
  void ConvertSiteIndex(int i, RealVector &Coords);
   
  // calculate vector from polar angle
  void GetPolarVector(double M, double theta, RealVector& Rst);
    
  // calculate lattice vectors
  void GetLatticeVector(double i, double j, RealVector& Rst);
    
  // calculate reciprocal lattice vectors
  void GetReciprocalLatticeVector(double M, double theta, RealVector& Rst);
  
  // calculate real space term in Ewald summation
  double CoulombErfc(double m, double n, RealVector r);
  
  // calculate Fourier term in Ewald summation
  double FourierTerm(double u, double v, RealVector r);
    
  // calculate self-energy
  double SelfInteraction();
  
  // calculate correction due to non-neutral cell
  double BackgroundCorrection();
    
  // carries out Ewald summation
  void DoSummation();
        
  // obtain Amplitude for a given separation
  double GetAmplitude(int dx, int dy);

  // accessor routine to get the magnitude of the interaction - overriding function provided in interface
  // obtain Amplitude for a given separation
  // s = sublattice of the initial site
  // dR1 = displacement final site in unit cells along 1-direction
  // dR2 = displacement final site in unit cells along 2-direction
  // s2 = sublattice index of the final site
  double GetAmplitude(int s, int dR1, int dR2, int s2);

private:
  int GetSiteIndex(int dx, int dy) {return dx+this->N_sites1*dy;} // return the site index, given a displacement in the unit cell


  RealVector TmpVector; // a 2D vector used for temporary storage in subroutines
};

#endif // TIGHTBINDINGINTERACTIONCOULOMBFORHOFSTADTERLATTICE_H
