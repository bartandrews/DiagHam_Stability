////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2017 Nicolas Regnault                  //
//                                                                            //
//                         class Author: Gunnar Möller                        //
//                                                                            //
//              class providing interface for interactions on lattice         //
//                                                                            //
//                        last modification : 16/01/2017                      //
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


#ifndef TIGHTBINDINGINTERACTIONCOULOMB2DEWALD_H
#define TIGHTBINDINGINTERACTIONCOULOMB2DEWALD_H


#include "config.h"

/* #include <iostream> */
/* using std::ostream; */
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"
#include "AbstractTightBindingModel.h"
#include "Abstract2DTightBindingModel.h"
#include "AbstractTightBindingInteraction.h"

class TightBindingInteractionCoulomb2DEwald : public AbstractTightBindingInteraction
{

 protected:
  Abstract2DTightBindingModel *TightBindingModel; // pointer to a tight-binding model
  double Alpha;           // parameter for Ewald sum
  double HalfInvAlpha;    // inverse parameter for Ewald sum
  double Ratio;           // ratio of length of e2 to length of e1

  double SimulationCellArea;   // size of the simulation cell
  double SelfInteraction;      // calculate self-interaction for a site with itself
  double BackgroundCorrection; // background correction due to non-neutral cell

  int NbrCells1;          // no. of sites in sim. cell in 1 direction
  int NbrCells2;          // no. of sites in sim. cell in 2 direction
  int NbrImages;             // number of images used for summation
    
  RealVector BasisVector1;         // first basis vector
  RealVector BasisVector2;         // second basis vector
    
  RealVector K1;          // first reciprocal lattice vector
  RealVector K2;          // second reciprocal lattice vector
    
  RealMatrix U;           // matrix storing the potential terms, format: sublattice1; linearizedIndex2

 public:
  
  TightBindingInteractionCoulomb2DEwald(Abstract2DTightBindingModel *tightBinding, double alpha=1.0, int nbrImages=100);
  
  // default destructor
  virtual ~TightBindingInteractionCoulomb2DEwald(){}

  // accessor routine to get the magnitude of the interaction
  // obtain Amplitude for a given separation
  // s = sublattice of the initial site
  // dR1 = displacement final site in unit cells along 1-direction
  // dR2 = displacement final site in unit cells along 2-direction
  // s2 = sublattice index of the final site
  double GetAmplitude(int s, int dR1, int dR2, int s2);
               
  // calculate lattice vectors
  void GetLatticeVector(double i, double j, RealVector& Rst);
    
  // calculate reciprocal lattice vectors
  void GetReciprocalLatticeVector(double M, double theta, RealVector& Rst);
  
  // calculate real space term in Ewald summation
  double CoulombErfc(double m, double n, RealVector &r);
  
  // calculate Fourier term in Ewald summation
  double FourierTerm(double u, double v, RealVector &r);
      
    
  // carries out Ewald summation
  void PrecalculateInteractions();
        
  // obtain Amplitude for a given separation
  double GetAmplitude(int dx, int dy);


private:
  int GetSiteIndex(int dx, int dy) {return dx+this->NbrCells1*dy;} // return the site index, given a displacement in the unit cell


  RealVector TmpVector; // a 2D vector used for temporary storage in subroutines
};

#endif // TIGHTBINDINGINTERACTIONCOULOMB2DEWALD_H
