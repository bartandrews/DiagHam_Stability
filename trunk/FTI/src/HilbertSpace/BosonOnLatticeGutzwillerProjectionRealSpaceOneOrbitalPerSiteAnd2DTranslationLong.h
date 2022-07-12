////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of bosons on lattice                          //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 11/09/2014                      //
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


#ifndef BOSONONLATTICEGUTZWILLERPROJECTIONONEORBITALPERSITEREALSPACEAND2DTRANSLATIONLONG_H
#define BOSONONLATTICEGUTZWILLERPROJECTIONONEORBITALPERSITEREALSPACEAND2DTRANSLATIONLONG_H

#include "config.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong.h"

#include <iostream>



class BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong : public  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong
{

 protected:
   
   int Lx;
   int Ly;
	
 public:

  // default constructor		
  // 
  BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong ();

  // basic constructor
  // 
  // nbrBosons = number of fermions
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong (int nbrBosons, int lx, int ly, int xMomentum, int maxXMomentum,
										int yMomentum, int maxYMomentum, unsigned long memory = 10000000);
  
  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong (const  BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong& bosons);
  
  // destructor
  //
  ~BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong & operator = (const BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong & bosons);
  
  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // compute sum of positions in the x and y direction for lattice class
  //
  // index = index of the state in the basis whose position sums are to be computed
  // positionX = reference on the sum of positions in the x direction
  // positionY = reference on the sum of positions in the y direction
  virtual void GetPositionSum(int index,int & positionX, int & positionY);
  virtual void GetPositionSum(unsigned long * monomial, int & positionX, int & positionY);

  //  virtual void GetCompositeFermionWavefunction(ComplexVector & trialState, ComplexMatrix & jastrowEigenVecs,ComplexMatrix & cFEigenVecs, double phaseTranslationX);
};

// compute sum of positions in the x and y direction for lattice class
//
// index = index of the state in the basis whose position sums are to be computed
// positionX = reference on the sum of positions in the x direction
// positionY = reference on the sum of positions in the y direction

inline void BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong::GetPositionSum(int index, int & positionX, int & positionY)
{
  //  cout <<"inside inline void BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::GetPositionSum(int index, int & positionX, int & positionY)"<<endl;
  positionX =0;
  positionY =0;
  int Nx = this->Lx / this->MaxXMomentum;
  int Ny = this->Ly / this->MaxYMomentum;
  
//cout <<"Nx = "<<Nx << endl;
//  cout <<"Ny = "<<Ny << endl;
//  cout <<"this->StateXShift = "<<this->StateXShift<<endl;
//  cout <<"this->StateYShift = "<<this->StateYShift<<endl;
//  cout <<"this->StateDescription[index] = "<<this->StateDescription[index]<<endl;
  
  for (int i = 0; i < this->NbrSite; i++)
    { 	
      if ( ((this->StateDescription[index]>>   ((ULONGLONG) i) ) & ((ULONGLONG) 0x1ul) ) ==  ((ULONGLONG) 0x1ul) )
	{
	  positionX +=  (i/this->StateXShift) *  Nx;
	  //          cout <<"positionX = "<<positionX <<endl;
	  positionY += ((i%this->StateXShift)/this->StateYShift) * Ny;
//          cout <<"positionY = "<<positionY <<endl;
	  positionX += (i % Nx); 
//          cout <<"positionX = "<<positionX <<endl;
	  positionY += ((i%this->StateYShift)/ Nx) ;
//          cout <<"positionY = "<<positionY <<endl;
	}
    }
}

// compute sum of positions in the x and y direction for lattice class
//
// index = index of the state in the basis whose position sums are to be computed
// positionX = reference on the sum of positions in the x direction
// positionY = reference on the sum of positions in the y direction

inline void BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong::GetPositionSum(unsigned long * monomial, int & positionX, int & positionY)
{
//  cout <<"inside inline void BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::GetPositionSum(unsigned long * monomial, int & positionX, int & positionY)"<<endl;
  positionX =0;
  positionY =0;
  int Nx = this->Lx / this->MaxXMomentum;
  int Ny = this->Ly / this->MaxYMomentum;
  
  for (int i = 0; i < this->NbrBosons; i++)
    { 	
	 positionX +=  ( monomial[i]/this->StateXShift) *  Nx;
         positionY += ((monomial[i]%this->StateXShift)/this->StateYShift) * Ny;
         positionX += (monomial[i] % Nx);
         positionY += ((monomial[i]%this->StateYShift)/ Nx) ;
    }
}

#endif
