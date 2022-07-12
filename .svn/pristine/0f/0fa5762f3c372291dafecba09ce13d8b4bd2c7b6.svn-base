////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the Checkerboard lattice       //
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


#ifndef TIGHTBINDINGMODELHOFSTADTERFINITECYLINDER_H
#define TIGHTBINDINGMODELHOFSTADTERFINITECYLINDER_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelHofstadterFiniteCylinder : public Abstract2DTightBindingModel
{

 protected:

  // number of sites in cell in x-direction
  int NbrSiteX;
  // number of sites in cell in y-direction
  int NbrSiteY;

  // number of flux quanta in cell (cancelled by opposite flux)
  int NbrFluxQuanta;
  
  // auxiliary variables:
  // flux density:
  double FluxDensity;
  double FluxInserted;
  double * TunnelElementX;
  double * TunnelElementY;
  // magnetic translation phases;
  Complex LxTranslationPhase;
  Complex LyTranslationPhase;
  bool TorusGeometry;

 public:

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
  TightBindingModelHofstadterFiniteCylinder(int nbrSiteX, int nbrSiteY, int nbrFlux, double * tunnelElementX, double * tunnelElementY, char axis,double gammaX, double gammaY,  AbstractArchitecture* architecture,   double fluxInserted = 0, bool storeOneBodyMatrices = true, bool torusGeometry = false);

  // destructor
  //
  ~TightBindingModelHofstadterFiniteCylinder();


  // get the tight binding hamiltonian in real space 
  // 
  // return value = tight binding hamiltonian
  HermitianMatrix GetRealSpaceTightBindingHamiltonian();

  HermitianMatrix  BuildTightBindingHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes);

  // get the index of the real space tight binding model from the real space coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndex(int x, int y) const;
  int  GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y) const;  

  double GetFluxDensity() const{return FluxDensity;}
  int GetNumberSiteX() const{return  NbrSiteX;}
  int GetNumberSiteY() const{return  NbrSiteY;}


 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

  // code set of quantum numbers posx, posy into a single integer
  // posx = position along x-direction
  // posy = position along y-direction
  // KX = current momentum in x-direction
  // translationPhase = phase factor associated with any crossings of unit cell boundary
  //
  int EncodeSublatticeIndex(int posx, int posy, double KX, Complex &translationPhase);


 // code set of quantum numbers posx, posy into a single integer
 // posx = position along x-direction
 // posy = position along y-direction
 // numXTranslations = number of translation in the x direction to get back to the unit cell 
 // numXTranslations = number of translation in the y direction to get back to the unit cell
 //
 int  EncodeSublatticeIndex(int posx, int posy,int & numXTranslations, Complex &translationPhase);
  

 int  GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int & numXTranslations);

};


// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// numXTranslations = number of translation in the x direction to get back to the unit cell 

inline int TightBindingModelHofstadterFiniteCylinder::EncodeSublatticeIndex(int posx, int posy,int & numXTranslations, Complex &translationPhase) 
{
  numXTranslations=0;

  while (posx<0)
    {
      posx=0;
      ++numXTranslations;      
    }
  while (posx>0)
    {
      posx=0;
      --numXTranslations;
    }

  if (posy < 0)
    {
      posy+=NbrSiteY;
    }
  while (posy >= NbrSiteY)
    {
      posy-=NbrSiteY;
    }


  Complex tmpPhase(1.0,0.0);
  Complex tmpPhase2;
  translationPhase=tmpPhase;
  if (numXTranslations>0)
    tmpPhase2=LxTranslationPhase;
  else
    tmpPhase2=Conj(LxTranslationPhase);
  for (int i=0; i<abs(numXTranslations); ++i)
    tmpPhase*=tmpPhase2;
  tmpPhase=1.0;
  return posy;
}
  


// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// return value = linearized index  

inline int  TightBindingModelHofstadterFiniteCylinder::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int & numXTranslations)
{
  numXTranslations=0;

  if(x >= this->NbrSiteX)
  {
    x -=  this->NbrSiteX;
    numXTranslations--;
  }
  if (x < 0)
  {
    x +=  this->NbrSiteX;
    numXTranslations++;
  }

  return this->GetRealSpaceTightBindingLinearizedIndex(x, y); 
}

// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// return value = linearized index  

inline int  TightBindingModelHofstadterFiniteCylinder::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y) const
{
  if(x >= this->NbrSiteX)
  {
    x -=  this->NbrSiteX;
  }
  if (x < 0)
  {
    x +=  this->NbrSiteX;
  }
  return this->GetRealSpaceTightBindingLinearizedIndex(x, y); 
}


// get the index of the real space tight binding model from the real space coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int TightBindingModelHofstadterFiniteCylinder::GetRealSpaceTightBindingLinearizedIndex(int x, int y) const
{
  return (y  + x * this->NbrSiteY); 
}


#endif
