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



#ifndef TIGHTBINDINGMODELHOFSTADTERSQUARE_H
#define TIGHTBINDINGMODELHOFSTADTERSQUARE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelHofstadterSquare : public Abstract2DTightBindingModel
{

 protected:


  // axis of Landau gauge:
  char LandauGaugeAxis;
  // number of sites in cell in x-direction
  
  // number of sites in cell in y-direction
  int UnitCellX;
  int UnitCellY;

  double tTwo;
  double tThree;
  double Alpha;

  // number of flux quanta in cell (cancelled by opposite flux)
  int NbrFluxQuanta;

  //double MixingAngle;
  //double MixingPhase;

  // auxiliary variables:
  // flux density:
  double FluxDensity;
  // magnetic translation phases;
  Complex LxTranslationPhase;
  Complex LyTranslationPhase;

 public:
  int GetXSitesInUC();
  int GetYSitesInUC();
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
  TightBindingModelHofstadterSquare(int nbrCellX, int nbrCellY, int unitCellX, int unitCellY, int nbrFlux, char axis,
				       double gammaX, double gammaY,
				       AbstractArchitecture* architecture,
                bool storeOneBodyMatrices = true,
                 bool useEmbedding = false,
                  double ttwo=0.0,
                   double tthree=0.0,
                   double alpha=1.0);

  // destructor
  //
  ~TightBindingModelHofstadterSquare();
  ComplexMatrix GetRealSpaceTightBindingEigenstates();

 protected :
  
  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

  // initialize number of flux quanta
  // nbrFluxQuanta = number of quanta of flux piercing the unit cell
  void SetNbrFluxQuanta(int nbrFluxQuanta);

  // code set of quantum numbers posx, posy into a single integer
  // posx = position along x-direction
  // posy = position along y-direction
  // KX = current momentum in x-direction
  // KY = current momentum in y-direction
  // translationPhase = phase factor associated with any crossings of unit cell boundary
  //
  int EncodeSublatticeIndex(int posx, int posy, double KX, double KY, Complex &translationPhase);
  void GetOneBodyBasisFromPython();


};


#endif
