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


#ifndef TIGHTBINDINGMODELOFLNORBITALTRIANGULARLATTICE_H
#define TIGHTBINDINGMODELOFLNORBITALTRIANGULARLATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelOFLNOrbitalTriangularLattice : public Abstract2DTightBindingModel
{

 protected:

  int NbrStep;
  int NbrInternalDegree;
  double LaserStrength;
  double InvMomentum;

  int ChernNumber;


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
  TightBindingModelOFLNOrbitalTriangularLattice(double laserStrength, int nbrInternalDegree, int nbrSiteX,  int nbrSiteY, int chernNumber, double gammaX, double gammaY, AbstractArchitecture* architecture, int cutOFF, bool storeOneBodyMatrices = true);

  // constructor from a binary file
  //
  // fileName = name of the binary file that contains the band structure information
  TightBindingModelOFLNOrbitalTriangularLattice(char* fileName);

  // destructor
  //
  ~TightBindingModelOFLNOrbitalTriangularLattice();
  
  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);
 
  // get linearized indices
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  int GetIntermediateLinearizedIndices(int xMomentum, int yMomentum,int spin); 

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
//  virtual void CoreComputeBandStructure(long minStateIndex = 0l, long nbrStates= 0l);

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex = 0l, long nbrStates= 0l);
  
  // write an header that describes the tight binding model
  // 
  // output = reference on the output stream
  // return value  = reference on the output stream
  virtual ofstream& WriteHeader(ofstream& output);
  
};


inline int TightBindingModelOFLNOrbitalTriangularLattice::GetIntermediateLinearizedIndices(int xMomentum, int yMomentum,int spin)
{
  int TmpXMomentum = xMomentum;
  int TmpYMomentum = yMomentum;
  int TmpSpin = spin;
 
 if(spin <0)
    {
      TmpSpin+=this->NbrInternalDegree;
      TmpXMomentum--;
    }
 if(spin >= this->NbrInternalDegree)
   {
     TmpSpin-=this->NbrInternalDegree;
     TmpXMomentum++;
   }
 
 if (TmpXMomentum < 0)
   TmpXMomentum += this->NbrStep;
 if (TmpYMomentum < 0)
   TmpYMomentum += this->NbrStep;
 
 return this->NbrInternalDegree*((TmpXMomentum%this->NbrStep)*this->NbrStep+(TmpYMomentum%this->NbrStep))+TmpSpin;
}


#endif
