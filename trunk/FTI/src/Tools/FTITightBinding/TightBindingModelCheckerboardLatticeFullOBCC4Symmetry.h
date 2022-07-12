////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the checkerboard lattice       //
//                with full open boundary conditions compatible               //
//                         with the C4 symmetry implentation                  //
//                                                                            //
//                        last modification : 28/02/2018                      //
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


#ifndef TIGHTBINDINGMODELCHECKERBOARDLATTICEFULLOBCC4SYMMETRY_H
#define TIGHTBINDINGMODELCHECKERBOARDLATTICEFULLOBCC4SYMMETRY_H


#include "config.h"
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLatticeFullOBC.h"


class TightBindingModelCheckerboardLatticeFullOBCC4Symmetry : public TightBindingModelCheckerboardLatticeFullOBC
{

 protected:

  // map from the canonical linearized index to the C4 compatible one
  int* LinearizedIndexMap;
  // map from the C4 compatible linearized index to the canonical one
  int* InvertLinearizedIndexMap;

 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x (or y) direction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // mus = sublattice chemical potential on A sites
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelCheckerboardLatticeFullOBCC4Symmetry(int nbrSiteX, double t1, double t2, double mus, 
							AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // constructor with an additional confining potential
  //
  // nbrSiteX = number of sites in the x (or y) direction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // mus = sublattice chemical potential on A sites
  // confiningPotentialXCoordinates = x coordiantes of the confining potential
  // confiningPotentialYCoordinates = y coordiantes of the confining potential
  // confiningPotentialAmplitudes = amplitudes of the confining potential on each sites
  // nbrConfiningPotentials = number of sites where there the confining potential has a non-zero amplitude
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelCheckerboardLatticeFullOBCC4Symmetry(int nbrSiteX, double t1, double t2, double mus,
							int* confiningPotentialXCoordinates, int* confiningPotentialYCoordinates, 
							double* confiningPotentialAmplitudes, int nbrConfiningPotentials,
							AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // destructor
  //
  ~TightBindingModelCheckerboardLatticeFullOBCC4Symmetry();

  // get the index of the real space tight binding model from the real space coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndex(int x, int y);
  
  // get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // return value = linearized index  (negative if non valid)
  virtual int GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y);

  // get the real space coordinates from the index of the real space tight binding model
  //
  // index = linearized index of the real space tight binding model
  // x = reference on the x coordinate of the unit cell
  // y = reference on the y coordinate of the unit cell
  virtual void GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y);

 protected :

  // create the linearized index map to go from the canonical linearized index to the C4 compatible one
  //
  virtual void CreateLinearizedIndexMap();

};

// get the index of the real space tight binding model from the real space coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int TightBindingModelCheckerboardLatticeFullOBCC4Symmetry::GetRealSpaceTightBindingLinearizedIndex(int x, int y)
{  
  return this->LinearizedIndexMap[y  + (x * this->NbrSiteY)]; 
}

// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// return value = linearized index (negative if non valid)

inline int TightBindingModelCheckerboardLatticeFullOBCC4Symmetry::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y)
{
  if ((x < 0) || (x >= this->NbrSiteX))
    {
      return -1;
    }
  if ((y < 0) || (y >= this->NbrSiteY))
    {
      return -1;
    }
  return this->GetRealSpaceTightBindingLinearizedIndex(x, y); 
}

// get the real space coordinates from the index of the real space tight binding model
//
// index = linearized index of the real space tight binding model
// x = reference on the x coordinate of the unit cell
// y = reference on the y coordinate of the unit cell

inline void TightBindingModelCheckerboardLatticeFullOBCC4Symmetry::GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y)
{
  index = this->InvertLinearizedIndexMap[index];
  y = index % this->NbrSiteY;
  x = index / this->NbrSiteY;
}


#endif
