////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the checkerboard lattice       //
//                      with full open boundary conditions                    //
//                                                                            //
//                        last modification : 17/02/2018                      //
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


#ifndef TIGHTBINDINGMODELCHECKERBOARDLATTICEFULLOBC_H
#define TIGHTBINDINGMODELCHECKERBOARDLATTICEFULLOBC_H


#include "config.h"
#include "Tools/FTITightBinding/AbstractTightBindingModel.h"


class TightBindingModelCheckerboardLatticeFullOBC : public AbstractTightBindingModel
{

 protected:

   // number of sites in the x direction
  int NbrSiteX;
   // number of sites in the y direction
  int NbrSiteY;

  // hoping amplitude between neareast neighbor sites
  double NNHopping;
  // hoping amplitude between next neareast neighbor sites
  double NextNNHopping;
  
  // four times the sublattice staggered chemical potential 
  double MuS;

  // linearized coordiantes of the confining potential
  int* ConfiningPotentialCoordinates;
  // amplitudes of the confining potential on each sites
  double* ConfiningPotentialAmplitudes;
  // number of sites where there the confining potential has a non-zero amplitude
  int NbrConfiningPotentials;


 public:

  // default constructor
  //
  TightBindingModelCheckerboardLatticeFullOBC();

  // basic constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction 
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // mus = sublattice chemical potential on A sites
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelCheckerboardLatticeFullOBC(int nbrSiteX, int nbrSiteY, double t1, double t2, double mus, 
					      AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // constructor with an additional confining potential
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction 
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // mus = sublattice chemical potential on A sites
  // confiningPotentialXCoordinates = x coordiantes of the confining potential
  // confiningPotentialYCoordinates = y coordiantes of the confining potential
  // confiningPotentialAmplitudes = amplitudes of the confining potential on each sites
  // nbrConfiningPotentials = number of sites where there the confining potential has a non-zero amplitude
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelCheckerboardLatticeFullOBC(int nbrSiteX, int nbrSiteY, double t1, double t2, double mus,
					      int* confiningPotentialXCoordinates, int* confiningPotentialYCoordinates, 
					      double* confiningPotentialAmplitudes, int nbrConfiningPotentials,
					      AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // destructor
  //
  ~TightBindingModelCheckerboardLatticeFullOBC();

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

  // get the size (length / area / volume ... ) of the unit cell
  //
  // return value =  size
  virtual double GetUnitCellSize();

  // get the energy at a given momentum of the band structure
  //
  // bandIndex = index of the band to consider
  // momentumIndex = linearized momentum
  // return value = energy
  virtual double GetEnergy(int bandIndex, int momentumIndex);

  // ask if the one body transformation matrices are available
  //
  // return value = true if the one body transformation matrices are available
  virtual bool HaveOneBodyBasis();

  // get  the one body transformation matrix corresponding to a given momentum of the band structure
  //
  // momentumIndex = linearized momentum
  // return value = reference on the one body transformation matrix
  virtual ComplexMatrix& GetOneBodyMatrix(int momentumIndex);

  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

 protected :

  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

};

// get the index of the real space tight binding model from the real space coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int TightBindingModelCheckerboardLatticeFullOBC::GetRealSpaceTightBindingLinearizedIndex(int x, int y)
{
  return (y  + (x * this->NbrSiteY)); 
}

// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// return value = linearized index (negative if non valid)

inline int TightBindingModelCheckerboardLatticeFullOBC::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y)
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

inline void TightBindingModelCheckerboardLatticeFullOBC::GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y)
{
  y = index % this->NbrSiteY;
  x = index / this->NbrSiteY;
}

// get the size (length / area / volume ... ) of the unit cell
//
// return value =  size

inline double TightBindingModelCheckerboardLatticeFullOBC::GetUnitCellSize()
{
  return 1.0;
}


// get the energy at a given momentum of the band structure
//
// bandIndex = index of the band to consider
// momentumIndex = linearized momentum
// return value = energy

inline double TightBindingModelCheckerboardLatticeFullOBC::GetEnergy(int bandIndex, int momentumIndex)
{
  return this->EnergyBandStructure[bandIndex][0];
}

// ask if the one body transformation matrices are available
//
// return value = true if the one body transformation matrices are available

inline bool TightBindingModelCheckerboardLatticeFullOBC::HaveOneBodyBasis()
{
  if (this->OneBodyBasis != 0)
    {
      return true;
    }
  else
    {
      return false;
    }
}

// get  the one body transformation matrix corresponding to a given momentum of the band structure
//
// momentumIndex = linearized momentum
// return value = reference on the one body transformation matrix

inline ComplexMatrix& TightBindingModelCheckerboardLatticeFullOBC::GetOneBodyMatrix(int momentumIndex)
{
  return this->OneBodyBasis[0];
}


#endif
