////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the checkerboard lattice       //
//      with full open boundary conditions and direct implentation of the C4  //
//  symmetry implentation (WANRING: not compatible with may-body hamitonians  //
//  use TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry)         //
//                                                                            //
//                        last modification : 02/03/2018                      //
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


#ifndef TIGHTBINDINGMODELSQUARELATTICEFULLOBCWITHFULLC4SYMMETRY_H
#define TIGHTBINDINGMODELSQUARELATTICEFULLOBCWITHFULLC4SYMMETRY_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract1DTightBindingModel.h"


class TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry : public AbstractTightBindingModel
{

 protected:

   // number of sites in the y (or x) direction
  int NbrSiteY;
   // half the number of sites in the y (or x) direction
  int HalfNbrSiteY;

 public:

  // default constructor
  //
  TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry();
  
  // destructor
  //
  ~TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry();

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

  // get the index of the real space tight binding model from the real space coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // sector = reference to the additional discrete symmetry sector
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetry(int x, int y, int& sector);

  // get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // sector = reference to the additional discrete symmetry sector
  // return value = linearized index (negative if non valid)  
  virtual int GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(int x, int y, int& sector);

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

  // find the orbitals connected to those located at the origin unit cell in a given discrete symmetry sector
  // 
  virtual void FindConnectedOrbitals(int sector);

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};

// get the index of the real space tight binding model from the real space coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// return value = linearized index  

inline int TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::GetRealSpaceTightBindingLinearizedIndex(int x, int y)
{  
  return (y  + (x * this->HalfNbrSiteY)); 
}

// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// return value = linearized index (negative if non valid)

inline int TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y)
{
  if ((x < 0) || (x >= this->NbrSiteY))
    {
      return -1;
    }
  if ((y < 0) || (y >= this->NbrSiteY))
    {
      return -1;
    }
  return this->GetRealSpaceTightBindingLinearizedIndex(x, y); 
}

// get the index of the real space tight binding model from the real space coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// sector = reference to the additional discrete symmetry sector
// return value = linearized index  

inline int TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetry(int x, int y, int& sector)
{  
  if (y >= this->HalfNbrSiteY)
    {
      if (x >= this->HalfNbrSiteY)
	{
	  sector = 2;
	  y = this->NbrSiteY - 1 - y;
	  x = this->NbrSiteY - 1 - x;
	}      
      else
	{
	  int Tmp = y;
	  y = x;
	  x = this->NbrSiteY - 1 - Tmp;
	  sector = 3;
	}
    }
  else
    {
      if (x >= this->HalfNbrSiteY)
	{
	  int Tmp = y;
	  y = this->NbrSiteY - 1 - x;
	  x = Tmp;
	  sector = 1;
	}  
      else
	{
	  sector = 0;
	}    
    }
  return (y  + (x * this->HalfNbrSiteY)); 
}

// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// sector = reference to the additional discrete symmetry sector
// return value = linearized index (negative if non valid)

inline int TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(int x, int y, int& sector)
{
  if ((x < 0) || (x >= this->NbrSiteY))
    {
      return -1;
    }
  if ((y < 0) || (y >= this->NbrSiteY))
    {
      return -1;
    }
  return this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetry(x, y, sector); 
}

// get the real space coordinates from the index of the real space tight binding model
//
// index = linearized index of the real space tight binding model
// x = reference on the x coordinate of the unit cell
// y = reference on the y coordinate of the unit cell

inline void TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y)
{
  y = index % this->HalfNbrSiteY;
  x = index / this->HalfNbrSiteY;
}


// get the size (length / area / volume ... ) of the unit cell
//
// return value =  size

inline double TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::GetUnitCellSize()
{
  return 1.0;
}


// get the energy at a given momentum of the band structure
//
// bandIndex = index of the band to consider
// momentumIndex = linearized momentum
// return value = energy

inline double TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::GetEnergy(int bandIndex, int momentumIndex)
{
  return this->EnergyBandStructure[bandIndex][momentumIndex];
}

// ask if the one body transformation matrices are available
//
// return value = true if the one body transformation matrices are available

inline bool TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::HaveOneBodyBasis()
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

inline ComplexMatrix& TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::GetOneBodyMatrix(int momentumIndex)
{
  return this->OneBodyBasis[momentumIndex];
}


#endif
