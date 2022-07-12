////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 3D tight binding model                  //
//                                                                            //
//                        last modification : 01/05/2012                      //
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


#ifndef ABSTRACT3DTIGHTBINDINGMODEL_H
#define ABSTRACT3DTIGHTBINDINGMODEL_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


using std::cout;
using std::endl;


class Abstract3DTightBindingModel : public Abstract2DTightBindingModel
{

 protected:

   // number of sites in the z direction
  int NbrSiteZ;
  // number of sites in the direction perpendicular to X
  int NbrSiteYZ;

  // numerical factor for momentum along z
  double KzFactor;

  // boundary condition twisting angle along z
  double GammaZ;

  // embedding of sublattices relative to the unit cell reference point along z
  RealVector EmbeddingZ;

  // lattice vector along the third direction
  RealVector LatticeVector3;
  
  // third coordinate of the first spanning vector for a tilted lattice
  int Nz1;
  // third coordinate of the second spanning vector for a tilted lattice
  int Nz2;
  // first coordinate of the third spanning vector for a tilted lattice
  int Nx3;
  // second coordinate of the third spanning vector for a tilted lattice
  int Ny3;
  // second coordinate of the third spanning vector for a tilted lattice
  int Nz3;

  // the number of sites per unit cell  assuming the model is 2D, i.e. the actual number of orbitals times the number of orbitals in the z direction
  int NbrBandsAs2DModel;
  
  
 public:

  // default constructor
  //
  Abstract3DTightBindingModel();

  // destructor
  //
  ~Abstract3DTightBindingModel();

  // get the size (length / area / volume ... ) of the unit cell
  //
  // return value =  size
  virtual double GetUnitCellSize();

  // get the linearized momentum index
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // ky = momentum along the z direction
  // return value = inearized momentum index
  virtual int GetLinearizedMomentumIndex(int kx, int ky, int kz);

  // get momentum value from a linearized momentum index
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // kz = reference on the momentum along the z direction
  // return value = inearized momentum index
  virtual void GetLinearizedMomentumIndex(int index, int& kx, int& ky, int& kz);

  // get the linearized momentum index, without assuming k to be in the first BZ
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // ky = momentum along the z direction
  // return value = inearized momentum index
  virtual int GetLinearizedMomentumIndexSafe(int kx, int ky, int kz);

  // get momentum value from a linearized momentum index, without assuming k to be in the first BZ
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // kz = reference on the momentum along the z direction
  // return value = inearized momentum index
  virtual void GetLinearizedMomentumIndexSafe(int index, int& kx, int& ky, int& kz);

  // get the index of the real space tight binding model from the real space coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // z = z coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndex(int x, int y, int z, int orbitalIndex);
  
  // get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // z = z coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int z, int orbitalIndex);

  // get the real space coordinates from the index of the real space tight binding model
  //
  // index = linearized index of the real space tight binding model
  // x = reference on the x coordinate of the unit cell
  // y = reference on the y coordinate of the unit cell
  // z = reference on the z coordinate of the unit cell
  // orbitalIndex = reference on the index of the orbital / site within the unit cell
  virtual void GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y, int& z, int& orbitalIndex);

  //Computes value of projected momentum along the lattice directions
  //
  // kx = first coordinate of the given point in the Brillouin zone
  // ky = second coordinate of the given point in the Brillouin zone
  // kz = third coordinate of the given point in the Brillouin zone
  // latticeComponent = index of the lattice vector along which the projection is to be performed
  // return value = projected momentum
  virtual double GetProjectedMomentum(int kx, int ky, int kz, int latticeComponent);
    
  // get the index of the real space tight binding model from the real space coordinates, assuming the model is 2D
  // (the number of sites per unit cell being the actual number of orbitals times the number of orbitals in the z direction)
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndex(int x, int y, int orbitalIndex);
  
  // get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates, assuming the model is 2D
  // (the number of sites per unit cell being the actual number of orbitals times the number of orbitals in the z direction)
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int orbitalIndex);

  // get the real space coordinates from the index of the real space tight binding model, assuming the model is 2D
  // (the number of sites per unit cell being the actual number of orbitals times the number of orbitals in the z direction)
  //
  // index = linearized index of the real space tight binding model
  // x = reference on the x coordinate of the unit cell
  // y = reference on the y coordinate of the unit cell
  // orbitalIndex = reference on the index of the orbital / site within the unit cell
  //  virtual void GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y, int& orbitalIndex);

  // get the number of sites in the z direction
  //
  // return value = number of sites in the z direction
  int GetNbrSiteZ();

  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

  // compute the index in real space lattice starting from the cartesian coordinates
  //
  // i = cartesian coordinate in the x direction of the Bravais lattice
  // j = cartesian coordinate in the y direction of the Bravais lattice
  // k = cartesian coordinate in the z direction of the Bravais lattice
  // p = reference on the first lattice index
  // q = reference on the second lattice index
  // r = reference on the third lattice index
  virtual void GetRealSpaceIndex (int i, int j, int k, int& p, int& q, int& r);

  // compute the index in real space lattice starting from the cartesian coordinates
  //
  // i = cartesian coordinate in the x direction of the Bravais lattice
  // j = cartesian coordinate in the y direction of the Bravais lattice
  // k = cartesian coordinate in the z direction of the Bravais lattice
  // latticeIndices = array containing the lattice indices
  virtual void GetRealSpaceIndex (int i, int j, int k, int* latticeIndices);

  // compute the band structure at a single point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // ky = momentum along the y axis
  // kz = momentum along the z axis
  // energies = array where the energies will be stored
  virtual void ComputeBandStructureSinglePoint(double kx, double ky, double kz, double* energies);

  // compute the Bloch hamiltonian at a point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // ky = momentum along the y axis
  // kz = momentum along the z axis
  // return value = Bloch hamiltonian
  virtual HermitianMatrix ComputeBlochHamiltonian(double kx, double ky, double kz);

  protected:

  // write an header that describes the tight binding model
  // 
  // output = reference on the output stream
  // return value  = reference on the output stream
  virtual ofstream& WriteHeader(ofstream& output);

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

  //computes all the values of the momentum projected and stores them in a double array
  //
  virtual void ComputeAllProjectedMomenta();
  
  // build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
  //
  // nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
  // orbitalIndices = array that gives the orbital indices of the connected orbitals
  // spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
  // hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
  // return value = tight binding hamiltonian in real space 
  virtual HermitianMatrix BuildTightBindingHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes);

  // build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions but without assuming its hermiticiy
  //
  // nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
  // orbitalIndices = array that gives the orbital indices of the connected orbitals
  // spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
  // hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
  // return value = tight binding hamiltonian in real space 
  virtual ComplexMatrix BuildTightBindingNonHermitianHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes);

  // build the tight binding hamiltonian in recirpocal space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
  //
  // kx = momentum along the x direction (in 2pi /N_x unit) for which the hamiltonian in recirpocal space has to be computed
  // ky = momentum along the y direction (in 2pi /N_y unit) for which the hamiltonian in recirpocal space has to be computed
  // ky = momentum along the z direction (in 2pi /N_z unit) for which the hamiltonian in recirpocal space has to be computed
  // nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
  // orbitalIndices = array that gives the orbital indices of the connected orbitals
  // spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
  // hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
  // return value = tight binding hamiltonian in real space 
  virtual HermitianMatrix BuildTightBindingHamiltonianReciprocalSpace(int kx, int ky, int kz, int* nbrConnectedOrbitals, int** orbitalIndices, 
								      int** spatialIndices, Complex** hoppingAmplitudes);
  
};

// get the linearized momentum index
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// kz = momentum along the z direction
// return value = linearized momentum index

inline int Abstract3DTightBindingModel::GetLinearizedMomentumIndex(int kx, int ky, int kz)
{
  return (((kx * this->NbrSiteY) + ky) * this->NbrSiteZ + kz);
}

// get momentum value from a linearized momentum index
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// kz = reference on the momentum along the z direction
// return value = inearized momentum index

inline void Abstract3DTightBindingModel::GetLinearizedMomentumIndex(int index, int& kx, int& ky, int& kz)
{
  kx = index / this->NbrSiteYZ;
  ky = index % this->NbrSiteYZ;
  kz = ky % this->NbrSiteZ;
  ky /= this->NbrSiteZ;
}

// get the linearized momentum index, without assuming k to be in the first BZ
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// kz = momentum along the z direction
// return value = linearized momentum index

inline int Abstract3DTightBindingModel::GetLinearizedMomentumIndexSafe(int kx, int ky, int kz)
{
  while (kx < 0)
      kx += this->NbrSiteX;
  kx %= this->NbrSiteX;
  while (ky < 0)
      ky += this->NbrSiteY;
  ky %= this->NbrSiteY;
  while (kz < 0)
      kz += this->NbrSiteZ;
  kz %= this->NbrSiteZ;
  return (((kx * this->NbrSiteY) + ky) * this->NbrSiteZ + kz);
}

// get momentum value from a linearized momentum index, without assuming k to be in the first BZ
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// kz = reference on the momentum along the z direction
// return value = inearized momentum index

inline void Abstract3DTightBindingModel::GetLinearizedMomentumIndexSafe(int index, int& kx, int& ky, int& kz)
{
  int n = this->NbrSiteX * this->NbrSiteYZ;
  while (index < 0)
      index += n;
  index %= n;
  kx = index / this->NbrSiteYZ;
  ky = index % this->NbrSiteYZ;
  kz = ky % this->NbrSiteZ;
  ky /= this->NbrSiteZ;
}

// get the index of the real space tight binding model from the real space coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// z = z coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int Abstract3DTightBindingModel::GetRealSpaceTightBindingLinearizedIndex(int x, int y, int z, int orbitalIndex)
{
  return (orbitalIndex + ((((y + (x * this->NbrSiteY)) * this->NbrSiteZ) + z) * this->NbrBands)); 
}

// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// z = z coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int Abstract3DTightBindingModel::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int z, int orbitalIndex)
{
  orbitalIndex %= this->NbrBands;
  if (orbitalIndex < 0)
    orbitalIndex +=  this->NbrBands;
  x %= this->NbrSiteX;
  if (x < 0)
    x +=  this->NbrSiteX;
  y %= this->NbrSiteY;
  if (y < 0)
    y +=  this->NbrSiteY;
  z %= this->NbrSiteZ;
  if (z < 0)
    z +=  this->NbrSiteZ;
  return this->GetRealSpaceTightBindingLinearizedIndex(x, y, z, orbitalIndex); 
}

// get the real space coordinates from the index of the real space tight binding model
//
// index = linearized index of the real space tight binding model
// x = reference on the x coordinate of the unit cell
// y = reference on the y coordinate of the unit cell
// z = reference on the z coordinate of the unit cell
// orbitalIndex = reference on the index of the orbital / site within the unit cell

inline void Abstract3DTightBindingModel::GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y, int& z, int& orbitalIndex)
{
  orbitalIndex = index % this->NbrBands;
  index /= this->NbrBands;
  x = index / this->NbrSiteYZ;
  y = index % this->NbrSiteYZ;
  z = y % this->NbrSiteZ;
  y /= this->NbrSiteZ;
}

// get the number of sites in the z direction
//
// return value = number of sites in the z direction

inline int Abstract3DTightBindingModel::GetNbrSiteZ()
{
  return this->NbrSiteZ;
}

// Computes value of projected momentum along the lattice directions
//
// kx = first coordinate of the given point in the Brillouin zone
// ky = second coordinate of the given point in the Brillouin zone
// kz = third coordinate of the given point in the Brillouin zone
// latticeComponent = index of the lattice vector along which the projection is to be performed
// return value = projected momentum
   
inline double Abstract3DTightBindingModel::GetProjectedMomentum(int kx, int ky, int kz, int latticeComponent)
{
  return this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky, kz)][latticeComponent];
}

// compute the index in real space lattice starting from the cartesian coordinates
//
// i = cartesian coordinate in the x direction of the Bravais lattice
// j = cartesian coordinate in the y direction of the Bravais lattice
// k = cartesian coordinate in the z direction of the Bravais lattice
// p = reference on the first lattice index
// q = reference on the second lattice index
// r = reference on the third lattice index

inline void Abstract3DTightBindingModel::GetRealSpaceIndex (int i, int j, int k, int& p, int& q, int& r)
{
  p = i - this->OffsetReal * j;
  q = j;
  r = k;
}

// compute the index in real space lattice starting from the cartesian coordinates, assuming the model is 2D
// (the number of sites per unit cell being the actual number of orbitals times the number of orbitals in the z direction)
//
// i = cartesian coordinate in the x direction of the Bravais lattice
// j = cartesian coordinate in the y direction of the Bravais lattice
// k = cartesian coordinate in the z direction of the Bravais lattice
// latticeIndices = array containing the lattice indices

inline void Abstract3DTightBindingModel::GetRealSpaceIndex (int i, int j, int k, int* latticeIndices)
{
  latticeIndices[0] = i - this->OffsetReal * j;
  latticeIndices[1] = j;
  latticeIndices[2] = k;  
}

// get the index of the real space tight binding model from the real space coordinates, assuming the model is 2D
// (the number of sites per unit cell being the actual number of orbitals times the number of orbitals in the z direction)
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int Abstract3DTightBindingModel::GetRealSpaceTightBindingLinearizedIndex(int x, int y, int orbitalIndex)
{
  return (orbitalIndex + ((y  + x * this->NbrSiteY) * this->NbrBandsAs2DModel)); 
}

// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates, assuming the model is 2D
// (the number of sites per unit cell being the actual number of orbitals times the number of orbitals in the z direction)
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int Abstract3DTightBindingModel::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int orbitalIndex)
{
  orbitalIndex %= this->NbrBandsAs2DModel;
  if (orbitalIndex < 0)
    orbitalIndex +=  this->NbrBandsAs2DModel;
  x %= this->NbrSiteX;
  if (x < 0)
    x +=  this->NbrSiteX;
  y %= this->NbrSiteY;
  if (y < 0)
    y +=  this->NbrSiteY;
  return this->GetRealSpaceTightBindingLinearizedIndex(x, y, orbitalIndex); 
}

// get the real space coordinates from the index of the real space tight binding model, assuming the model is 2D
// (the number of sites per unit cell being the actual number of orbitals times the number of orbitals in the z direction)
//
// index = linearized index of the real space tight binding model
// x = reference on the x coordinate of the unit cell
// y = reference on the y coordinate of the unit cell
// orbitalIndex = reference on the index of the orbital / site within the unit cell

// inline void Abstract3DTightBindingModel::GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y, int& orbitalIndex)
// {
//   orbitalIndex = index % this->NbrBandsAs2DModel;
//   index /= this->NbrBandsAs2DModel;
//   y = index % this->NbrSiteY;
//   x = index / this->NbrSiteY;
// }


#endif
