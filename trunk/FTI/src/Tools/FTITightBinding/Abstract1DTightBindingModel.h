////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 1D tight binding model                  //
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


#ifndef ABSTRACT1DTIGHTBINDINGMODEL_H
#define ABSTRACT1DTIGHTBINDINGMODEL_H


#include "config.h"
#include "Tools/FTITightBinding/AbstractTightBindingModel.h"


class Abstract1DTightBindingModel : public AbstractTightBindingModel
{

 protected:

   // number of sites in the x direction
  int NbrSiteX;

  // numerical factor for momentum along x
  double KxFactor;

  // boundary condition twisting angle along x
  double GammaX;

  // lattice vector along the 1-direction
  RealVector LatticeVector1;

  // embedding of sublattices relative to the unit cell reference point along x
  RealVector EmbeddingX;

 public:

  // default constructor
  //
  Abstract1DTightBindingModel();

  // destructor
  //
  ~Abstract1DTightBindingModel();

  // get the linearized momentum index
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // return value = linearized momentum index
  virtual int GetLinearizedMomentumIndex(int kx, int ky);

  // get momentum value from a linearized momentum index
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // return value = linearized momentum index
  virtual void GetLinearizedMomentumIndex(int index, int& kx, int& ky);

  // get the linearized momentum index
  //
  // kx = momentum along the x direction
  // return value = linearized momentum index
  virtual int GetLinearizedMomentumIndex(int kx);

  // get momentum value from a linearized momentum index
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // return value = linearized momentum index
  virtual void GetMomentumFromLinearizedIndex(int index, int& kx);

  // get the linearized momentum index, without assuming k to be in the first BZ
  //
  // kx = momentum along the x direction
  // return value = linearized momentum index
  virtual int GetLinearizedMomentumIndexSafe(int kx);

  // get momentum value from a linearized momentum index, without assuming k to be in the first BZ
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // return value = inearized momentum index
  virtual void GetMomentumFromLinearizedIndexSafe(int index, int& kx);

  // get the index of the real space tight binding model from the real space coordinates
  //
  // x = x coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndex(int x, int orbitalIndex);
  
  // get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
  //
  // x = x coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndexSafe(int x, int orbitalIndex);

  // get the real space coordinates from the index of the real space tight binding model
  //
  // index = linearized index of the real space tight binding model
  // x = reference on the x coordinate of the unit cell
  // orbitalIndex = reference on the index of the orbital / site within the unit cell
  virtual void GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& orbitalIndex);

  // compute the index in real space lattice starting from the cartesian coordinates
  //
  // i = cartesian coordinate in the x direction of the Bravais lattice
  // p = reference on the first lattice index
  virtual void GetRealSpaceIndex (int i, int& p);
  
  // get the energy at a given momentum of the band structure
  //
  // bandIndex = index of the band to consider
  // momentumIndex = linearized momentum
  // return value = energy
  virtual double GetEnergy(int bandIndex, int momentumIndex);

  // get all the energies of a band, sorted from the smallest to the largest
  //
  // energies = reference to the array where the energies will be stored (the allocation is done by the method)
  // momenta = reference to the array where the linearized momenta associated to each energy will be stored (the allocation is done by the method)
  // bandIndex = index of the band to consider
  virtual void GetEnergies(double*& energies, int*& momenta, int bandIndex);

  // get all the energies, sorted from the smallest to the largest
  //
  // energies = reference to the array where the energies will be stored (the allocation is done by the method)
  // momenta = reference to the array where the linearized momentum associated to each energy will be stored (the allocation is done by the method)
  // bandIndices = reference to the array where the band index associated to each energy will be stored (the allocation is done by the method)
  virtual void GetAllEnergies(double*& energies, int*& momenta, int*& bandIndices);

  // Computes value of projected momentum along the lattice directions
  //
  // kx = first coordinate of the given point in the Brillouin zone
  // return value = projected momentum
  virtual double GetProjectedMomentum(int kx);
  
  // ask if the one body transformation matrices are available
  //
  // return value = true if the one body transformation matrices are available
  virtual bool HaveOneBodyBasis();

  // get  the one body transformation matrix corresponding to a given momentum of the band structure
  //
  // momentumIndex = linearized momentum
  // return value = reference on the one body transformation matrix
  virtual ComplexMatrix& GetOneBodyMatrix(int momentumIndex);
  
  // get the number of sites in the x direction
  //
  // return value = number of sites in the x direction
  int GetNbrSiteX();

  // get the position of a sublattice site
  //
  // position = reference on a vector where the answer is supplied
  // sublatticeIndex = index of the sub-lattice position
  virtual void GetSublatticeVector(RealVector &position, int sublatticeIndex);

  // get the lattice vector for translation along the fundamental lattice directions
  //
  // latticeVector[out] = reference on a vector where the answer is supplied
  // numTranslations = vector of the number of translations along each lattice direction, in units of unit cell size
  virtual void GetLatticeVector(RealVector &position, RealVector &numTranslations);

  // get the elementary lattice vector for translation along the n-th fundamental lattice directions
  //
  // latticeVector[out] = reference on a vector where the answer is supplied
  // dimensionIdx = index of lattice dimension, labelled from 0, ..., d-1
  virtual void GetLatticeVector(RealVector &position, int dimensionIdx = 0);

  // get the reciprocal lattice vector for the n-th fundamental lattice direction
  //
  // latticeVector[out] = reference on a vector where the answer is supplied
  // dimensionIdx = index of lattice dimension, labeled from 0, ..., d-1
  virtual void GetReciprocalLatticeVector(RealVector &position, int dimensionIdx = 0);

  // get the size (length / area / volume ... ) of the unit cell
  //
  // return value =  size
  virtual double GetUnitCellSize();

  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

  // flatten the energy spectrum
  //
  // flatteningFactor = flattening factor applied to each band (the band is rescale with respect to its average value)
  // bandShift = shift each band by bandShift times the band index
  virtual void FlattenBands(double flatteningFactor, double bandShift);

  // evaluate the two point correlation function in a given region
  //
  // maxX = x coordinate of the region upper right corner 
  // maxY = y coordinate of the region upper right corner 
  // occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
  // nbrOccupiedMomenta = number of occupied momenta
  // bandIndex = index of the band to consider
  // return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)
  virtual HermitianMatrix EvaluateFullTwoPointCorrelationFunction(int maxX, int maxY, int* occupiedMomenta, int nbrOccupiedMomenta, int bandIndex);

  // evaluate the two point correlation function in a given region
  //
  // maxX = x coordinate of the region upper right corner 
  // maxY = y coordinate of the region upper right corner 
  // occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
  // bandIndices = array that gives the band index of each occupied state
  // nbrOccupiedMomenta = number of occupied momenta
  // return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)
  virtual HermitianMatrix EvaluateFullTwoPointCorrelationFunction(int maxX, int maxY, int* occupiedMomenta, int* bandIndices, int nbrOccupiedMomenta);

  // evaluate the mixed two point correlation function in a given region, assuming translation invariance along one direction
  //
  // maxX = length along the borken translation direction of the region 
  // ky = momentum along the translation invariant direction
  // occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
  // bandIndices = array that gives the band index of each occupied state
  // nbrOccupiedMomenta = number of occupied momenta
  // return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)
  virtual HermitianMatrix EvaluateFullMixedTwoPointCorrelationFunctionWithK(int maxX, int ky, int* occupiedMomenta, int* bandIndices, int nbrOccupiedMomenta);

  // compute the band structure at a single point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // energies = array where the energies will be stored
  virtual void ComputeBandStructureSinglePoint(double kx, double* energies);

  // compute the Bloch hamiltonian at a point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // return value = Bloch hamiltonian
  virtual HermitianMatrix ComputeBlochHamiltonian(double kx);

 protected:

  // write an header that describes the tight binding model
  // 
  // output = reference on the output stream
  // return value  = reference on the output stream
  virtual ofstream& WriteHeader(ofstream& output);

  // read the header that describes the tight binding model
  // 
  // return value = size of header that was read
  virtual int ReadHeader(ifstream& input);

  // read the eigenvalues and eigenstates from a band structure file
  // input = input stream from which header was read
  // HeaderSize = size of header that was read
  // return value  = size of band structure that was
  virtual void ReadEigensystem(ifstream& input, int HeaderSize, unsigned long FileSize = 0);

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
  virtual ComplexMatrix BuildTightBindingNonHermitianHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, 
									  Complex** hoppingAmplitudes);

  // build the tight binding hamiltonian in recirpocal space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
  //
  // kx = momentum along the x direction (in 2pi /N_x unit) for which the hamiltonian in recirpocal space has to be computed
  // nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
  // orbitalIndices = array that gives the orbital indices of the connected orbitals
  // spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
  // hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
  // return value = tight binding hamiltonian in real space 
  virtual HermitianMatrix BuildTightBindingHamiltonianReciprocalSpace(int kx, int* nbrConnectedOrbitals, int** orbitalIndices, 
								      int** spatialIndices, Complex** hoppingAmplitudes);

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};

// get the energy at a given momentum of the band structure
//
// bandIndex = index of the band to consider
// momentumIndex = linearized momentum
// return value = energy

inline double Abstract1DTightBindingModel::GetEnergy(int bandIndex, int momentumIndex)
{
  return this->EnergyBandStructure[bandIndex][momentumIndex];
}

// get  the one body transformation matrix corresponding to a given momentum of the band structure
//
// momentumIndex = linearized momentum
// return value = reference on the one body transformation matrix

inline ComplexMatrix&  Abstract1DTightBindingModel::GetOneBodyMatrix(int momentumIndex)
{
  if (this->OneBodyBasis!=NULL)
    return this->OneBodyBasis[momentumIndex];
  else
    {
      std::cerr << "Error: requested one-body matrix, but TightBindingModel did not store the single-particle eigenstates."<<std::endl;
      exit(1);
    }
}

// get the number of sites in the x direction
//
// return value = number of sites in the x direction

inline int Abstract1DTightBindingModel::GetNbrSiteX()
{
  return this->NbrSiteX;
}

// ask if the one body transformation matrices are available
//
// return value = true if the one body transformation matrices are available

inline bool Abstract1DTightBindingModel::HaveOneBodyBasis()
{
  return (this->OneBodyBasis != 0);
}

// get the linearized momentum index
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// return value = linearized momentum index

inline int Abstract1DTightBindingModel::GetLinearizedMomentumIndex(int kx, int ky)
{
  return kx;
}

// get momentum value from a linearized momentum index
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// return value = linearized momentum index

inline void Abstract1DTightBindingModel::GetLinearizedMomentumIndex(int index, int& kx, int& ky)
{
  kx = index;
  ky = 0;
}

// get the linearized momentum index
//
// kx = momentum along the x direction
// return value = linearized momentum index

inline int Abstract1DTightBindingModel::GetLinearizedMomentumIndex(int kx)
{
  return kx;
}

// get the linearized momentum index, without assuming k to be in the first BZ
//
// kx = momentum along the x direction
// return value = linearized momentum index

inline int Abstract1DTightBindingModel::GetLinearizedMomentumIndexSafe(int kx)
{
  return kx;
}

// get momentum value from a linearized momentum index, without assuming k to be in the first BZ
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// return value = inearized momentum index

inline void  Abstract1DTightBindingModel::GetMomentumFromLinearizedIndex(int index, int& kx)
{
  kx = index; 
}

// get momentum value from a linearized momentum index, without assuming k to be in the first BZ
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// return value = inearized momentum index

inline void Abstract1DTightBindingModel::GetMomentumFromLinearizedIndexSafe(int index, int& kx)
{
  index %= this->NbrSiteX;
  if (index < 0)
    index +=  this->NbrSiteX;
  kx = index; 
}

// get the index of the real space tight binding model from the real space coordinates
//
// x = x coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int  Abstract1DTightBindingModel::GetRealSpaceTightBindingLinearizedIndex(int x, int orbitalIndex)
{
  return (orbitalIndex + (x * this->NbrBands)); 
}
  
// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int Abstract1DTightBindingModel::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int orbitalIndex)
{
  orbitalIndex %= this->NbrBands;
  if (orbitalIndex < 0)
    orbitalIndex +=  this->NbrBands;
  x %= this->NbrSiteX;
  if (x < 0)
    x +=  this->NbrSiteX;
  return (orbitalIndex + (x * this->NbrBands)); 
}

// get the real space coordinates from the index of the real space tight binding model
//
// index = linearized index of the real space tight binding model
// x = reference on the x coordinate of the unit cell
// orbitalIndex = reference on the index of the orbital / site within the unit cell

inline void Abstract1DTightBindingModel::GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& orbitalIndex)
{
  orbitalIndex = index % this->NbrBands;
  index /= this->NbrBands;
  x = index;
}

// Computes value of projected momentum along the lattice directions
//
// kx = first coordinate of the given point in the Brillouin zone
// return value = projected momentum

inline double Abstract1DTightBindingModel::GetProjectedMomentum(int kx)
{
  return (this->KxFactor * ((double) kx));
}
  
// compute the index in real space lattice starting from the cartesian coordinates
//
// i = cartesian coordinate in the x direction of the Bravais lattice
// p = reference on the first lattice index

inline void Abstract1DTightBindingModel::GetRealSpaceIndex (int i, int& p)
{
  p = i;
}
  
#endif
