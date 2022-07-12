////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of abstract tight binding model                    //
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


#ifndef ABSTRACTTIGHTBINDINGMODEL_H
#define ABSTRACTTIGHTBINDINGMODEL_H


#include "config.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>


using std::ostream;


class AbstractTightBindingModel
{

 friend class FTIComputeBandStructureOperation;

 protected:

  // number of bands (including any internal degree of freedom such as spin)
  int NbrBands;

  // number of states per band
  long NbrStatePerBand;

  // pointer to the architecture
  AbstractArchitecture* Architecture;

  //  array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
  int* NbrConnectedOrbitals;
  // array that gives the orbital indices of the connected orbitals
  int** ConnectedOrbitalIndices;
  // array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
  int** ConnectedOrbitalSpatialIndices;
  // array that gives the hopping amplitudes for each pair of connected orbitals
  Complex** ConnectedOrbitalHoppingAmplitudes;

  // One body eigenstate basis associated to each point of the band structure, the array index corresponds to the linearized momentum
  ComplexMatrix* OneBodyBasis;

  // energy spectrum of the band structure, first index is the band index, second index is the linearized momentum
  double** EnergyBandStructure;

 
 public:

  // default constructor
  //
  AbstractTightBindingModel();

  // destructor
  //
  virtual ~AbstractTightBindingModel();

  // get the number of bands
  //
  // return value =  number of bands
  virtual int GetNbrBands();

  // get the number of unit cells
  //
  // return value =  number of unit cells
  virtual long GetNbrCells() {return this->GetNbrStatePerBand();}

  // get the number of sites
  //
  // return value = number of sites in the simulation cell
  virtual long GetNbrSites() {return this->GetNbrStatePerBand() * this->GetNbrBands();}

  // get the number of states per band
  //
  // return value =  number of states per band
  virtual long GetNbrStatePerBand();

  // get the size (length / area / volume ... ) of the unit cell
  //
  // return value =  size
  virtual double GetUnitCellSize() = 0;

  // get the energy at a given momentum of the band structure
  //
  // bandIndex = index of the band to consider
  // momentumIndex = linearized momentum
  // return value = energy
  virtual double GetEnergy(int bandIndex, int momentumIndex) = 0;

  // ask if the one body transformation matrices are available
  //
  // return value = true if the one body transformation matrices are available
  virtual bool HaveOneBodyBasis() = 0;

  // get  the one body transformation matrix corresponding to a given momentum of the band structure
  //
  // momentumIndex = linearized momentum
  // return value = reference on the one body transformation matrix
  virtual ComplexMatrix& GetOneBodyMatrix(int momentumIndex) = 0;

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
  // dimensionIdx = index of lattice dimension, labeled from 0, ..., d-1
  virtual void GetLatticeVector(RealVector &position, int dimensionIdx = 0);

  // get the reciprocal lattice vector for the n-th fundamental lattice direction
  //
  // latticeVector[out] = reference on a vector where the answer is supplied
  // dimensionIdx = index of lattice dimension, labeled from 0, ..., d-1
  virtual void GetReciprocalLatticeVector(RealVector &position, int dimensionIdx = 0);

  // get the lattice vector for translation along the fundamental lattice directions
  //
  // latticeVector[out] = reference on a vector where the answer is supplied
  // numTranslations = coordinates in terms of lattice vectors
  // sublatticeIndex = index of the sub-lattice position
  virtual void GetSitePosition(RealVector &position, RealVector &numTranslations, int sublatticeIndex);

  // get the tight binding hamiltonian in real space 
  // 
  // return value = tight binding hamiltonian
  virtual HermitianMatrix GetRealSpaceTightBindingHamiltonian();

  // get the tight binding hamiltonian in real space , without assuming its hermiticity
  // 
  // return value = tight binding hamiltonian
  virtual ComplexMatrix GetRealSpaceTightBindingNonHermitianHamiltonian();

  // test the tight binding hamiltonian in real space, checking its hermiticity
  // 
  // return value = true if the tight binding hamiltonian is hermitian
  bool TestRealSpaceTightBindingHamiltonian();

  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

  // write the full band structure information in a binary file
  //
  // fileName = name of the output file 
  // return value = true if no error occured  
  virtual bool WriteBandStructure(char* fileName);

  // write the full band structure information in an ASCII file
  //
  // fileName = name of the output file 
  // return value = true if no error occured  
  virtual bool WriteBandStructureASCII(char* fileName);

  // compute the spread of a given band
  // 
  // band = band index
  // return value = band spread 
  virtual double ComputeBandSpread(int band);

  // compute the direct gap between two bands
  // 
  // band1 = index of the lower band 
  // band2 = index of the upper band (-1 if it has to be set to band1 + 1)
  // return value =  direct band gap
  virtual double ComputeDirectBandGap(int band1, int band2 = -1);
  
  // compute the indirect gap between two bands
  // 
  // band1 = index of the lower band 
  // band2 = index of the upper band (-1 if it has to be set to band1 + 1)
  // return value =  direct band gap
  virtual double ComputeIndirectBandGap(int band1, int band2 = -1);
  
  // compute the ground state energy for a number of fermions filling the band 
  // 
  // nbrFermions = number of particles in the band structure
  // bands = number of bands used in groundstate configuration
  // return value =  total groundstate energy
  virtual double ComputeGroundstateEnergy(int nbrFermions, int &bands, bool verbose = false);
  
  // return the energy of the lowest energy single-particle state
  // 
  virtual double SingleParticleGroundstateEnergy();
  
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

 protected:

  // write a header that describes the tight binding model
  // 
  // output = reference on the output stream
  // return value  = reference on the output stream
  virtual ofstream& WriteHeader(ofstream& output);
  
  // read the header that describes the tight binding model
  // 
  // return value = size of header that was read (negative if unsuccessful)
  virtual int ReadHeader(ifstream& input);

  // write the eigenvalues and eigenstates (if available) to a band structure file
  //
  // output = reference on the output stream
  // return value  = reference on the output stream
  virtual ostream& WriteEigensystem(ofstream& output);

  // read the eigenvalues and eigenstates from a band structure file
  //
  // input = input stream from which header was read
  // HeaderSize = size of header that was read
  // return value  = size of band structure that was
  virtual void ReadEigensystem(ifstream& input, int HeaderSize, unsigned long FileSize = 0);

  // write an ASCII header that describes the tight binding model
  // 
  // output = reference on the output stream
  // commentChar = optional ASCII character to add in front of each header line
  // return value  = reference on the output stream
  virtual ostream& WriteASCIIHeader(ostream& output, char commentChar = '\0');

  // compute the band structure
  //
  virtual void ComputeBandStructure();

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

};

// get the number of bands
//
// return value =  number of bands

inline int AbstractTightBindingModel::GetNbrBands()
{
  return this->NbrBands;
}

// get the number of states per band
//
// return value =  number of states per band

inline long AbstractTightBindingModel::GetNbrStatePerBand()
{
  return this->NbrStatePerBand;
}

#endif
