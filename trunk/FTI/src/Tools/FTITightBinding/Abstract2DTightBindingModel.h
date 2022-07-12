////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 2D tight binding model                  //
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


#ifndef ABSTRACT2DTIGHTBINDINGMODEL_H
#define ABSTRACT2DTIGHTBINDINGMODEL_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract1DTightBindingModel.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"



class Abstract2DTightBindingModel : public Abstract1DTightBindingModel
{

 protected:

   // number of sites in the y direction
  int NbrSiteY;

  // numerical factor for momentum along y
  double KyFactor;

  // boundary condition twisting angle along y
  double GammaY;

  // lattice vector along the 2-direction
  RealVector LatticeVector2;

  // embedding of sublattices relative to the unit cell reference point along y
  RealVector EmbeddingY;

  // angle between the two primitive vectors
  double TwistAngle;

  // Chern number of each band
  int* Chern;

  // Berry curvature of each band over the BZ
  RealMatrix* Curvature;

  double *LLLGammaX;
  double *LLLGammaY;
    
  //first coordinate of the first spanning vector for a tilted lattice
  int Nx1;
  //second coordinate of the first spanning vector for a tilted lattice
  int Ny1;
  //first coordinate of the second spanning vector for a tilted lattice
  int Nx2;
  //second coordinate of the second spanning vector for a tilted lattice
  int Ny2;
  //array of projected momenta
  double** ProjectedMomenta;
  //second coordinate in momentum space of the second spanning vector of the reciprocal lattice for a tilted lattice
  int Offset;
  //second coordinate in real space of the second spanning vector of the direct lattice for a tilted lattice
  int OffsetReal;

  // unitary & involutory matrix where the (a, b) element Iab appears in: I|x,y,a> = sum_b Iab |-x-d_ax,-y-d_ay,b>
  ComplexMatrix Inversion;

 public:

  // default constructor
  //
  Abstract2DTightBindingModel();

  // destructor
  //
  ~Abstract2DTightBindingModel();

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
  void GetLatticeVector(RealVector &position, int dimensionIdx);

  // get the reciprocal lattice vector along the n-th fundamental lattice direction
  //
  // latticeVector[out] = reference on a vector where the answer is supplied
  // dimensionIdx = index of lattice dimension, labeled from 0, ..., d-1
  virtual void GetReciprocalLatticeVector(RealVector &position, int dimensionIdx = 0);

  // get the size (length / area / volume ... ) of the unit cell
  //
  // return value =  size
  virtual double GetUnitCellSize();

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

  // get the linearized momentum index, without assuming k to be in the first BZ
  //
  // kx = momentum along the x direction
  // ky = momentum along the y direction
  // return value = linearized momentum index
  virtual int GetLinearizedMomentumIndexSafe(int kx, int ky);

  // get momentum value from a linearized momentum index, without assuming k to be in the first BZ
  //
  // index = linearized momentum index
  // kx = reference on the momentum along the x direction
  // ky = reference on the momentum along the y direction
  // return value = inearized momentum index
  virtual void GetLinearizedMomentumIndexSafe(int index, int& kx, int& ky);
  
  // get the index of the real space tight binding model from the real space coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndex(int x, int y, int orbitalIndex);
  
  // get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
  //
  // x = x coordinate of the unit cell
  // y = y coordinate of the unit cell
  // orbitalIndex = index of the orbital / site within the unit cell
  // return value = linearized index  
  virtual int GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int orbitalIndex);

  // get the real space coordinates from the index of the real space tight binding model
  //
  // index = linearized index of the real space tight binding model
  // x = reference on the x coordinate of the unit cell
  // y = reference on the y coordinate of the unit cell
  // orbitalIndex = reference on the index of the orbital / site within the unit cell
  virtual void GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y, int& orbitalIndex);

  // get the angle between the two primitive lattice vectors
  //
  // return value = angle between the two primitive lattice vectors
  virtual double GetTwistAngle();

  // get the number of sites in the y direction
  //
  // return value = number of sites in the y direction
  int GetNbrSiteY();

  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);
  
  // write the energy spectrum in an ASCII file in a single column
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrumColumn(char* fileName);

  // write the energy spectrum in an ASCII file, focusing on lines connecting the high symmetry points
  //
  // fileName = name of the ASCII file 
  // nbrSteps = number of steps between two consecutive high symmetry points
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrumAlongHighSymmetryPoints(char* fileName, int nbrSteps);

  // write the full band structure information in an ASCII file
  //
  // fileName = name of the output file 
  // return value = true if no error occured  
  virtual bool WriteBandStructureASCII(char* fileName);

  // compute the exponentiated, unitary Abelian connection
  //
  // kx = momentum along x
  // ky = momentum along y
  // qx = momentum transfer along x
  // qy = momentum transfer along y
  // band = band index
  // return value = Phase of < u(k) | u(k+q) >
  virtual Complex GetAbelianConnection(int kx, int ky, int qx, int qy, int band);

  // compute the exponentiated, unitary Abelian connection times the quantum distance
  //
  // kx = momentum along x
  // ky = momentum along y
  // qx = momentum transfer along x
  // qy = momentum transfer along y
  // band = band index
  // return value = < u(k) | u(k+q) >
  virtual Complex GetAbelianConnectionQuantumDistance(int kx, int ky, int qx, int qy, int band);

  // compute the unitary Abelian Wilson loop
  //
  // ky = momentum along y
  // band = band index
  // return value = value of the Wilson loop
  virtual Complex GetAbelianWilsonLoopX(int ky, int band);

  // compute the unitary Abelian Wilson loop
  //
  // kx = momentum along x
  // band = band index
  // return value = value of the Wilson loop
  virtual Complex GetAbelianWilsonLoopY(int kx, int band);

  // get the total curvature over the plaquette specified by the lower-left corner. C > 0 corresponds to Bz < 0 magnetic field (holomorphic z = x + iy)
  //
  // kx = momentum along x
  // ky = momentum along y
  // band = band index
  // return value = Arg(<k|k+x><k+x|k+x+y><k+x+y|k+y><k+y|k>) / (2 Pi), range = (-0.5, 0.5]
  virtual double GetCurvature(int kx, int ky, int band);

  // get the Chern number of a specified band
  //
  // band = band index
  // return = Chern number
  virtual int GetChernNumber(int band);

  // compute the Chern number of several bands
  //
  // bands = band indices
  // nbrBands = number of bands that have to be taken into account
  // return value = Chern number
  virtual double ComputeChernNumber(int* bands, int nbrBands);

  // get the LLL twist angle along x for Bloch construction gauge fixing
  //
  // band = band index
  // return = gamma_x
  virtual double GetLLLGammaX(int band);

  // get the LLL twist angle along y for Bloch construction gauge fixing
  //
  // band = band index
  // return = gamma_y
  virtual double GetLLLGammaY(int band);

  // compute the total curvature over the plaquette specified by the lower-left corner. C > 0 corresponds to Bz < 0 magnetic field (holomorphic z = x + iy)
  //
  // kx = momentum along x
  // ky = momentum along y
  // band = band index
  // return value = Arg(<k|k+x><k+x|k+x+y><k+x+y|k+y><k+y|k>) / (2 Pi), range = (-0.5, 0.5]
  virtual double ComputeCurvatureSinglePlaquette(int kx, int ky, int band);

  // compute the stream function for the part of Berry connections that accounts for curvature fluctuations
  //
  // band = band index
  // phi = reference to the stream function over linearized BZ
  // vx = reference to the linear coefficent in phi along ky with a minus sign
  // vy = reference to the linear coefficent in phi along kx
  // return = 0 if succeed, otherwise fail
  virtual int ComputeStreamFunction(int band, RealVector& phi, double& vx, double& vy);

  // build the gauge transform such that gauge(k) * |k⟩_lat is in the "Γ"-shaped parallel-transport gauge
  // 
  // band = band index
  // gauge = reference to the gauge transform
  virtual void BuildParallelTransportGauge(int band, ComplexMatrix& gauge);

  // build the gauge transform such that gauge(k) * |k⟩_lat is in the generalized π/2-rotated Landau gauge
  //
  // band = band index
  // gauge = reference to the gauge transform
  // return = 0 if succeed, otherwise fail
  virtual int BuildGeneralizedLandauGauge(int band, ComplexMatrix& gauge);

  // compute the curvature over each plaquette in the BZ, and also Chern number
  //
  // band = band index
  // return = 0 if succeed, otherwise fail
  virtual void ComputeCurvature();

  // compute the Chern number of a given band
  //
  // band = band index
  // return value = Chern number
  virtual double ComputeChernNumber(int band);

  // compute the Berry curvature  of a given band
  //
  // band = band index
  // fileName = name of the output file 
  // return value = Chern number
  virtual double ComputeBerryCurvature(int band, char* fileName);
  
  
  //compute the complex eigenvalues of the D(ky) matrix (in order to compute the Z2 invariant)
  //
  //bandIndex = band index (corresponds to two bands that are related by time reversal symmetry)
  //nbrOccupiedBands = dimension of the D matrix
  //DMatrixEigenvalues = array of eigenvalues of the D Matrix, for all values of ky
  //kyMin = minimal value of ky for which the D matrix has to be diagonalized
  //kyMax = maximal value of ky for which the D matrix has to be diagonalized
  //nbrKy = number of ky values for which the D matrix has to be diagonalized
  //return value = array of eigenvalues of the D Matrix
  virtual Complex** ComputeDMatrixEigenvalues(int nbrOccupiedBands, int kyMin, int kyMax, int nbrKy);
  
  // write the eigenvalues of the D matrix in an ASCII file
  //
  // fileName = name of the ASCII file 
  //nbrOccupiedBands = nbr of occupied bands
  // return value = true if no error occured
  virtual bool WriteAsciiDMatrixEigenValues(char* fileName, int nbrOccupiedBands);
  
  //compute the Z2 topological invariant for a system with time reversal symmetry
  //
  //nbrOccupiedBands = number of occupied bands
  //return value = Z2 invariant
  virtual int ComputeZ2Invariant(int nbrOccupiedBands);
  
  //Computes value of projected momentum along the lattice directions
  //
  //kx = first coordinate of the given point in the Brillouin zone
  //ky = second coordinate of the given point in the Brillouin zone
  //latticeComponent = index of the lattice vector along which the projection is to be performed
  //return value = projected momentum
  virtual double GetProjectedMomentum(int kx, int ky, int latticeComponent);
  
  //get the value of the embedding vectors
  //
  virtual void GetEmbedding(RealVector& embeddingX, RealVector& embeddingY);

  //get the value of the embedding for a given sublattice
  //
  virtual void GetEmbedding(int sublattice, double &embeddingX, double &embeddingY);

  // get phase for embedding a given sublattice for the given momenta
  double GetEmbeddingPhase(int subl, double K1, double K2) {return K1*this->EmbeddingX[subl] + K2*this->EmbeddingY[subl];}
  
  //set the value of the embedding vectors
  //
  virtual void SetEmbedding(RealVector embeddingX, RealVector embeddingY);

  //set the value of the embedding vectors to zero
  //
  void SetNoEmbedding();
  
  //set the value of the embedding vectors from an external ascii file
  //
  //embeddingFileName = name of the ascii file that defines the embedding
  //
  bool SetEmbeddingFromAsciiFile(char* embeddingFileName);

  // code set of quantum numbers posx, posy into a single integer
  // posx = position along x-direction
  // posy = position along y-direction
  // numXTranslations = number of translation in the x direction to get back to the unit cell 
  // numXTranslations = number of translation in the y direction to get back to the unit cell
  //
  virtual int EncodeSublatticeIndex(int posx, int posy, int & numXTranslations, int &numYTranslations, Complex &translationPhase);

  // decode single integer for sublattice index into set of quantum numbers/positions posx, posy
  // index = sublattice index
  // [out] posx = position along x-direction
  // [out] posy = position along y-direction
  //
  virtual void DecodeSublatticeIndex(int index, int &posx, int &posy);

  // obtain dimensions of magnetic unit cell for case of Hofstadter model
  // numX = number of unit cells within MUC along x-direction
  // numY = number of unit cells within MUC along y-direction
  virtual void GetMUCDimensions(int &numX, int &numY);
  
  // compute the index in real space lattice starting from the cartesian coordinates
  //
  // i = cartesian coordinate in the x direction of the Bravais lattice
  // j = cartesian coordinate in the y direction of the Bravais lattice
  // p = reference on the first lattice index
  // q = reference on the second lattice index
  virtual void GetRealSpaceIndex (int i, int j, int& p, int& q);


  // generate a tight-binding Density-Density interaction in real space for the current Tight-Binding Model, given neighbourship relations of orbitals
  // nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  // interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  // interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  // interactingOrbitalsPotentials = intensity of each density-density term 
  RealSymmetricMatrix GenerateDensityDensityInteraction(int *NbrInteractingOrbitals, int **InteractingOrbitalsOrbitalIndices, int **InteractingOrbitalsSpatialIndices,  double **InteractingOrbitalsPotentials);
  
  // compute the band structure at a single point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // ky = momentum along the x axis
  // energies = array where the energies will be stored
  virtual void ComputeBandStructureSinglePoint(double kx, double ky, double* energies);

  // compute the Bloch hamiltonian at a point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // ky = momentum along the x axis
  // return value = Bloch hamiltonian
  virtual HermitianMatrix ComputeBlochHamiltonian(double kx, double ky);

  // get the high symmetry points 
  //
  // pointNames = name of each high symmetry point
  // pointCoordinates = coordinates in the first Brillouin zone of the high symmetry points
  // return value = number of high symmetry points
  virtual int GetHighSymmetryPoints(char**& pointNames, double**& pointCoordinates);
  
  // compute the distance between two points in the first Brillouin zone, changing the coordinates the second one by a reciprocal lattice vector if needed
  //
  // kx1 = momentum of the first point along the x axis
  // ky1 = momentum of the first point along the y axis
  // kx2 = reference on the momentum of the second point along the x axis
  // ky2 = reference on the momentum of the second point along the y axis
  // return value = distance between the two points
  virtual double GetDistanceReciprocalSpace(double kx1, double ky1, double& kx2, double& ky2);

  // evaluate the two point correlation function 
  //
  // x = linearized position index of the first point
  // y = linearized position index of the second point
  // occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
  // nbrOccupiedMomenta = number of occupied momenta
  // bandIndex = index of the band to consider
  // return value = value of the two point correlation function 
  virtual Complex EvaluateTwoPointCorrelationFunction(int x, int y, int* occupiedMomenta, int nbrOccupiedMomenta, int bandIndex);

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
  // bandIndices = indices of the band corresponding ot each occupied state
  // nbrOccupiedStates = number of occupied states
  // return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)
  virtual HermitianMatrix EvaluateFullTwoPointCorrelationFunction(int maxX, int maxY, int* occupiedMomenta, int* bandIndices, int nbrOccupiedStates);

  // evaluate the mixed two point correlation function in a given region, assuming translation invariance along one direction
  //
  // maxX = length along the borken translation direction of the region 
  // ky = momentum along the translation invariant direction
  // occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
  // bandIndices = array that gives the band index of each occupied state
  // nbrOccupiedMomenta = number of occupied momenta
  // return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)
  virtual HermitianMatrix EvaluateFullMixedTwoPointCorrelationFunctionWithK(int maxX, int ky, int* occupiedMomenta, int* bandIndices, int nbrOccupiedMomenta);
  
   // compute the form factor for the density operator 
  // 
  // kx = momentum along x of annihilation operator
  // ky = momentum along y of creation operator
  // qx = momentum transfer along x direction
  // qy = momentum transfer along y direction
  // valleyIndex = valley index of density operator
  virtual Complex ComputeDensityFormFactor(int kx, int ky, int qx, int qy, int valleyIndex);
  
  // evaluate the norm of a momentum space vector
  //
  // kx = component of momentum along first Bravais vector
  // ky = component of momentum along second Bravais vector
  // return value = norm of vector
  virtual double EvaluateNormQ(int kx, int ky);
  
 protected:

  // write an header that describes the tight binding model
  // 
  // output = reference on the output stream
  // return value  = reference on the output stream
  virtual ofstream& WriteHeader(ofstream& output);

  // read the header that describes the tight binding model
  // 
  // return value = size of header that was read (negative if unsuccessful)
  virtual int ReadHeader(ifstream& input);
   
  //computes all the values of the momentum projected and stores them in a double array
  //
  virtual void ComputeAllProjectedMomenta();
  
  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

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
  // nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
  // orbitalIndices = array that gives the orbital indices of the connected orbitals
  // spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
  // hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
  // return value = tight binding hamiltonian in real space 
  virtual HermitianMatrix BuildTightBindingHamiltonianReciprocalSpace(int kx, int ky, int* nbrConnectedOrbitals, int** orbitalIndices, 
								      int** spatialIndices, Complex** hoppingAmplitudes);
};

// get the linearized momentum index
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// return value = linearized momentum index

inline int Abstract2DTightBindingModel::GetLinearizedMomentumIndex(int kx, int ky)
{
  return ((kx * this->NbrSiteY) + ky);
}

// get momentum value from a linearized momentum index
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// return value = inearized momentum index

inline void Abstract2DTightBindingModel::GetLinearizedMomentumIndex(int index, int& kx, int& ky)
{
  kx = index / this->NbrSiteY;
  ky = index % this->NbrSiteY;
}

// get the linearized momentum index, without assuming k to be in the first BZ
//
// kx = momentum along the x direction
// ky = momentum along the y direction
// return value = linearized momentum index

inline int Abstract2DTightBindingModel::GetLinearizedMomentumIndexSafe(int kx, int ky)
{
  while (kx < 0)
      kx += this->NbrSiteX;
  kx %= this->NbrSiteX;
  while (ky < 0)
      ky += this->NbrSiteY;
  ky %= this->NbrSiteY;
  return ((kx * this->NbrSiteY) + ky);
}

// get momentum value from a linearized momentum index, without assuming k to be in the first BZ
//
// index = linearized momentum index
// kx = reference on the momentum along the x direction
// ky = reference on the momentum along the y direction
// return value = inearized momentum index

inline void Abstract2DTightBindingModel::GetLinearizedMomentumIndexSafe(int index, int& kx, int& ky)
{
  int n = this->NbrSiteX * this->NbrSiteY;
  while (index < 0)
      index += n;
  index %= n;
  kx = index / this->NbrSiteY;
  ky = index % this->NbrSiteY;
}


// get the index of the real space tight binding model from the real space coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int Abstract2DTightBindingModel::GetRealSpaceTightBindingLinearizedIndex(int x, int y, int orbitalIndex)
{
  return (orbitalIndex + ((y  + x * this->NbrSiteY) * this->NbrBands)); 
}

// get the index of the real space tight binding model from the real space coordinates, without assumption on the coordinates
//
// x = x coordinate of the unit cell
// y = y coordinate of the unit cell
// orbitalIndex = index of the orbital / site within the unit cell
// return value = linearized index  

inline int Abstract2DTightBindingModel::GetRealSpaceTightBindingLinearizedIndexSafe(int x, int y, int orbitalIndex)
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
  return this->GetRealSpaceTightBindingLinearizedIndex(x, y, orbitalIndex); 
}

// get the real space coordinates from the index of the real space tight binding model
//
// index = linearized index of the real space tight binding model
// x = reference on the x coordinate of the unit cell
// y = reference on the y coordinate of the unit cell
// orbitalIndex = reference on the index of the orbital / site within the unit cell

inline void Abstract2DTightBindingModel::GetRealSpaceTightBindingLinearizedIndex(int index, int& x, int& y, int& orbitalIndex)
{
  orbitalIndex = index % this->NbrBands;
  index /= this->NbrBands;
  y = index % this->NbrSiteY;
  x = index / this->NbrSiteY;
}

// get the angle between the two primitive lattice vectors
//
// return value = angle between the two primitive lattice vectors

inline double Abstract2DTightBindingModel::GetTwistAngle()
{
  return this->TwistAngle;
}

// get the number of sites in the y direction
//
// return value = number of sites in the y direction

inline int Abstract2DTightBindingModel::GetNbrSiteY()
{
  return this->NbrSiteY;
}

// get the total curvature over the plaquette specified by the lower-left corner. C > 0 corresponds to Bz < 0 magnetic field (holomorphic z = x + iy)
//
// kx = momentum along x
// ky = momentum along y
// band = band index
// return value = Arg(<k|k+x><k+x|k+x+y><k+x+y|k+y><k+y|k>) / (2 Pi), range = (-0.5, 0.5]

inline double Abstract2DTightBindingModel::GetCurvature(int kx, int ky, int band)
{
    if (this->Curvature == NULL)
        this->ComputeCurvature();
    return this->Curvature[band][ky][kx];
}

// get the Chern number of a specified band
//
// band = band index
// return = Chern number

inline int Abstract2DTightBindingModel::GetChernNumber(int band)
{
    if (this->Chern == NULL)
        this->ComputeCurvature();
    return this->Chern[band];
}

// get the LLL twist angle along x for Bloch construction gauge fixing
//
// band = band index
// return = gamma_x

inline double Abstract2DTightBindingModel::GetLLLGammaX(int band)
{
    if (this->LLLGammaX == NULL)
        this->ComputeCurvature();
    return this->LLLGammaX[band];
}

// get the LLL twist angle along y for Bloch construction gauge fixing
//
// band = band index
// return = gamma_y
inline double Abstract2DTightBindingModel::GetLLLGammaY(int band)
{
    if (this->LLLGammaY == NULL)
        this->ComputeCurvature();
    return this->LLLGammaY[band];
}

// compute the curvature sum over the plaquette specified by the lower-left corner. C > 0 corresponds to Bz < 0 magnetic field (holomorphic z = x + iy)
//
// kx = momentum along x
// ky = momentum along y
// band = band index
// return value = Arg(<k|k+x><k+x|k+x+y><k+x+y|k+y><k+y|k>) / (2 Pi), range = (-0.5, 0.5]

inline double Abstract2DTightBindingModel::ComputeCurvatureSinglePlaquette(int kx, int ky, int band)
{
    Complex W = this->GetAbelianConnection(kx, ky, 1, 0, band) * this->GetAbelianConnection(kx + 1, ky, 0, 1, band);
    W *= this->GetAbelianConnection(kx + 1, ky + 1, -1, 0, band) * this->GetAbelianConnection(kx, ky + 1, 0, -1, band);
    return Arg(W) / (2 * M_PI);
}

// Computes value of projected momentum along the lattice directions
//
// kx = first coordinate of the given point in the Brillouin zone
// ky = second coordinate of the given point in the Brillouin zone
// latticeComponent = index of the lattice vector along which the projection is to be performed
// return value = projected momentum
   
  
inline double Abstract2DTightBindingModel::GetProjectedMomentum(int kx, int ky, int latticeComponent)
{
  return this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky)][latticeComponent];
}

// get the value of the embedding vectors
//

inline void Abstract2DTightBindingModel::GetEmbedding(RealVector& embeddingX, RealVector& embeddingY)
 {
   embeddingX = this->EmbeddingX;
   embeddingY = this->EmbeddingY;
 }

// get the value of the embedding for a given sublattice
//

inline void Abstract2DTightBindingModel::GetEmbedding(int sublattice, double &embeddingX, double &embeddingY)
{
  embeddingX = this->EmbeddingX[sublattice];
  embeddingY = this->EmbeddingY[sublattice];
}
 
// set the value of the embedding vectors
//

inline void Abstract2DTightBindingModel::SetEmbedding(RealVector embeddingX, RealVector embeddingY)
{
  this->EmbeddingX = embeddingX;
  this->EmbeddingY = embeddingY;
}

//set no/trivial embedding, i.e. all sublattices at the origin
//

inline void Abstract2DTightBindingModel::SetNoEmbedding()
{
  this->EmbeddingX.ResizeAndClean(NbrBands);
  this->EmbeddingY.ResizeAndClean(NbrBands);
}

 
// compute the index in real space lattice starting from the cartesian coordinates
//
// i = cartesian coordinate in the x direction of the Bravais lattice
// j = cartesian coordinate in the y direction of the Bravais lattice
// p = reference on the first lattice index
// q = reference on the second lattice index

inline void Abstract2DTightBindingModel::GetRealSpaceIndex (int i, int j, int& p, int& q)
{
  p = i - this->OffsetReal * j;
  q = j;
}

// code set of quantum numbers posx, posy into a single integer
// posx = position along x-direction
// posy = position along y-direction
// numXTranslations = number of translation in the x direction to get back to the unit cell 
// numXTranslations = number of translation in the y direction to get back to the unit cell
//

inline int  Abstract2DTightBindingModel::EncodeSublatticeIndex(int posx, int posy,int & numXTranslations,int &numYTranslations, Complex &translationPhase)
{
  std::cout <<"using dummy Abstract2DTightBindingModel::EncodeSublatticeIndex(int posx, int posy,int & numXTranslations,int &numYTranslations, Complex &translationPhase)"<<std::endl;
  return -1;
}


// decode single integer for sublattice index into set of quantum numbers/positions posx, posy
//
// index = sublattice index
// [out] posx = position along x-direction
// [out] posy = position along y-direction
//

inline void Abstract2DTightBindingModel::DecodeSublatticeIndex(int index, int &posx, int &posy)
{
  std::cout <<"using dummy Abstract2DTightBindingModel::DecodeSublatticeIndex(int index, int &posx, int &posy)"<<std::endl;
  posx = 0;
  posy = 0;
}

   
// obtain dimensions of magnetic unit cell for case of Hofstadter model
//
// numX = number of unit cells within MUC along x-direction
// numY = number of unit cells within MUC along y-direction

inline void Abstract2DTightBindingModel::GetMUCDimensions(int &numX, int &numY) 
{
  numX = 1;
  numY = 1;
}

// evaluate the norm of a momentum space vector
//
// kx = component of momentum along first Bravais vector
// ky = component of momentum along second Bravais vector
// return value = norm of vector

inline double Abstract2DTightBindingModel::EvaluateNormQ(int kx, int ky)
{
    double Kx = this->KxFactor * ((double) kx);
    double Ky = this->KyFactor * ((double) ky);
    return sqrt(Kx*Kx + Ky*Ky);
}

#endif
