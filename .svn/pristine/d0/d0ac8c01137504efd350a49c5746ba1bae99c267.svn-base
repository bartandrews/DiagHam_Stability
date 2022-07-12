////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of particle on a torus                       //
//                                                                            //
//                        last modification : 18/07/2002                      //
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


#ifndef PARTICLEONTORUS_H
#define PARTICLEONTORUS_H


#include "config.h"
#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include <iostream>


using std::ostream;


class Matrix;
class SparseRealMatrix;
class SparseComplexMatrix;


class ParticleOnTorus :  public ParticleOnSphere
{

 protected:

  friend class FQHETorusApplyCNRotationOperation;
  friend class FQHETorusSymmetrizeU1U1StateOperation;

 public:

  // virtual destructor
  //
  virtual ~ParticleOnTorus ();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic() = 0;

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnTorus* targetSpace);

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphere* targetSpace);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

  // return matrix representation of the annihilation operator a_i
  //
  // i = operator index
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  virtual Matrix& A (int i, Matrix& M) = 0;

  // return matrix representation of the creation operator a^+_i
  //
  // i = operator index
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  virtual Matrix& Ad (int i, Matrix& M) = 0;

  // apply an annihilation operator a_i and return the index in the target space
  //
  // i = state index
  // m = index of annihilation operator
  // return value = index in the target space
  virtual int A (int i, int m, double &coefficient) 
  { std::cout << "Need to implement operator A" << std::endl; return this->GetHilbertSpaceDimension(); }

  // apply a creation operator a_i and return the index in the target space
  //
  // i = state index
  // m = index of annihilation operator
  // return value = index in the target space
  virtual int Ad (int i, int m, double &coefficient)
  { std::cout << "Need to implement operator Ad" << std::endl; return this->GetHilbertSpaceDimension();}

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, RealVector& groundState);
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, ComplexVector& groundState);

  // remove part of each Fock state, discarding component if the Fock state does not a given pattern
  //
  // inputVector = state to truncate
  // reducedSpace = Hilbert space where the truncated state will lie
  // pattern = array describing the pattern 
  // patternSize = pattern size
  // patternShift = indicate where the pattern has to be applied
  // return value = trucated state
  virtual RealVector TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnTorus* reducedSpace, int* pattern, int patternSize, int patternShift = 0);

  // apply the C4 rotation to a given state assumin it is an eigenstate of both kx and ky
  //
  // inputState = reference on the state that has to be rotated
  // inputSpace = Hilbert space associated to the input state
  // architecture = pointer to the architecture
  // clockwise = the rotation is done clockwise
  // return value = rotated state
  virtual ComplexVector C4Rotation (ComplexVector& inputState, ParticleOnTorus* inputSpace, bool clockwise = false, AbstractArchitecture* architecture = 0);

  // apply a magnetic translation along x to a given state
  //
  // index = state index 
  // return value = translated state index
  virtual int ApplyXMagneticTranslation(int index);

  // apply a magnetic translation along x to a given state
  //
  // index = state index 
  // sign = additional sign due to the particle statistics
  // return value = translated state index
  virtual int ApplyXMagneticTranslation(int index, double& sign);

  // transform a state expressed on a torus with a given angle to a state expressed on the same trous but a different angle
  //
  // inputVector = reference on the input vector
  // inputAngle = angle (in radian) between the two vectors that span the torus on which the input state is defined
  // inputAspectRatio = length ratio of the two vectors that span the torus on which the input state is defined
  // outputAngle = angle (in radian) between the two vectors that span the torus on which the output state is defined
  // outputAspectRatio = length ratio of the two vectors that span the torus on which the output state is defined
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  // return value = transformed state
  virtual ComplexVector ChangeTorusAngle (ComplexVector& inputVector, double inputAngle, double inputAspectRatio, double outputAngle, double outputAspectRatio,
					  unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrized a product of two uncoupled states 
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // return value = symmetrized state
  virtual RealVector SymmetrizeU1U1State (RealVector& leftVector, RealVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag = false, AbstractArchitecture* architecture = 0);
  
  // symmetrized a product of two uncoupled states 
  //
  // outputVector = reference on the vector which will contain the symmetrozed state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // return value = symmetrized state
  virtual ComplexVector SymmetrizeU1U1State (ComplexVector& leftVector, ComplexVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag = false, AbstractArchitecture* architecture = 0);
  

  // symmetrized a product of several uncoupled states 
  //
  // inputStates = states which will be symmetrized
  // inputSpaces = Hilbert spaces attached to each states
  // nbrStates = number of states to symmetrize
  // architecture = pointer to the architecture 
  // return value = symmetrized state
  virtual RealVector SymmetrizeU1U1State (RealVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, double precision = MACHINE_PRECISION,  AbstractArchitecture* architecture = 0);
  
  // symmetrize a product of several uncoupled states 
  //
  // inputStates = states which will be symmetrized
  // inputSpaces = Hilbert spaces attached to each states
  // nbrStates = number of states to symmetrize
  // architecture = pointer to the architecture 
  // return value = symmetrized state
  virtual ComplexVector SymmetrizeU1U1State (ComplexVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, double precision = MACHINE_PRECISION,  AbstractArchitecture* architecture = 0);
  
  // symmetrize a vector by grouping neighbouring orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // nbrOrbitals = number of orbitals to group together
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
  // kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
  // architecture = pointer to the architecture
  // return value = symmetrized state
  int SymmetrizeSingleStateGroupingNeighbouringOrbitals (RealVector& inputVector, int nbrOrbitals, RealVector*& symmetrizedVectors, int*& kySectors, AbstractArchitecture* architecture, double precision = MACHINE_PRECISION);

  // symmetrize a vector by grouping neighbouring orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // nbrOrbitals = number of orbitals to group together
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
  // kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
  // architecture = pointer to the architecture
  // return value = symmetrized state  
  virtual int SymmetrizeSingleStateGroupingNeighbouringOrbitals (ComplexVector& inputVector, int nbrOrbitals, ComplexVector*& symmetrizedVectors, 
								 int*& kySectors, AbstractArchitecture* architecture, double precision = MACHINE_PRECISION);

  // symmetrize a vector by grouping distant and equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // nbrOrbitals = number of orbitals to group together
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
  // kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
  // architecture = pointer to the architecture
  // return value = symmetrized state
  int SymmetrizeSingleStateGroupingDistantOrbitals (RealVector& inputVector, int nbrOrbitals, RealVector*& symmetrizedVectors, int*& kySectors, AbstractArchitecture* architecture, double precision = MACHINE_PRECISION);

  // symmetrize a vector by grouping distant and equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // nbrOrbitals = number of orbitals to group together
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
  // kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
  // architecture = pointer to the architecture
  // return value = symmetrized state
  int SymmetrizeSingleStateGroupingDistantOrbitals (ComplexVector& inputVector, int nbrOrbitals, ComplexVector*& symmetrizedVectors, int*& kySectors, AbstractArchitecture* architecture, double precision = MACHINE_PRECISION, bool twistedTorus = false);

  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // periodicity = momentum periodicity 
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest number of particles to the largest 
  //                      number of particles and the smallest Ky to the largest Ky
  // nbrParticlesSectors = reference on the array on the particle number sectors that have been generated through the symmetrization procedure
  // kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
  // architecture = pointer to the architecture
  // return value = number of states that have been generated through the symmetrization procedure
  virtual int SymmetrizeSingleStatePeriodicSubsetOrbitals (RealVector& inputVector, int firstOrbitalIndex, int periodicity, 
							   RealVector*& symmetrizedVectors, int*& nbrParticlesSectors, int*& kySectors, 
							   AbstractArchitecture* architecture);

  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // periodicity = momentum periodicity 
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest number of particles to the largest 
  //                      number of particles and the smallest Ky to the largest Ky
  // nbrParticlesSectors = reference on the array on the particle number sectors that have been generated through the symmetrization procedure
  // kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
  // architecture = pointer to the architecture
  // return value = number of states that have been generated through the symmetrization procedure
  virtual int SymmetrizeSingleStatePeriodicSubsetOrbitals (ComplexVector& inputVector, int firstOrbitalIndex, int periodicity, 
							   ComplexVector*& symmetrizedVectors, int*& nbrParticlesSectors, int*& kySectors, 
							   AbstractArchitecture* architecture);

  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // periodicity = momentum periodicity 
  // phase = an optional phase (in pi units) that can be added for each kept and occupied orbital
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest number of particles to the largest 
  //                      number of particles and the smallest Ky to the largest Ky
  // nbrParticlesSectors = reference on the array on the particle number sectors that have been generated through the symmetrization procedure
  // kySectors = reference on the array on twice the Ky sectors that have been generated through the symmetrization procedure
  // architecture = pointer to the architecture
  // return value = number of states that have been generated through the symmetrization procedure
  virtual int SymmetrizeSingleStatePeriodicSubsetOrbitals (ComplexVector& inputVector, int firstOrbitalIndex, int periodicity, double phase,
							   ComplexVector*& symmetrizedVectors, int*& nbrParticlesSectors, int*& kySectors, 
							   AbstractArchitecture* architecture);

  // create a state from its MPS description
  //
  // bMatrices = array that gives the B matrices 
  // twistMatrix = reference on the twist matrix to insert in the trace
  // state = reference to vector that will contain the state description
  // mPSSumIndices = diagonal indices that have to be kept in the trace
  // nbrMPSSumIndices = number of diagonal indices that have to be kept in the trace
  // memory = amount of memory that can be use to precompute matrix multiplications  
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, SparseRealMatrix& twistMatrix, RealVector& state, 
					      int* mPSSumIndices, int nbrMPSSumIndices,
					      long memory, long initialIndex, long nbrComponents);

  // create a state from its MPS description
  //
  // bMatrices = array that gives the B matrices 
  // twistMatrix = reference on the twist matrix to insert in the trace
  // state = reference to vector that will contain the state description
  // mPSSumIndices = diagonal indices that have to be kept in the trace
  // nbrMPSSumIndices = number of diagonal indices that have to be kept in the trace
  // memory = amount of memory that can be use to precompute matrix multiplications  
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (SparseComplexMatrix* bMatrices, SparseRealMatrix& twistMatrix, ComplexVector& state, 
					      int* mPSSumIndices, int nbrMPSSumIndices,
					      long memory, long initialIndex, long nbrComponents);

 protected:

  // core part of the C4 rotation
  //
  // inputState = reference on the state that has to be rotated
  // inputSpace = Hilbert space associated to the input state
  // outputState = reference on the rotated state
  // minIndex = minimum index that has to be computed
  // nbrIndices = number of indices that have to be computed
  // clockwise = the rotation is done clockwise
  // return value = reference on the rotated state
  virtual ComplexVector& CoreC4Rotation (ComplexVector& inputState, ParticleOnTorus* inputSpace, ComplexVector& outputState, int minIndex, int nbrIndices, bool clockwise);

  // symmetrized a product of two uncoupled states 
  //
  // outputVector = reference on the vector which will contain the symmetrozed state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // return value = symmetrized state
  virtual void SymmetrizeU1U1StateCore (RealVector& symmetrizedVector, RealVector& leftVector, RealVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents);

 // symmetrized a product of two uncoupled states 
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // return value = symmetrized state
  virtual void SymmetrizeU1U1StateCore (ComplexVector& symmetrizedVector, ComplexVector& leftVector, ComplexVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrized a product of several uncoupled states 
  //
  // outputState = reference on the output state
  // inputStates = states which will be symmetrized
  // inputSpaces = Hilbert spaces attached to each states
  // nbrStates = number of states to symmetrize
  // firstComponent = first component to symmetrize within the first Hilbert space of inputSpaces
  // nbrComponents = number of components to symmetrize within the first Hilbert space of inputSpaces
  virtual void SymmetrizeU1U1StateCore (RealVector& outputState, RealVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrized a product of several uncoupled states 
  //
  // outputState = reference on the output state
  // inputStates = states which will be symmetrized
  // inputSpaces = Hilbert spaces attached to each states
  // nbrStates = number of states to symmetrize
  // firstComponent = first component to symmetrize within the first Hilbert space of inputSpaces
  // nbrComponents = number of components to symmetrize within the first Hilbert space of inputSpaces
  virtual void SymmetrizeU1U1StateCore (ComplexVector& outputState, ComplexVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrize a vector by grouping neighbouring orbitals, core part
  //
  // inputVector = reference on the vector to symmetrize
  // nbrOrbitals = number of orbitals to group together
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
  // first component = index of the first vector component 
  // last component = index of the last component
  virtual void SymmetrizeSingleStateGroupingNeighbouringOrbitalsCore (ComplexVector& inputVector, ComplexVector* symmetrizedVectors, int nbrOrbitals, 
								unsigned long firstComponent, unsigned long nbrComponents);

  // symmetrize a vector by grouping distant and equally separated orbitals, core part
  //
  // inputVector = reference on the vector to symmetrize
  // nbrOrbitals = number of orbitals to group together
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
  // first component = index of the first vector component 
  // last component = index of the last component
  virtual void SymmetrizeSingleStateGroupingDistantOrbitalsCore (ComplexVector& inputVector, ComplexVector* symmetrizedVectors, int nbrOrbitals, unsigned long firstComponent, unsigned long nbrComponents, bool twistedTorus = false);

  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // symmetrizedVectors = array on the symmetrize states ranging from the smallest Ky to the largest Ky
  // periodicity = momentum periodicity (should be a multiple of the number of orbitals)
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  // return value = symmetrized state
  virtual void SymmetrizeSingleStatePeriodicSubsetOrbitalCore (ComplexVector& inputVector, ComplexVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, 
							       unsigned long firstComponent, unsigned long nbrComponents);

  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // symmetrizedVectors = array on the symmetrize states ranging from the smallest Ky to the largest Ky
  // periodicity = momentum periodicity (should be a multiple of the number of orbitals)
  // phase = an optional phase (in pi units) that can be added for each kept and occupied orbital
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  // return value = symmetrized state
  virtual void SymmetrizeSingleStatePeriodicSubsetOrbitalCore (ComplexVector& inputVector, ComplexVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, double phase, 
							       unsigned long firstComponent, unsigned long nbrComponents);


};

#endif


