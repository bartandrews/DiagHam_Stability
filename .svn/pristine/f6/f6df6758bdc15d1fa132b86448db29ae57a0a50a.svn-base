////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 10/10/2007                      //
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


#ifndef BOSONONSPHERESHORT_H
#define BOSONONSPHERESHORT_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "MathTools/LongRational.h"
#include "Vector/LongRationalVector.h"
#include "Matrix/LongRationalMatrix.h"

#include <iostream>
#include <map>

using std::map;

using std::cout;
using std::endl;
using std::dec;
using std::hex;

class FermionOnSphereWithSpin;
class BosonOnSphereFullShort;

class BosonOnSphereShort :  public ParticleOnSphere
{

  friend class FermionOnSphere;

  friend class BosonOnSphereHaldaneBasisShort;
  friend class BosonOnSphereSymmetricBasisShort;
  friend class BosonOnSphereHaldaneSymmetricBasisShort;
  friend class BosonOnSphereHaldaneHugeBasisShort;

  friend class BosonOnSphereFullShort;
  friend class FermionOnSphereTwoLandauLevels;
  friend class FermionOnSphereThreeLandauLevels;
  friend class FermionOnSphereFourLandauLevels;
	
  friend class FermionOnSphereWithSpin;
  friend class FermionOnSphereWithSpinHaldaneBasis;	
  friend class FermionOnSphereWithSpinHaldaneLzSzSymmetry;	

  friend class BosonOnTorusShort;

  friend class BosonOnSphereTwoLandauLevels;

  friend class BosonOnSphereWithSpin;
  friend class BosonOnSphereWithSU2Spin;
  friend class BosonOnSphereWithSU3Spin;
  friend class BosonOnSphereWithSU4Spin;
  friend class BosonOnSphereWithSU4SpinAllEntanglement;
  
  friend class BosonOnCP2TzSymmetry;

 protected:

  // the fermionic Hilbert space associated to the bosonic one
  FermionOnSphere* FermionBasis;

  // number of bosons
  int NbrBosons;
  // number of bosons plus 1
  int IncNbrBosons;
  // twice the momentum total value
  int TotalLz;
  // momentum total value shifted by LzMax / 2 * NbrBosons
  int ShiftedTotalLz;
  // twice the maximum Lz value reached by a boson
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;

  // pointer to the target space when an index is require after applying basic operation
  BosonOnSphereShort* TargetSpace;

  // pointer to an integer which indicate which coordinates are kept for the next time step iteration
  int* KeptCoordinates;
  // minors of permanents used for the time coherent wave function evaluation
  Complex** Minors;

  // temporary state used when applying operators
  unsigned long* TemporaryState;
  int TemporaryStateLzMax;
  // temporary state used when applying ProdA operator
  unsigned long* ProdATemporaryState;
  int ProdATemporaryStateLzMax;

 public:

  // default constructor
  //
  BosonOnSphereShort ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a boson
  BosonOnSphereShort (int nbrBosons, int totalLz, int lzMax);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereShort(const BosonOnSphereShort& bosons);

  // destructor
  //
  virtual ~BosonOnSphereShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereShort& operator = (const BosonOnSphereShort& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphere* targetSpace);

  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
  //
  // index = index of the state on which the operator has to be applied
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient);

  // apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int nbrIndices);

  // apply Prod_i a_mi operator to the state produced using ProdA method (without destroying it)
  // use double when calculating normalization factors to avoid overflow
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 

  virtual double ProdAL (int index, int* n, int nbrIndices);

  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int nbrIndices, double& coefficient);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  // use double when calculating normalization factors to avoid overflow
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 

  virtual int ProdAdL (int* m, int nbrIndices, double& coefficient);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (int index, int m);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (long index, int m);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient);
  
  // apply a_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A (int index, int n);
  
  // apply a^+_m operator to the state produced using the A or Ad method (without destroying it)
  //
  // m = first index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad (int m, double& coefficient);

  // apply a_n  operator to a given state. 
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value =  index of the resulting state 
  virtual int A (int index, int n, double& coefficient);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, long state);

  // print a given State using the monomial notation, with one column per particle (using space as a seperator)
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintColumnFormattedStateMonomial (ostream& Str, long state);

  // print a given state using the most compact notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintCompactState (ostream& Str, long state);

  // convert a fermionic state to its monomial representation
  //
  // index = index of the fermionic state
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void GetMonomial(long index, unsigned long*& finalState);

  // convert a state to its occupation number representation
  //
  // index = index of the state
  // finalState = reference on the array where the occupation number representation has to be stored
  virtual void GetOccupationNumber(long index, unsigned long*& finalState);

  // evaluate wave function in real space using a given basis and only for a given range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, int firstComponent, int nbrComponent);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = entanglement matrix of the subsytem
  virtual LongRationalMatrix EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int lzSector, LongRationalVector& groundState);
  
// reconstruct a state that contains only a certain subset of Schmidt eigenvalues of the given ground state
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// eigenvalueCut = only keep Schmidt levels that are larger than e^{-eigenvalueCut}
// groundState = reference on the total system ground state
// rebuiltSchmitGroundState = reference on the final state
// diagonalizedDensityMatrix = set of density matrix (Schmidt) eigenvalues
// transformationMatrix = the "truncation" matrix that connects the coefficients  (in the N-body basis) of the ground state and the final (truncated) state
// return value = reconstructed ground state vector
   virtual RealVector&  EvaluatePartialSchmidtDecomposition (int subsytemSize, int nbrBosonSector, int lzSector, double eigenvalueCut, RealVector& groundState, RealVector& rebuiltSchmidtGroundState, RealDiagonalMatrix& diagonalizedDensityMatrix, RealMatrix& transformationMatrix);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, ComplexVector& groundState);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particle. The geometrical cut is a stripe.
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax+shitedCut to -Lzmax+shitedCut+subsytemSize-1)
  // shiftedCut = first orbital belonging to the subsystem (with angular momentum -Lzmax+shitedCut)
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluateShiftedPartialDensityMatrix (int subsytemSize, int nbrShiftedOrbitals, int nbrBosonSector, int lzSector, RealVector& groundState);
	 
  // Compute the row and column dimension of the orbital entanglement matrix of 2 cuts counting only those rows/columns that are not completely zero
  // Also returns from the set of indices in the reduced density matrix corresponding to rows with atleast one non-zero entry
  // Columns contain the hilbert space of B and C which are traced out
  // SizeB = number of orbitals in part B, i.e. in the cap around Lzmax/2.
  // SizeA = number of orbitals in the bulk of the sphere 
  // NbrBosonsA = number of particles that belong to A
  // groundState = reference on the total system ground state
  // LzA = Lz sector of A in which the density matrix has to be evaluated as measured on a sphere with only A
  // return value = pointer with the 1st element being the row dimension
  // 2nd element is the column dimension of the oem 
  // 3rd element onward gives the positions of the rows in the oem/reduced density matrix which are not filled with 0's
  // (returns 0 if there is a probem/there is no hilbert space) 
  long* Compute2CutEntanglementMatrixDimensions (int SizeB, int SizeA, int NbrBosonsA, int LzA, RealVector& groundState);
  
  // evaluate a density matrix with 2 cuts of the whole system described by the RealVector groundState. The reduced density matrix is evaluated for a given Lz sector and number of particles
  //
  // SizeB = number of orbitals in part B, i.e. in the cap around Lzmax/2.
  // SizeA = number of orbitals in the bulk of the sphere 
  // NbrBosonsA = number of particles that belong to A
  // groundState = reference on the total system ground state
  // LzA = Lz sector of A in which the density matrix has to be evaluated as measured on a sphere with only A
  // return value = density matrix of the subsytem (return a zero dimension matrix if the density matrix is equal to zero)  
  RealSymmetricMatrix Evaluate2CutReducedDensityMatrix (int SizeB, int SizeA, int NbrBosonsA, int LzA, RealVector& groundState);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
								  RealVector& groundState, RealSymmetricMatrix* densityMatrix);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
								  ComplexVector& groundState, HermitianMatrix* densityMatrix);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using real space partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
  // thetaTop =  inclination angle defining one edge of the cut in degrees
  // thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  RealSymmetricMatrix EvaluatePartialDensityMatrixRealSpacePartition (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // core part of the evaluation density matrix real space partition calculation
  // 
  // minIndex = first index to consider in complementary Hilbert space
  // nbrIndex = number of indices to consider in complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // incompleteBetaThetaTop = pointer to the array where the top part coefficients are stored
  // incompleteBetaThetaBotton = pointer on the pointer to the array where the bottom part coefficients are stored
  // phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
  // return value = number of components that have been added to the density matrix
  long EvaluatePartialDensityMatrixRealSpacePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
							   RealVector& groundState,  RealSymmetricMatrix* densityMatrix, double* incompleteBetaThetaBottom, double* incompleteBetaThetaTop, double phiRange);

  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, bool removeBinomialCoefficient = false);

  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and both A and B are resticted to a given number of orbitals
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated
  // nbrOrbitalA = number of orbitals that have to be kept for the A part (starting from the leftmost orbital)
  // nbrOrbitalA = number of orbitals that have to be kept for the B part (starting from the rightmost orbital)
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int lzSector, int nbrOrbitalA, int nbrOrbitalB, 
									RealVector& groundState, bool removeBinomialCoefficient);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrOrbitalA = number of orbitals that have to be kept for the A part
  // weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
  // nbrOrbitalB = number of orbitals that have to be kept for the B part
  // weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixGenericRealSpacePartition (int nbrFermionSector, int lzSector, int nbrOrbitalA, double* weightOrbitalA, 
										     int nbrOrbitalB, double* weightOrbitalB, RealVector& groundState, 
										     AbstractArchitecture* architecture = 0);

  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition. The entanglement matrix is only evaluated in a given Lz sector.
  // and computed from precalculated particle entanglement matrix
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
  // thetaTop =  inclination angle defining one edge of the cut in degrees
  // thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
  // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
   virtual  RealMatrix& EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrBosonSector, int lzSector, 
												    double thetaTop, double thetaBottom, double phiRange, RealMatrix& entanglementMatrix);

  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrOrbitalA = number of orbitals that have to be kept for the A part
  // weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
  // nbrOrbitalB = number of orbitals that have to be kept for the B part
  // weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
  // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
   virtual RealMatrix& EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrBosonSector, int lzSector, 
													  int nbrOrbitalA, double* weightOrbitalA, 
													  int nbrOrbitalB, double* weightOrbitalB, RealMatrix& entanglementMatrix);

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has to be normalized to 1
  // symmetryFactor = if true also remove the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // convert a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // symmetryFactor = if true also add the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // convert a state such that its components are now expressed in the normalized basis, without applying the global normalization to the final state
  //
  // state = reference to the state to convert
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomialNoGlobalNormalization(RealVector& state);
 
  // convert a state such that its components, given in the conformal limit,  are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // return value = converted state
  virtual RealVector& ConvertFromConformalLimit(RealVector& state, long reference);

  // normalize Jack with respect to cylinder basis
  //
  // state = reference to the Jack state to normalize
  // aspect = aspect ratio of cylinder
  // return value = normalized state
  virtual RealVector& NormalizeJackToCylinder(RealVector& state, double aspect);

  // fuse two states which belong to different Hilbert spaces 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // leftVector = reference on the vector whose Hilbert space will be fuse to the left
  // rightVector = reference on the vector whose Hilbert space will be fuse to the right
  // padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
  // leftSpace = point to the Hilbert space that will be fuse to the left
  // rightSpace = point to the Hilbert space that will be fuse to the right
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // coefficient = optional multiplicative factor to apply to the fused state 
  // return value = reference on the fused state
  virtual RealVector& FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
				  ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace, bool symmetrizedFlag = false, double coefficient = 1.0);

  // fuse two states which belong to different Hilbert spaces 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // leftVector = reference on the vector whose Hilbert space will be fuse to the left
  // rightVector = reference on the vector whose Hilbert space will be fuse to the right
  // padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
  // leftSpace = point to the Hilbert space that will be fuse to the left
  // rightSpace = point to the Hilbert space that will be fuse to the right
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // coefficient = optional multiplicative factor to apply to the fused state 
  // return value = reference on the fused state
  virtual LongRationalVector& FuseStates (LongRationalVector& outputVector, LongRationalVector& leftVector, LongRationalVector& rightVector, int padding, 
					  ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace, bool symmetrizedFlag, LongRational& coefficient);

  // fuse multiple states which belong to different Hilbert spaces 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // nbrInputVectors = number of input vectors
  // inputVectors = input vectors whose Hilbert space will be fuse from  left to right
  // paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
  // inputSpaces = point to the Hilbert space that will be fuse to the left
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // coefficient = optional multiplicative factor to apply to the fused state 
  // return value = reference on the fused state
  virtual RealVector& FuseMultipleStates (RealVector& outputVector, int nbrInputVectors, RealVector* inputVectors, int* paddings, 
					  ParticleOnSphere** inputSpaces, bool symmetrizedFlag, double coefficient);

  // fuse multiple states which belong to different Hilbert spaces 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // nbrInputVectors = number of input vectors
  // inputVectors = input vectors whose Hilbert space will be fuse from  left to right
  // paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
  // inputSpaces = point to the Hilbert space that will be fuse to the left
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // coefficient = optional multiplicative factor to apply to the fused state 
  // return value = reference on the fused state
  virtual LongRationalVector& FuseMultipleStates (LongRationalVector& outputVector, int nbrInputVectors, LongRationalVector* inputVectors, int* paddings, 
						  ParticleOnSphere** inputSpaces, bool symmetrizedFlag, LongRational& coefficient);

  // use product rule to produce part of the components of a system from a smaller one
  //
  // outputVector = reference on the vector which will contain the product rule state  (without zeroing components which do not occur in the fusion)
  // inputVector = reference on the vector associated to the smaller system
  // inputSpace = pointer to the Hilbert space of the smaller system
  // commonPattern = array describing the shared leftmost pattern between the n-body states in both the smaller and larger system sizes
  // commonPatterSize = number of elements in the commonPattern array
  // addedPattern = array describing the pattern that has to be inserted to go from the smaller system to the larger one
  // addedPatterSize = number of elements in the addedPattern array
  // coefficient = multiplicqtive fqctor to go fron the component of the smaller system to the larger one
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  // return value = reference on the product rule state
  virtual RealVector& ProductRules (RealVector& outputVector, RealVector& inputVector, ParticleOnSphere* inputSpace, 
				    int* commonPattern, int commonPatterSize, int* addedPattern, int addedPatterSize,
				    double coefficient, bool symmetrizedFlag);

  // compute part of the Jack polynomial square normalization in a given range of indices
  //
  // state = reference on the unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  virtual double JackSqrNormalization (RealVector& outputVector, long minIndex = 0l, long nbrComponents = 0l);

  // compute part of the Jack polynomial square normalization in a given range of indices
  //
  // state = reference on the unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  virtual LongRational JackSqrNormalization (LongRationalVector& outputVector, long minIndex = 0l, long nbrComponents = 0l);

  // compute part of the Jack polynomial scalar product in a given range of indices
  //
  // state1 = reference on the first unnormalized Jack polynomial
  // state2 = reference on the second unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  virtual double JackScalarProduct (RealVector& state1, RealVector& state2, long minIndex = 0l, long nbrComponents = 0l);

  // compute part of the Jack polynomial square normalization in a given range of indices
  //
  // state1 = reference on the first unnormalized Jack polynomial
  // state2 = reference on the second unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  virtual LongRational JackScalarProduct (LongRationalVector& state1, LongRationalVector& state2, long minIndex = 0l, long nbrComponents = 0l);
  
  // symmetrized a product of two uncoupled states 
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // return value = symmetrized state
  virtual RealVector SymmetrizeU1U1State (RealVector& leftVector, RealVector& rightVector, BosonOnSphereShort* leftSpace, BosonOnSphereShort* rightSpace, bool unnormalizedBasisFlag = false, AbstractArchitecture* architecture = 0);
  

  // symmetrized a product of two uncoupled states 
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // return value = symmetrized state
  virtual void SymmetrizeU1U1StateCore (RealVector& symmetrizedVector, RealVector& leftVector, RealVector& rightVector, BosonOnSphereShort* leftSpace, BosonOnSphereShort* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents);

  // symmetrized a product of two uncoupled states, using rational input vectors
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // return value = symmetrized state
  virtual LongRationalVector SymmetrizeU1U1State (LongRationalVector& leftVector, LongRationalVector& rightVector, 
						  BosonOnSphereShort* leftSpace, BosonOnSphereShort* rightSpace, 
						  AbstractArchitecture* architecture = 0);
  

  // symmetrized a product of two uncoupled states, using rational input vectors
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // return value = symmetrized state
  virtual void SymmetrizeU1U1StateCore (LongRationalVector& symmetrizedVector, LongRationalVector& leftVector, LongRationalVector& rightVector, 
					BosonOnSphereShort* leftSpace, BosonOnSphereShort* rightSpace, 
					unsigned long firstComponent, unsigned long nbrComponents);

  // find state index from a string
  //
  // stateDescription = string describing the state
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(char* stateDescription);

  // find state index from an array of occupied orbitals
  //
  // stateDescription = array describing the state (stored as k1,k2,k3,...)
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(int* stateDescription);
  
  // find state index from unsigned long representation
  //
  // stateDescription = integer describint the state
  // lzMax = the lzmax of the state
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(unsigned long int stateDescription, int lzMax);

  // find state index from its occupation number description
  //
  // stateDescription = array that describes the state in the occupation number basis
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndexFromOccupationNumber(unsigned long* stateDescription);

  // get Lz component of a component
  //
  // j = index of the component in Hilbert space
  // return value = twice the Lz component
  virtual int GetLzValue(int j=0);

  // compute all Kostka coefficients for a given Schur polynomial 
  //
  // index = index of the partition that describe the Schur polynomial 
  // kostkaVector = vector where kostka numbers will be stored
  // return value = number of kostka coefficients
  virtual long KostkaForOneSchur(long index, RealVector& kostkaVector);
    
  // divide a set of fermionic vectors by a Jastrow factor to get their bosonic counterpart
  //
  // sourceVector = array of fermionic statesc be divided by a jastrow factor
  // outputVector = array of where bosonic states will be stored
  // firstComponent = index of the first component to transform 
  // nbrComponent = number of components to compute
  // nbrStates  = number of states to handle
  virtual void DivideByJastrow(RealVector * sourceVector, RealVector* OutputVector, long firstComponent, long nbrComponent, int nbrStates);

  virtual void FuseParticlesInState(RealVector& firstState, RealVector& outputVector, BosonOnSphereShort* finalSpace, long minIndex = 0l, long nbrComponents = 0l);
	
  // transform a vector belonging to this vector space in the lz->-lz
  //
  // finalSpace = the space obtained after the lz->-lz operation
  // initialVector = vector on which the operation will be apply
  // return value = vector resulting of the operation
  virtual RealVector GetLzSymmetricVector(ParticleOnSphere* finalSpace, RealVector& initialVector);

  // transform a vector belonging to this vector space in the lz->-lz
  //
  // finalSpace = the space obtained after the lz->-lz operation
  // initialVector = vector on which the operation will be apply
  // return value = vector resulting of the operation
  virtual LongRationalVector GetLzSymmetricVector(ParticleOnSphere* finalSpace, LongRationalVector& initialVector);
	
  // Compute the product of two states that belong to the same Hilbert Space
  //
  // firstState = reference on one of the states whose product will be computed
  // secondState = reference on the other state whose product will be computed
  // OutputVector = reference on the vector where the result will be stored
  // FinalSpace = pointer on the Hilbert Space whose the final state belong
  virtual void BosonicStateTimeBosonicState(RealVector& firstState, RealVector& secondState, RealVector& outputVector, int minIndex, int nbrComponents, BosonOnSphereShort * finalSpace);

  // Compute the product of two states that belong to different Hilbert Spaces
  //
  // firstState = reference on one of the states whose product will be computed
  // secondState = reference on the other state whose product will be computed
  // OutputVector = reference on the vector where the result will be stored
  // SecondSpace = pointer on the Hilbert Space whose the second state belong
  // minIndex = first computed component
  // nbrComponents = Nomber of computed components
  // FinalSpace = pointer on the Hilbert Space whose the final state belong
  virtual void BosonicStateTimeBosonicState(RealVector& firstState, RealVector& secondState, RealVector& outputVector, BosonOnSphereShort * secondSpace, int minIndex, int nbrComponents, BosonOnSphereShort * finalSpace);
  
  // Compute the product of two Monomials
  //
  // firstState = array where the monomial representation of the first state is stored
  // secondState = array where the monomial representation of the second state is stored
  // finalStates = reference on the array where the fermionic representation of the states product will be stored
  // weigth = reference on the array where the coefficient of the states product will be stored
  // FinalSpace = pointer on the Hilbert Space whose the final monomials belong
  // return value = number of different states generated by the product
  virtual unsigned long ProductOfTwoMonomials (unsigned long* firstState,int * equalPowerIndex,const int nbrEqualPower,unsigned long* secondState, unsigned long * & finalStates, long * & weigth,BosonOnSphereShort * finalSpace);
  
  // Compute the coefficient of a monomial in the decomposition of a product of two monomials
  //
  // state = array where the monomial representation of the state whose coefficient is to be computed is stored
  // firstState = array where the monomial representation of one of the states whose product is computed, is stored
  // return value = coefficient of the monomial stored in state
  virtual long ComputeCoefficient(unsigned long * state,const unsigned long * firstState);
  
  // Compute the product of a bosonic state and a fermionic state
  //
  // firstState = reference on the bosonic state
  // secondState = reference on the fermionic state
  // outputVector = reference on the vector where the result will be stored
  // fermionSpace = pointer on the fermionic Hilbert Space whose the second state belong
  // minIndex = first computed component
  // nbrComponents = Nomber of computed components
  // finalSpace = pointer on the Hilbert Space whose the final state belong to
  virtual void BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, FermionOnSphere * fermionSpace, int minIndex, int nbrComponents, FermionOnSphere * finalSpace);
  
  // Compute the product of a Monomial and a Slater determinant
  //
  // firstState = array where the monomial representation of the first state is stored
  // secondState = array where the monomial representation of the second state is stored
  // finalStates = reference on the array where the fermionic representation of the states product will be stored
  // weigth = reference on the array where the coefficient of the states product will be stored
  // finalSpace = pointer to the destination Hilbert space
  // return value = number of different states generated by the product
  virtual unsigned long MonomialsTimesSlater (unsigned long* slater, unsigned long* monomial, unsigned long * & finalStates, long * & weigth, FermionOnSphere * finalSpace);
  
  // Compute the product of two fermionic states that belong to different Hilbert Spaces
  //
  // firstState = reference on one of the states whose product will be computed
  // secondState = reference on the other state whose product will be computed
  // OutputVector = reference on the vector where the result will be stored
  // fermionSpace1 = pointer on the Hilbert Space whose the first state belong to
  // fermionSpace2 = pointer on the Hilbert Space whose the second state belong to 
  // minIndex = first computed component
  // nbrComponents = Nomber of computed components
  virtual void FermionicStateTimeFermionicState(RealVector& fermionState1, RealVector& fermionState2, RealVector& outputVector, FermionOnSphere * fermionSpace1, FermionOnSphere * fermionSpace2, int minIndex, int nbrComponents);
  
  // Compute the product of two Slater determinants
  //
  // slater1 = array where the monomial representation of the first slater is stored
  // slater2 = array where the monomial representation of the second slater is stored
  // finalStates = reference on the array where the fermionic representation of the states product will be stored
  // weigth = reference on the array where the coefficient of the states product will be stored
  virtual unsigned long SlaterTimesSlater (unsigned long* slater1,unsigned long* slater2, unsigned long * & finalStates, long * & weigth);
  
  // Compute the product of two fermionic states that belong to the same Hilbert Spaces
  //
  // fermionState1 = reference on one of the states whose product will be computed
  // fermionState2 = reference on the other state whose product will be computed
  // OutputVector = reference on the vector where the result will be stored
  // fermionSpace = pointer on the fermionic Hilbert Space
  // minIndex = first computed component
  // nbrComponents = Nomber of computed components
  virtual void FermionicStateTimeFermionicState(RealVector& fermionState1, RealVector& fermionState2, RealVector& outputVector, FermionOnSphere * fermionSpace, int minIndex, int nbrComponents);

  // compute the projection of the product of a bosonic state and the halperin 110 state
  //
  // bosonState = real vector where the bosonic state is stored
  // outputVector = real vector where the result has to be stored
  // fermionSpace = pointer to the fermionic Hilbert space
  // finalSpace = pointer to the final Hilbert space
  // firstComponent = first component to be computed
  // nbrComponent = number of components to be computed	
  virtual void BosonicStateTimePolarizedSlaters(RealVector& bosonState, RealVector& outputVector, FermionOnSphere * fermionSpace , FermionOnSphereWithSpin* finalSpace, int firstComponent,int nbrComponent);
	
  
  virtual void EvaluatePartialDensityMatrixMultipartiteParticlePartition(ParticleOnSphere * spaceA, ParticleOnSphere * spaceB, ParticleOnSphere * spaceC,  RealVector groundstate, RealSymmetricMatrix* densityMatrix, AbstractArchitecture* architecture = 0);
 
  // symmetrize a vector with even number of orbitals 
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector to be symmetrized
  // leftSpace = pointer to the Hilbert space
  // return value = symmetrized state
  virtual RealVector SymmetrizeU1U1SingleState (RealVector& leftVector, BosonOnSphereShort* leftSpace, bool unnormalizedBasisFlag, AbstractArchitecture* architecture = 0);
  
  // symmetrize a vector with even number of orbitals and rational coefficients
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector to be symmetrized
  // leftSpace = pointer to the Hilbert space
  // return value = symmetrized state
  virtual LongRationalVector SymmetrizeU1U1SingleState (LongRationalVector& leftVector, BosonOnSphereShort* leftSpace, bool unnormalizedBasisFlag, AbstractArchitecture* architecture = 0);
  
  // symmetrize a vector with even number of orbitals
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // leftSpace = pointer to the Hilbert space of the first color
  // return value = symmetrized state
  virtual void SymmetrizeU1U1SingleStateOneInTwoCore (RealVector& symmetrizedVector, RealVector& leftVector, BosonOnSphereShort* leftSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrize a vector with even number of orbitals and rational components
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // leftSpace = pointer to the Hilbert space of the first color
  // return value = symmetrized state
  virtual void SymmetrizeU1U1SingleStateOneInTwoCore (LongRationalVector& symmetrizedVector, LongRationalVector& leftVector, BosonOnSphereShort* leftSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrize a vector by grouping several orbitals into a single one
  //
  // inputVector = reference on the vector to symmetrize
  // nbrOrbitals = number of orbitals to group together
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Lz to the largest Lz
  // lzSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
  // return value = number of states that have been generated through the symmetrization procedure
  virtual int SymmetrizeSingleStateOneIntoManyOrbital (LongRationalVector& inputVector, int nbrOrbitals, LongRationalVector*& symmetrizedVectors, int*& lzSectors);
  
  // symmetrize a vector by grouping several orbitals that are related by a periodicity condition on their momenta
  //
  // inputVector = reference on the vector to symmetrize
  // periodicity = momentum periodicity (should be a multiple of the number of orbitals)
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Lz to the largest Lz
  // lzSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
  // return value = number of states that have been generated through the symmetrization procedure
  virtual int SymmetrizeSingleStatePeriodicOrbitals (LongRationalVector& inputVector, int periodicity, LongRationalVector*& symmetrizedVectors, int*& lzSectors);

  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // periodicity = momentum periodicity 
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest number of particles to the largest 
  //                      number of particles and the smallest Lz to the largest Lz
  // nbrParticlesSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
  // lzSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
  // return value = number of states that have been generated through the symmetrization procedure
  virtual int SymmetrizeSingleStatePeriodicSubsetOrbitals (LongRationalVector& inputVector, int firstOrbitalIndex, int periodicity, 
							   LongRationalVector*& symmetrizedVectors, int*& nbrParticlesSectors, int*& lzSectors);
  
  // symmetrize a vector by grouping several orbitals into a single one
  //
  // inputVector = reference on the vector to symmetrize
  // nbrOrbitals = number of orbitals to group together
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Lz to the largest Lz
  // lzSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
  // return value = number of states that have been generated through the symmetrization procedure
  virtual int SymmetrizeSingleStateOneIntoManyOrbital (RealVector& inputVector, int nbrOrbitals, RealVector*& symmetrizedVectors, bool unnormalizedBasisFlag, int*& lzSectors);
  
  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // periodicity = momentum periodicity 
  // symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest number of particles to the largest 
  //                      number of particles and the smallest Lz to the largest Lz
  // nbrParticlesSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
  // lzSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
  // return value = number of states that have been generated through the symmetrization procedure
  virtual int SymmetrizeSingleStatePeriodicSubsetOrbitals (RealVector& inputVector, int firstOrbitalIndex, int periodicity, bool unnormalizedBasisFlag,
							   RealVector*& symmetrizedVectors, int*& nbrParticlesSectors, int*& lzSectors);

  // Compute the product of two bsonic states, automatically dealing with reverse flux attachement
  //
  // bosonicState1 = reference on the first bosonic state
  // bosonicState2 = reference on the second fermionic state
  // outputVector = reference on the vector where the result will be stored
  // bosonicSpace1 = pointer on the Hilbert Space associated to the first bosonic state
  // bosonicSpace2 = pointer on the Hilbert Space associated to the second bosonic state
  // minIndex = first component to compute (refering to the bosonic state)
  // nbrComponents = number of components to compute (refering to the first bosonic state)
  // unnormalizedFlag = true if the state should be written in the unnormalized basis
  // architecture = pointer to the architecture
  virtual void BosonicStateTimeBosonicState(RealVector& bosonicState1, RealVector& bosonicState2, RealVector& outputVector, 
					    BosonOnSphereShort* bosonicSpace1, BosonOnSphereShort* bosonicSpace2,
					    int minIndex, int nbrComponents, bool unnormalizedFlag, AbstractArchitecture* architecture);

 protected:

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initialbosonic  state is stored
  // initialStateLzMax = reference on the initial bosonic state maximum Lz value
  // return value = corresponding fermionic state
  unsigned long BosonToFermion(unsigned long*& initialState, int& initialStateLzMax);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax = initial fermionic state maximum Lz value
  // finalState = reference on the array where the bosonic state has to be stored
  // finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored
  void FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax);

  // convert a fermionic state into its bosonic counterpart
  //
  // initialState = initial fermionic state
  // finalState = reference on the array where the bosonic state has to be stored
  void FermionToBoson(unsigned long initialState, unsigned long*& finalState);

  // check if a state satisfies a maximum occupation contraint
  //
  // initialState =  initial fermionic state
  // maximumOccupation = orbital maximum occupation
  // return value = true if the maximum occupation contraint is satisfied
  bool CheckMaximumOccupation(unsigned long initialState, unsigned int maximumOccupation);

  // convert a bosonic state to its monomial representation
  //
  // initialState = initial  bosonic state
  // initialStateLzMax = initial bosonic state maximum Lz value
  // finalState = reference on the array where the monomial representation has to be stored
  void ConvertToMonomial(unsigned long* initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a bosonic state to its monomial representation
  //
  // initialState = initial bosonic state in its fermionic representation
  // initialStateLzMax = initial bosonic state maximum Lz value
  // finalState = reference on the array where the monomial representation has to be stored
  void ConvertToMonomial(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a bosonic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // return value = bosonic state in its fermionic representation
  unsigned long ConvertFromMonomial(unsigned long* initialState);

  // convert a bosonic state from its monomial representation, assuming the array is not sorted
  //
  // initialState = array where the monomial representation is stored
  // return value = bosonic state in its fermionic representation
  unsigned long ConvertFromUnsortedMonomial(unsigned long* initialState);
 
  // convert a bosonic state from its monomial representation, assuming the array is not sorted and find the index of the last occupied orbital
  //
  // initialState = array where the monomial representation is stored
  // maximumLzMax = reference on the index of the last occupied orbital
  // return value = bosonic state in its fermionic representation
  unsigned long ConvertFromUnsortedMonomial(unsigned long* initialState, int& maximumLzMax);

  // fuse particles two by two in a given monomial
  //
  // index = monomial index
  // finalSpace = space where the fused state lies
  // weigthVector = weigths of each fused component
  // indicesVector = indices of each fused component
  // return value = number of generated components when fusing
  virtual int FuseParticlesInMonomial(long index, BosonOnSphereShort* finalSpace, long* weigthVector, unsigned long* indicesVector);
	
  virtual int GeneratePairs(unsigned long* Monomial, long* weigthVector, unsigned long* indicesVector, unsigned long* FinalMonomial, int reste, int compteur, BosonOnSphereShort * finalSpace);

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual bool HasPauliExclusions(int index, int pauliK, int pauliR);

  // core part of multiple state fuse 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // nbrInputVectors = number of input vectors
  // inputVectors = input vectors whose Hilbert space will be fuse from  left to right
  // paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
  // inputSpaces = point to the Hilbert space that will be fuse to the left
  // currentPosition = index of the current space to fuse
  // currentState = current fermionic state obtained by fusing previous states
  // currentCoefficient = current multiplicative coefficient
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  virtual void CoreFuseMultipleStates (RealVector& outputVector, int nbrInputVectors, RealVector* inputVectors, int* paddings, BosonOnSphereShort** inputSpaces, int currentPosition, unsigned long currentState, int currentPadding, double currentCoefficient, bool symmetrizedFlag);

  // core part of multiple state fuse 
  //
  // outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
  // nbrInputVectors = number of input vectors
  // inputVectors = input vectors whose Hilbert space will be fuse from  left to right
  // paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
  // inputSpaces = point to the Hilbert space that will be fuse to the left
  // currentPosition = index of the current space to fuse
  // currentState = current fermionic state obtained by fusing previous states
  // currentCoefficient = current multiplicative coefficient
  // symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
  virtual void CoreFuseMultipleStates (LongRationalVector& outputVector, int nbrInputVectors, LongRationalVector* inputVectors, int* paddings, BosonOnSphereShort** inputSpaces, int currentPosition, unsigned long currentState, int currentPadding, LongRational currentCoefficient, bool symmetrizedFlag);
	
	
  // check that the product firstState*secondState is in the lexicographical order
  //
  // firstState = array where the monomial representation  of a state is stored
  // secondState = array where the monomial representation  of another state is stored
  // return value = true if the product firstState*secondState is in the lexicographical order
  virtual bool CheckLexiOrder(int * firstState,unsigned long* secondState,int TailleEgal);

  // symmetrize a vector by grouping several orbitals into a single one
  //
  // inputVector = reference on the vector to symmetrize
  // symmetrizedVectors = array on the symmetrized states ranging from the smallest Lz to the largest Lz
  // nbrOrbitals = number of orbitals to group together
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  virtual void SymmetrizeSingleStateOneIntoManyOrbitalCore (LongRationalVector& inputVector, LongRationalVector* symmetrizedVectors, int nbrOrbitals, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrize a vector by grouping several orbitals that are related by a periodicity condition on their momenta
  //
  // inputVector = reference on the vector to symmetrize
  // symmetrizedVectors = array on the symmetrize states ranging from the smallest Lz to the largest Lz
  // periodicity = momentum periodicity (should be a multiple of the number of orbitals)
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  // return value = symmetrized state
  virtual void SymmetrizeSingleStatePeriodicOrbitalCore (LongRationalVector& inputVector, LongRationalVector* symmetrizedVectors, int periodicity, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // symmetrizedVectors = array on the symmetrize states ranging from the smallest Lz to the largest Lz
  // periodicity = momentum periodicity (should be a multiple of the number of orbitals)
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  // return value = symmetrized state
  virtual void SymmetrizeSingleStatePeriodicSubsetOrbitalCore (LongRationalVector& inputVector, LongRationalVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, 
							       unsigned long firstComponent, unsigned long nbrComponents);

  // symmetrize a vector by grouping several orbitals into a single one
  //
  // inputVector = reference on the vector to symmetrize
  // symmetrizedVectors = array on the symmetrized states ranging from the smallest Lz to the largest Lz
  // nbrOrbitals = number of orbitals to group together
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  virtual void SymmetrizeSingleStateOneIntoManyOrbitalCore (RealVector& inputVector, RealVector* symmetrizedVectors, bool unnormalizedBasisFlag, int nbrOrbitals, unsigned long firstComponent, unsigned long nbrComponents);
  
  // symmetrize a vector by keeping only a subset of equally separated orbitals
  //
  // inputVector = reference on the vector to symmetrize
  // firstOrbitalIndex = index of the first orbital to keep
  // symmetrizedVectors = array on the symmetrize states ranging from the smallest Lz to the largest Lz
  // periodicity = momentum periodicity (should be a multiple of the number of orbitals)
  // firstComponent = first component of the input vector that has to be symmetrized
  // nbrComponents = number of components of the input vector that have to be symmetrized
  // return value = symmetrized state
  virtual void SymmetrizeSingleStatePeriodicSubsetOrbitalCore (RealVector& inputVector, RealVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, bool unnormalizedBasisFlag,
							       unsigned long firstComponent, unsigned long nbrComponents);

  // Compute the product of two symmetric monomials
  //
  // symmetricMonomial1 = first symmetric monomial
  // symmetricMonomial2 = second symmetric monomial
  // finalState = reference on the vector the produced state will be stored
  // threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
  virtual void SymmetricMonomialTimesSymmetricMonomial (unsigned long* symmetricMonomial1, unsigned long* symmetricMonomial2, RealVector& finalState, double** threeOrbitalOverlaps);

  // Compute the product of two symmetric monomials, assuming an unnormalized basis
  //
  // symmetricMonomial1 = first symmetric monomial
  // symmetricMonomial2 = second symmetric monomial
  // finalState = reference on the vector the produced state will be stored
  virtual void UnnormalizedSymmetricMonomialTimesSymmetricMonomial (unsigned long* symmetricMonomial1, unsigned long* symmetricMonomial2, RealVector& finalState);
  
  // Compute the product of two symmetric monomials, assuming an unnormalized basis
  //
  // symmetricMonomial1 = first symmetric monomial
  // symmetricMonomial2 = second symmetric monomial
  // finalStateConfigurations = fermionic configurations that are generated durig the multiplication
  // finalStateWeights = weights associated to each generated fermionic configurations
  // return value = number of generated fermionic configurations
  virtual int UnnormalizedSymmetricMonomialTimesSymmetricMonomial (unsigned long* symmetricMonomial1, unsigned long* symmetricMonomial2, 
								   unsigned long* finalStateConfigurations, double* finalStateWeights);

  // Compute the product of two symmetric monomials, assuming a reverse flux attachment for the first symmetric monomial
  //
  // symmetricMonomial1 = first symmetric monomial
  // symmetricMonomial2 = second symmetric monomial
  // finalState = reference on the vector the produced state will be stored
  // threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
  virtual void ReverseSymmetricMonomialTimesSymmetricMonomial (unsigned long* symmetricMonomial1, unsigned long* symmetricMonomial2, RealVector& finalState, double** threeOrbitalOverlaps);


};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnSphereShort::GetParticleStatistic()
{
  return ParticleOnSphere::BosonicStatistic;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initialbosonic  state is stored
// initialStateLzMax = reference on the initial bosonic state maximum Lz value
// return value = corresponding fermionic state

inline unsigned long BosonOnSphereShort::BosonToFermion(unsigned long*& initialState, int& initialStateLzMax)
{
  unsigned long TmpState = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= initialStateLzMax; ++i)
    {
      TmpState |= ((1ul << initialState[i]) - 1ul) << Shift;
      Shift += initialState[i];
      ++Shift;
    }
  return TmpState;
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateLzMax = initial fermionic state maximum Lz value
// finalState = reference on the array where the bosonic state has to be stored
// finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored

inline void BosonOnSphereShort::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax)
{
  finalStateLzMax = 0;
  while (initialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialState - 1ul) ^ (~initialState);
      TmpState &= ~(TmpState >> 1);
//      cout << hex << initialState << "  " << TmpState << dec << endl;
#ifdef  __64_BITS__
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
      TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
//      cout << TmpPower << endl;
      finalState[finalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState >>= TmpPower;
      ++finalStateLzMax;
      initialStateLzMax -= TmpPower;
    }
  --finalStateLzMax;
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// finalState = reference on the array where the bosonic state has to be stored

inline void BosonOnSphereShort::FermionToBoson(unsigned long initialState, unsigned long*& finalState)
{
  int FinalStateLzMax = 0;
  int InitialStateLzMax = this->FermionBasis->LzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialState - 1ul) ^ (~initialState);
      TmpState &= ~(TmpState >> 1);
//      cout << hex << initialState << "  " << TmpState << dec << endl;
#ifdef  __64_BITS__
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
      TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
//      cout << TmpPower << endl;
      finalState[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      ++FinalStateLzMax;
      initialState >>= TmpPower;
      InitialStateLzMax -= TmpPower;
    }
  while (FinalStateLzMax <= this->LzMax)
    {
      finalState[FinalStateLzMax] = 0x0ul;
      ++FinalStateLzMax;
    }
}

// check if a state satisfies a maximum occupation contraint
//
// initialState =  initial fermionic state
// maximumOccupation = orbital maximum occupation
// return value = true if the maximum occupation contraint is satisfied

inline bool BosonOnSphereShort::CheckMaximumOccupation(unsigned long initialState, unsigned int maximumOccupation)
{
  while (initialState != 0x0ul)
    {
      unsigned long TmpState = (~initialState - 1ul) ^ (~initialState);
      TmpState &= ~(TmpState >> 1);
#ifdef  __64_BITS__
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
      TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
      if (TmpPower > maximumOccupation)
	return false;
      ++TmpPower;
      initialState >>= TmpPower;
    }
  return true;
}

// convert a fermionic state to its monomial representation
//
// index = index of the fermionic state
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnSphereShort::GetMonomial(long index, unsigned long*& finalState)
{
  int Index = 0;
  unsigned long InitialState = this->FermionBasis->StateDescription[index];
  int InitialStateLzMax  = this->FermionBasis->LzMax;
  int TmpLz = InitialStateLzMax  - this->NbrBosons + 1;
  while (InitialStateLzMax >= 0)
    {
      while ((InitialStateLzMax >= 0) && (((InitialState >> InitialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalState[Index++] = (unsigned long) TmpLz;
	  --InitialStateLzMax;
	}
      while ((InitialStateLzMax >= 0) && (((InitialState >> InitialStateLzMax) & 0x1ul) == 0x0ul))
	{
	  --TmpLz;
	  --InitialStateLzMax;
	}
    }
}

// convert a state to its occupation number representation
//
// index = index of the state
// finalState = reference on the array where the occupation number representation has to be stored

inline void BosonOnSphereShort::GetOccupationNumber(long index, unsigned long*& finalState)
{
  int TmpLzMax = 0;
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], finalState, TmpLzMax);
  for (TmpLzMax += 1; TmpLzMax <= this->LzMax; ++TmpLzMax) 
    finalState[TmpLzMax] = 0x0ul; 
}



// convert a bosonic state to its monomial representation
//
// initialState = initial  bosonic state
// initialStateLzMax = initial bosonic state maximum Lz value
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnSphereShort::ConvertToMonomial(unsigned long* initialState, int initialStateLzMax, unsigned long*& finalState)
{
  int Index = 0;
  for (int i = initialStateLzMax; i >= 0; --i)
    for (unsigned long j = 0l; j < initialState[i]; ++j)
      finalState[Index++] = i;
}

// convert a bosonic state to its monomial representation
//
// initialState = initial bosonic state in its fermionic representation
// initialStateLzMax = initial bosonic state maximum Lz value
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnSphereShort::ConvertToMonomial(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState)
{
  int Index = 0;
  int TmpLz = initialStateLzMax - this->NbrBosons + 1;
  while (initialStateLzMax >= 0)
    {
      while ((initialStateLzMax >= 0) && (((initialState >> initialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalState[Index++] = TmpLz;
	  --initialStateLzMax;
	}
      while ((initialStateLzMax >= 0) && (((initialState >> initialStateLzMax) & 0x1ul) == 0x0ul))
	{
	  --TmpLz;
	  --initialStateLzMax;
	}
    }
}

// convert a bosonic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// return value = bosonic state in its fermionic representation

inline unsigned long BosonOnSphereShort::ConvertFromMonomial(unsigned long* initialState)
{
  unsigned long Tmp = 0x0ul;
  for (int i = 0; i < this->NbrBosons; ++i)
    Tmp |= 0x1ul << (initialState[i] + ((unsigned long) (this->NbrBosons - i)) - 1ul);
  return Tmp;
}

// convert a bosonic state from its monomial representation, assuming the array is not sorted
//
// initialState = array where the monomial representation is stored
// return value = bosonic state in its fermionic representation

inline unsigned long BosonOnSphereShort::ConvertFromUnsortedMonomial(unsigned long* initialState)
{
  for (int i = 0; i <= this->LzMax; ++i)
    this->TemporaryState[i] = 0x0ul;
  for (int i = 0; i < this->NbrBosons; ++i)
    this->TemporaryState[initialState[i]]++;
  unsigned long TmpState = 0x0ul;
  unsigned long Shift = 0;
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0x0ul)
    {
      --this->TemporaryStateLzMax;
    }
  for (int i = 0; i <= this->TemporaryStateLzMax; ++i)
     {
      TmpState |= ((1ul << this->TemporaryState[i]) - 1ul) << Shift;
      Shift += this->TemporaryState[i];
      ++Shift;
    }
  return TmpState;
}

// convert a bosonic state from its monomial representation, assuming the array is not sorted and find the index of the last occupied orbital
//
// initialState = array where the monomial representation is stored
// maximumLzMax = reference on the index of the last occupied orbital
// return value = bosonic state in its fermionic representation

inline unsigned long BosonOnSphereShort::ConvertFromUnsortedMonomial(unsigned long* initialState, int& maximumLzMax)
{
  for (int i = 0; i <= this->LzMax; ++i)
    this->TemporaryState[i] = 0x0ul;
  for (int i = 0; i < this->NbrBosons; ++i)
    this->TemporaryState[initialState[i]]++;
  unsigned long TmpState = 0x0ul;
  unsigned long Shift = 0x0ul;
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0x0ul)
    {
      --this->TemporaryStateLzMax;
    }
  maximumLzMax = this->TemporaryStateLzMax;
  maximumLzMax += this->NbrBosons;
  --maximumLzMax;
  for (int i = 0; i <= this->TemporaryStateLzMax; ++i)
     {
      TmpState |= ((1ul << this->TemporaryState[i]) - 1ul) << Shift;
      Shift += this->TemporaryState[i];
      ++Shift;
    }
  return TmpState;
}


// check that the product firstState*secondState is in the lexicographical order
//
// firstState = array where the monomial representation  of a state is stored
// secondState = array where the monomial representation  of another state is stored
// return value = true if the product firstState*secondState is in the lexicographical order

inline bool BosonOnSphereShort::CheckLexiOrder(int * egal,unsigned long* secondState,int TailleEgal)
{
  for (int index = 0; index < TailleEgal; ++index)
    {
      if (secondState[egal[index]]<secondState[egal[index]+1])
	{
	  return false;
	}
    }
  return true;
}

// print a given state using the most compact notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

inline ostream& BosonOnSphereShort::PrintCompactState (ostream& Str, long state)
{
  return this->FermionBasis->PrintCompactState(Str, state);
}

#endif


