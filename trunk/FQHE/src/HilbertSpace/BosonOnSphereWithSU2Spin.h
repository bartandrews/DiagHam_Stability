////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//      altenate version of the class for bosons on sphere with SU(2) spin    //
//                                                                            //
//                        last modification : 26/01/2012                      //
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


#ifndef BOSONONSPHEREWITHSU2SPIN_H
#define BOSONONSPHEREWITHSU2SPIN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"

#include <iostream>

using std::cout;
using std::endl;


class BosonOnSphereShort;
class FermionOnSphereWithSpin;
class FermionOnSphereWithSpinTwoLandauLevels;


class BosonOnSphereWithSU2Spin :  public ParticleOnSphereWithSpin
{


  friend class BosonOnSphereWithSU4Spin;
  friend class BosonOnSphereWithSU2SpinSzSymmetry;
  friend class BosonOnSphereWithSU2SpinLzSymmetry;
  friend class BosonOnSphereWithSU2SpinLzSzSymmetry;
  friend class FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation;

 protected:

  // number of bosons
  int NbrBosons;
  // number of bosons plus 1
  int IncNbrBosons;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a boson
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;
  // twice the total spin value
  int TotalSpin;
  // number of bosons with spin up
  int NbrBosonsUp;
  // number of bosons with spin down
  int NbrBosonsDown;

  // lzmax value for the fermionic states associated to the type up particles
  int NUpLzMax;
  // lzmax value for the fermionic states associated to the type down particles
  int NDownLzMax;
  // largest lzmax value for the fermionic states among N1LzMax, N2LzMax and N3LzMax
  int FermionicLzMax;

  // array that contains the state description, the first entry being StateDescriptionUp and the second entry being StateDescriptionDown
  unsigned long* StateDescriptionSigma[2];

  // temporay array describing the type up particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionUp;
  // temporay array describing the type down particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionDown;

  // sorted array that contains each unique configuration for the type up particles
  unsigned long* UniqueStateDescriptionUp;
  // number of time each unique configuration for the type up particles appears in StateDescriptionUp
  int* UniqueStateDescriptionSubArraySizeUp;
  // number of unique configurations for the type up-plus particles
  long NbrUniqueStateDescriptionUp;
  // first time a type up appears in the Hilbert space
  int* FirstIndexUniqueStateDescriptionUp;

  // temporary state used when applying operators for type up particles
  unsigned long* TemporaryStateUp;
  // temporary state used when applying operators for type down particles
  unsigned long* TemporaryStateDown;
  // array that contains the temporary state, the first entry being TemporaryStateUp and the second entry being TemporaryStateDown
  unsigned long* TemporaryStateSigma[2];

  // temporary state used when applying ProdA operator for type up particles
  unsigned long* ProdATemporaryStateUp;
  // temporary state used when applying ProdA operator for type down particles
  unsigned long* ProdATemporaryStateDown;
  // array that contains the temporary state used when applying ProdA operator, the first entry being ProdATemporaryStateUp and the second entry being ProdATemporaryStateDown
  unsigned long* ProdATemporaryStateSigma[2];

  // pointer to the Hilbert space where the result of any operator lies
  BosonOnSphereWithSU2Spin* TargetSpace;

 public:

  // default constructor
  // 
  BosonOnSphereWithSU2Spin ();

  // basic constructor without any constraint on Sz
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  BosonOnSphereWithSU2Spin (int nbrBosons, int totalLz, int lzMax);

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // totalSpin = twice the total spin value
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2Spin (int nbrBosons, int totalLz, int lzMax, int totalSpin, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  // 
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2Spin (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSU2Spin(const BosonOnSphereWithSU2Spin& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSU2Spin ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSU2Spin& operator = (const BosonOnSphereWithSU2Spin& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // get the number of orbitals
  //
  // return value = number of orbitals
  virtual int GetNbrOrbitals();

  // get the number of particles
  //
  // return value = number of particles
  virtual int GetNbrParticles();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphereWithSpin* targetSpace);

  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

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

  // apply a^+_m_d a_m_d operator to a given state (only spin down isospin plus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_dm a_m_dm
  virtual double AddAd (int index, int m);

  // apply a^+_m_u a_m_u operator to a given state  (only spin up isospin plus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_um a_m_um
  virtual double AduAu (int index, int m);

  // apply a^+_m_u a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAu (int index, int m, int n, double& coefficient);

  // apply a^+_m_u a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAd (int index, int m, int n, double& coefficient);

  // apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // sigma1 = SU(2) index for the first annihilation operator
  // sigma2 = SU(2) index for the second annihilation operator
  // return value =  multiplicative factor 
  virtual double AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AuAu (int index, int n1, int n2);

  // apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AuAd (int index, int n1, int n2);

  // apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AdAd (int index, int n1, int n2);

  // apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // sigma1 = SU(2) index for the first creation operator
  // sigma2 = SU(2) index for the second creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient);

  // apply a^+_m1_u a^+_m2_u operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdu (int m1, int m2, double& coefficient);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdd (int m1, int m2, double& coefficient);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAdd (int m1, int m2, double& coefficient);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int* spinIndices, int nbrIndices);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient);

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
 
   // convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // initialState = state to transform  
  // targetState = vector where the transformed state has to be stored
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // firstComponent = index of the first component to compute in initialState
  // nbrComponents = number of consecutive components to compute
  virtual void TransformOneBodyBasis(RealVector& initialState, RealVector& targetState, RealMatrix* oneBodyBasis, long firstComponent = 0l, long nbrComponents = 0l);

 // convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // initialState = state to transform  
  // targetState = vector where the transformed state has to be stored
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // firstComponent = index of the first component to compute in initialState
  // nbrComponents = number of consecutive components to compute
  virtual void TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, long firstComponent = 0l, long nbrComponents = 0l);

  // compute the transformation matrix from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // return value = transformation matrix
  virtual ComplexMatrix TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis);

  // compute the projection matrix from the SU(2) Hilbert space to an U(1) Hilbert space
  // 
  // targetSpace = pointer to the U(1) Hilbert space
  // type = type of particles that has to be kept (0 for type up, 1 for type up-minus, 2 for type down, 3 for type down-minus)
  // return value = projection matrix
  virtual ComplexMatrix TransformationMatrixSU2ToU1(BosonOnSphereShort* targetSpace, int type = 0);

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component as to be normalized to 1
  // symmetryFactor = if true also remove the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);    

  // convert a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // symmetryFactor = if true also add the symmetry factors
  // return value = converted state
  //  virtual RealVector& ConvertFromUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // normalize Jack with respect to cylinder basis
  //
  // state = reference to the Jack state to normalize
  // aspect = aspect ratio of cylinder
  // return value = normalized state
  virtual RealVector& NormalizeJackToCylinder(RealVector& state, double aspect);

  // normalize from the cylinder geometry to the Jack normalization
  //
  // state = reference to the state to unnormalize
  // aspect = cylinder aspect ratio
  // reference = set which component as to be normalized to 1
  // return value = unnormalized state
  virtual RealVector& NormalizeCylinderToJack(RealVector& state, double aspect, long reference = 0l);

  // normalize a state defined on the sphere geometry with respect to cylinder basis
  //
  // state = reference to the state to normalize
  // aspect = aspect ratio of cylinder
  // return value = normalized state
  virtual RealVector& NormalizeSphereToCylinder(RealVector& state, double aspect);

  // Compute the product of a spinful fermionic state with a Van der Monde determinant 
  //
  // fermionicState = reference on the spinful fermionic state
  // outputVector = reference on the vector where the result will be stored
  // fermionicSpace = pointer on the Hilbert Space associated to the spinful fermionic state
  // minIndex = first component to compute
  // nbrComponents = number of components to compute
  // unnormalizedFlag = true if the state should be written in the unnormalized basis
  // cylinderFlag = true if the state should be written on the cylinder geometry
  // cylinderPerimeter = cylinder perimeter
  // architecture = pointer to the architecture
  virtual void SlaterTimeSpinfulFermionicState(RealVector& fermionicState, RealVector& outputVector, FermionOnSphereWithSpin* fermionicSpace, 
					       int minIndex, int nbrComponents, bool unnormalizedFlag, bool cylinderFlag, double cylinderPerimeter, 
					       AbstractArchitecture* architecture);

  // Compute the product of a spinful fermionic state with a Van der Monde determinant 
  //
  // fermionicState = reference on the spinful fermionic state
  // outputVector = reference on the vector where the result will be stored
  // fermionicSpace = pointer on the Hilbert Space associated to the spinful fermionic state
  // minIndex = first component to compute
  // nbrComponents = number of components to compute
  // unnormalizedFlag = true if the state should be written in the unnormalized basis
  // cylinderFlag = true if the state should be written on the cylinder geometry
  // cylinderPerimeter = cylinder perimeter
  // architecture = pointer to the architecture
  virtual void SlaterTimeSpinfulFermionicState(RealVector& fermionicState, RealVector& outputVector, FermionOnSphereWithSpinTwoLandauLevels* fermionicSpace, 
					       int minIndex, int nbrComponents, bool unnormalizedFlag, bool cylinderFlag, double cylinderPerimeter,
					       AbstractArchitecture* architecture);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrParticleSector, int lzSector, int szSector, RealVector& groundState);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrParticleSector, int lzSector, int szSector, RealVector& groundState);

  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, RealVector& groundState, 
									 bool removeBinomialCoefficient = false, AbstractArchitecture* architecture = 0);
   
  // core part of the entanglement matrix evaluation for the particle partition
  // 
  // minIndex = first index to consider in the complementary Hilbert space
  // nbrIndex = number of indices to consider in the complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // entanglementMatrix = pointer to entanglement matrix
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = number of components that have been added to the entanglement matrix
  virtual long EvaluatePartialEntanglementMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace, 
								       ParticleOnSphere* destinationHilbertSpace, RealVector& groundState, RealMatrix* entanglementMatrix, 
								       bool removeBinomialCoefficient = false);
   
  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // nbrOrbitalA = number of orbitals that have to be kept for the A part
  // weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
  // weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
  // nbrOrbitalB = number of orbitals that have to be kept for the B part
  // weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
  // weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
  // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
  virtual RealMatrix& EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector,
													 int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
													 int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, RealMatrix& entanglementMatrix);

  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition breaking the momentum conservation. 
  // The entanglement matrix is computed from precalculated particle entanglement matrices in each momentum sector
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // nbrOrbitalA = number of orbitals that have to be kept for the A part
  // nbrConnectedOrbitalAUp = number of orbitals connected to a given one by the A part real space cut (for the spin up)
  // nbrConnectedOrbitalADown = number of orbitals connected to a given one by the A part real space cut (for the spin down)
  // connectedOrbitalAUp = orbitals taht connected to a given one by the A part real space cut (for the spin up)
  // connectedOrbitalADown = orbitals taht connected to a given one by the A part real space cut (for the spin down)
  // weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
  // weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
  // nbrOrbitalB = number of orbitals that have to be kept for the B part
  // nbrConnectedOrbitalBUp = number of orbitals connected to a given one by the B part real space cut (for the spin up)
  // nbrConnectedOrbitalBDown = number of orbitals connected to a given one by the B part real space cut (for the spin down)
  // connectedOrbitalBUp = orbitals taht connected to a given one by the B part real space cut (for the spin up)
  // connectedOrbitalBDown = orbitals taht connected to a given one by the B part real space cut (for the spin down)
  // weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
  // weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
  // nbrEntanglementMatrices = number of available entanglement matrices with a fixed moementum
  // entanglementMatrixLzSectors = momentum sector of each entanglement matrix
  // entanglementMatrices = array containing the entanglement matrices with a fixed moementum
  // return value = real space entanglement matrix
  virtual RealMatrix EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int szSector,
													int nbrOrbitalA, int* nbrConnectedOrbitalAUp, int* nbrConnectedOrbitalADown,
													int** connectedOrbitalAUp, int** connectedOrbitalADown, 
													double** weightOrbitalAUp, double** weightOrbitalADown, 
													int nbrOrbitalB, int* nbrConnectedOrbitalBUp, int* nbrConnectedOrbitalBDown, 
													int** connectedOrbitalBUp, int** connectedOrbitalBDown, 
													double** weightOrbitalBUp, double** weightOrbitalBDown, 
													int nbrEntanglementMatrices, int* entanglementMatrixLzSectors,
													RealMatrix* entanglementMatrices);

  // convert a given state from a generic basis from the current Sz subspace basis
  //
  // state = reference on the vector to convert
  // space = reference on the basis associated to state
  // return value = converted vector
  virtual RealVector ConvertToNbodyBasis(RealVector& state, ParticleOnSphereWithSpin* space);
  
  // convert a given state from a generic basis to the current Sz subspace basis
  //
  // state = reference on the vector to convert
  // space = reference on the basis associated to state
  // return value = converted vector
  virtual RealVector ConvertFromNbodyBasis(RealVector& state, ParticleOnSphereWithSpin* space);
  
  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // symmetrized a product of two decoupled states, core part
  //
  // outputVector = reference on the vector which will contain the symmetrized state
  // leftVector = reference on the vector associated to the first color
  // rightVector = reference on the vector associated to the second color
  // leftSpace = pointer to the Hilbert space of the first color
  // rightSpace = pointer to the Hilbert space of the second color
  // unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
  // firstComponent = index of the first component
  // nbrComponents = number of components to symmetrize
  // return value = symmetrized state
  virtual void SymmetrizeSU2SU2StateCore (RealVector& symmetrizedVector, RealVector& leftVector, RealVector& rightVector, 
					  ParticleOnSphereWithSpin* leftSpace, ParticleOnSphereWithSpin* rightSpace, 
					  bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents);

  // convert state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
  //
  // state = state that needs to be projected
  // su2space = SU(2) space with fixed sz of the input state
  // return value = input state expression in the SU(2) basis
  virtual RealVector SU2ToSU2AllSz(RealVector& state, ParticleOnSphereWithSpin* su2space);

  // convert state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
  //
  // state = state that needs to be projected
  // su2space = SU(2) space with fixed sz of the input state
  // return value = input state expression in the SU(2) basis
  virtual ComplexVector SU2ToSU2AllSz(ComplexVector& state, ParticleOnSphereWithSpin* su2space);

 protected:
  
  // read Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description is stored
  // return value = true if no error occured
  virtual bool ReadHilbertSpace (char* fileName);

  // find state index
  //
  // stateDescriptionUp = unsigned integer describing the fermionic state for type up particles
  // stateDescriptionDown = unsigned integer describing the fermionic state for type down particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // nbrNUp = number of particles with quantum number up
  // nbrNDown = number of particles with quantum number down
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int nbrNUp, int nbrNDown);

  // evaluate Hilbert space dimension without the Sz constraint
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMaxUp = momentum maximum value for a boson in the state with up
  // lzMaxDown = momentum maximum value for a boson in the state with down
  // totalLz = momentum total value
  // nbrNUp = number of particles with up
  // nbrNDown = number of particles with down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int lzMaxUp, int lzMaxDown, int totalLz, 
			      int nbrNUp, int nbrNDown, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // currentFermionicPositionUp = current fermionic position within the state description for the spin up
  // currentFermionicPositionDown = current fermionic position within the state description for the spin down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int lzMax, int totalLz, 
			      int currentFermionicPositionUp, int currentFermionicPositionDown, long pos);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initial bosonic state is stored
  // return value = corresponding fermionic state
  virtual unsigned long BosonToFermion(unsigned long*& initialState);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialStateUp = reference on the array where initial bosonic state for the type up particles is stored
  // initialStateDown = reference on the array where initial bosonic state for the type down particles is stored
  // finalStateUp = reference on the corresponding fermionic state for the type up particles
  // finalStateDown = reference on the corresponding fermionic state for the type down particles
  virtual void BosonToFermion(unsigned long*& initialStateUp, unsigned long*& initialStateDown, 
			      unsigned long& finalStateUp, unsigned long& finalStateDown);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax= maximum lz value reached by the fermionic state
  // finalState = reference on the array where the bosonic state for the type up particles has to be stored
  virtual void FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialStateUp = initial fermionic state for the type up particles
  // initialStateDown = initial fermionic state for the type down particles
  // finalStateUp = reference on the array where the bosonic state for the type up particles has to be stored
  // finalStateDown = reference on the array where the bosonic state for the type down particles has to be stored
  virtual void FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown, 
			      unsigned long*& finalStateUp, unsigned long*& finalStateDown);

  // convert a bosonic state to its monomial representation
  //
  // initialStateUp = initial spin up bosonic state in its fermionic representation 
  // initialStateDown = initial spin down bosonic state in its fermionic representation
  // finalStateUp = reference on the array where the monomial representation for the spin up has to be stored
  // finalStateDown = reference on the array where the monomial representation for the spin up has to be stored
  virtual void ConvertToMonomial(unsigned long initialStateUp, unsigned long initialStateDown, unsigned long*& finalStateUp, unsigned long*& finalStateDown);

  // convert a bosonic state from its monomial representation for single componentw
  //
  // initialStateUp = array where the monomial representation for the up spin is stored
  // initialStateDown = array where the monomial representation for the down spin is stored
  // finalStateUp = reference on the bosonic up state in its fermionic representation
  // finalStateDown = reference on the bosonic down state in its fermionic representation
  virtual void ConvertFromMonomial(unsigned long* initialStateUp, unsigned long* initialStateDown, 
				   unsigned long& finalStateUp, unsigned long& finalStateDown);

  // convert a bosonic state from its monomial representation for single component
  //
  // initialStateUp = array where the monomial representation for the up spin is stored
  // initialStateDown = array where the monomial representation for the down spin is stored
  // finalStateUp = reference on the bosonic up state
  // finalStateDown = reference on the bosonic down state
  virtual void ConvertFromMonomial(unsigned long* initialStateUp, unsigned long* initialStateDown, 
				   unsigned long*& finalStateUp, unsigned long*& finalStateDown);

  // convert a bosonic state from its monomial representation for single component
  //
  // initialState = array where the monomial representation for a single spin component is stored
  // finalState = reference on the bosonic single spin component state
  // nbrParticles = numebr of particles in the single spin component
  virtual void ConvertFromMonomial(unsigned long* initialState, unsigned long*& finalState, int nbrParticles);

  // apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // temporaryStatei= reference on the temporary array for the type i particles
  // temporaryStatej= reference on the temporary array for the type j particles
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state
  virtual int AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej, double& coefficient);

  // find state index
  //
  // stateDescriptionUp = array describing the bosonic state for type up particles
  // stateDescriptionDown = array describing the bosonic state for type down particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long*& stateDescriptionUp, unsigned long*& stateDescriptionDown);

  // recursive part of the convertion from a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // targetState = vector where the transformed state has to be stored
  // coefficient = current coefficient to assign
  // position = current particle consider in the n-body state
  // momentumIndices = array that gives the momentum partition of the initial n-body state
  // initialSU2Indices = array that gives the spin dressing the initial n-body state
  // currentSU2Indices = array that gives the spin dressing the current transformed n-body state
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
  // occupationCoefficientArray = array that provides 1/2 ln (N!)
  virtual void TransformOneBodyBasisRecursive(RealVector& targetState, double coefficient,
					      int position, int* momentumIndices, int* initialSU2Indices, int* currentSU2Indices, RealMatrix* oneBodyBasis, 
					      double occupationCoefficient, double* occupationCoefficientArray);
  
  // recursive part of the convertion from a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // targetState = vector where the transformed state has to be stored
  // coefficient = current coefficient to assign
  // position = current particle consider in the n-body state
  // momentumIndices = array that gives the momentum partition of the initial n-body state
  // initialSU2Indices = array that gives the spin dressing the initial n-body state
  // currentSU2Indices = array that gives the spin dressing the current transformed n-body state
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
  // occupationCoefficientArray = array that provides 1/2 ln (N!)
  virtual void TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
					      int position, int* momentumIndices, int* initialSU2Indices, int* currentSU2Indices, ComplexMatrix* oneBodyBasis, 
					      double occupationCoefficient, double* occupationCoefficientArray);

  // Compute the product of a spinful Slater determinant with a Van der Monde determinant
  //
  // slaterUp = monomial representation of the Slater spin up part
  // slaterDown = monomial representation of the Slater spin up part
  // finalState = reference on the vector the produced state will be stored
  // threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
  virtual void VanDerMondeTimesSlater (unsigned long* slaterUp, unsigned long* slaterDown, RealVector& finalState, double** threeOrbitalOverlaps);
  
  // Compute the product of a spinful Slater determinant with a Van der Monde determinant, assuming a reverse flux attachment
  //
  // slaterUp = monomial representation of the Slater spin up part
  // slaterDown = monomial representation of the Slater spin up part
  // finalState = reference on the vector the produced state will be stored
  // threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
  virtual void ReverseVanDerMondeTimesSlater (unsigned long* slaterUp, unsigned long* slaterDown, RealVector& finalState, double** threeOrbitalOverlaps);
  
  // Compute the product of a spinful Slater determinant with a Van der Monde determinant, assuming a reverse flux attachment and peforming only (N-1)! permutations
  //
  // slaterUp = monomial representation of the Slater spin up part
  // slaterDown = monomial representation of the Slater spin up part
  // finalState = reference on the vector the produced state will be stored
  // threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
  // position = perform a swap between the last element the one at given position 
  virtual void ReverseVanDerMondeTimesSlater (unsigned long* slaterUp, unsigned long* slaterDown, RealVector& finalState, 
					      double** threeOrbitalOverlaps, int position);

  // Compute the product of a spinful Slater determinant in two Landau levels with a Van der Monde determinant
  //
  // slaterLLLUp = monomial representation of the lowest Landau part of the Slater spin up part
  // slater2LLUp = monomial representation of the second Landau part of the Slater spin up part
  // slaterLLLDown = monomial representation of the lowest Landau part  of the Slater spin down part
  // slater2LLDown = monomial representation of the second Landau part of the Slater spin down part
  // nbrBosonsLLLUp - number of spin up bosons in the lowest Landau level
  // nbrBosonsLLLDown - number of spin down bosons in the lowest Landau level
  // finalState = reference on the vector the produced state will be stored
  // threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
  virtual void VanDerMondeTimesSlater (unsigned long* slaterLLLUp, unsigned long* slater2LLUp, unsigned long* slaterLLLDown, unsigned long* slater2LLDown, 
				       int nbrBosonsLLLUp, int nbrBosonsLLLDown, RealVector& finalState, double*** threeOrbitalOverlaps);
  
  // Compute the product of a spinful Slater determinant in two Landau levels with a Van der Monde determinant, assuming a reverse flux attachment
  //
  // slaterLLLUp = monomial representation of the lowest Landau part of the Slater spin up part
  // slater2LLUp = monomial representation of the second Landau part of the Slater spin up part
  // slaterLLLDown = monomial representation of the lowest Landau part  of the Slater spin down part
  // slater2LLDown = monomial representation of the second Landau part of the Slater spin down part
  // nbrBosonsLLLUp - number of spin up bosons in the lowest Landau level
  // nbrBosonsLLLDown - number of spin down bosons in the lowest Landau level
  // finalState = reference on the vector the produced state will be stored
  // threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
  virtual void ReverseVanDerMondeTimesSlater (unsigned long* slaterLLLUp, unsigned long* slater2LLUp, unsigned long* slaterLLLDown, unsigned long* slater2LLDown, 
					      int nbrBosonsLLLUp, int nbrBosonsLLLDown, RealVector& finalState, double*** threeOrbitalOverlaps);
  
  // Compute the product of a spinful Slater determinant in two Landau levels with a Van der Monde determinant, assuming a reverse flux attachment
  //
  // slaterLLLUp = monomial representation of the lowest Landau part of the Slater spin up part
  // slater2LLUp = monomial representation of the second Landau part of the Slater spin up part
  // slaterLLLDown = monomial representation of the lowest Landau part  of the Slater spin down part
  // slater2LLDown = monomial representation of the second Landau part of the Slater spin down part
  // nbrBosonsLLLUp - number of spin up bosons in the lowest Landau level
  // nbrBosonsLLLDown - number of spin down bosons in the lowest Landau level
  // finalState = reference on the vector the produced state will be stored
  // threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
  // position = perform a swap between the last element the one at given position 
  virtual void ReverseVanDerMondeTimesSlater (unsigned long* slaterLLLUp, unsigned long* slater2LLUp, unsigned long* slaterLLLDown, unsigned long* slater2LLDown, 
					      int nbrBosonsLLLUp, int nbrBosonsLLLDown, RealVector& finalState, double*** threeOrbitalOverlaps, int position);

  
  // compute the integral \int phi_{k1,0}^* phi_{k2,0} phi_{k3,0} on the cylinder geometry where phi_{k,n} in wavefunction with momentum k and int the n-th LL
  // 
  // k1 = momentum of the first orbital for the phi_{k1,n} wavefunction
  // l1 = square value of magnetic length for the phi_{k1,n} wavefunction
  // k2 = momentum of the first orbital for the phi_{k2,n} wavefunction
  // l2 = square value of magnetic length for the phi_{k2,n} wavefunction
  // k3 = momentum of the first orbital for the phi_{k3,n} wavefunction
  // l3 = square value of magnetic length for the phi_{k3,n} wavefunction
  // perimeter = cylinder perimeter
  // return value = corresponding integral
  virtual double ComputeIntegralPhi0Phi0Phi0OnCylinder(double k1, double l1, double k2, double l2, double k3, double l3, double perimeter);

  // compute the number of particles in a given state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of particles
  virtual int ComputeNbrParticles(unsigned long stateDescription);

  };

// get the number of orbitals
//
// return value = number of orbitals

inline int BosonOnSphereWithSU2Spin::GetNbrOrbitals()
{
  return this->NbrLzValue;
}

// get the number of particles
//
// return value = number of particles

inline int BosonOnSphereWithSU2Spin::GetNbrParticles()
{
  return this->NbrBosons;
}

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnSphereWithSU2Spin::GetParticleStatistic()
{
  return AbstractQHEParticle::BosonicStatistic;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initial bosonic state is stored
// return value = corresponding fermionic state

inline unsigned long BosonOnSphereWithSU2Spin::BosonToFermion(unsigned long*& initialState)
{
  unsigned long FinalState = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      FinalState |= ((1ul << initialState[i]) - 1ul) << Shift;
      Shift += initialState[i];
      ++Shift;
    }
  return FinalState;
}

// convert a bosonic state into its fermionic counterpart
//
// initialStateUp = reference on the array where initial bosonic state for the type up particles is stored
// initialStateDown = reference on the array where initial bosonic state for the type down particles is stored
// finalStateUp = reference on the corresponding fermionic state for the type up particles
// finalStateDown = reference on the corresponding fermionic state for the type down particles

inline void BosonOnSphereWithSU2Spin::BosonToFermion(unsigned long*& initialStateUp, unsigned long*& initialStateDown,
						     unsigned long& finalStateUp, unsigned long& finalStateDown)
{
  finalStateUp = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateUp |= ((1ul << initialStateUp[i]) - 1ul) << Shift;
      Shift += initialStateUp[i];
      ++Shift;
    }
  finalStateDown = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateDown |= ((1ul << initialStateDown[i]) - 1ul) << Shift;
      Shift += initialStateDown[i];
      ++Shift;
    }
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateLzMax= maximum lz value reached by the fermionic state
// finalState = reference on the array where the bosonic state for the type up particles has to be stored

inline void BosonOnSphereWithSU2Spin::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState)
{
  int FinalStateLzMax = 0;
  while ((initialStateLzMax >= 0) && ((initialState >> initialStateLzMax) == 0x0ul))
    --initialStateLzMax;
  while (initialStateLzMax >= 0)
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
      finalState[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState >>= TmpPower;
      ++FinalStateLzMax;
      initialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalState[FinalStateLzMax] = 0x0ul;
}


// convert a fermionic state into its bosonic  counterpart
//
// initialStateUp = initial fermionic state for the type up particles
// initialStateDown = initial fermionic state for the type down particles
// finalStateUp = reference on the array where the bosonic state for the type up particles has to be stored
// finalStateDown = reference on the array where the bosonic state for the type down particles has to be stored

inline void BosonOnSphereWithSU2Spin::FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown,
						     unsigned long*& finalStateUp, unsigned long*& finalStateDown)
{
  int FinalStateLzMax = 0;
  int InitialStateLzMax = this->NUpLzMax;
  while ((InitialStateLzMax >= 0) && ((initialStateUp >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateUp - 1ul) ^ (~initialStateUp);
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
      finalStateUp[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateUp >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalStateUp[FinalStateLzMax] = 0x0ul;

  FinalStateLzMax = 0;
  InitialStateLzMax = this->NDownLzMax;
  while ((InitialStateLzMax >= 0) && ((initialStateDown >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateDown - 1ul) ^ (~initialStateDown);
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
      finalStateDown[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateDown >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalStateDown[FinalStateLzMax] = 0x0ul;
}

// convert a bosonic state to its monomial representation
//
// initialStateUp = initial spin up bosonic state in its fermionic representation 
// initialStateDown = initial spin down bosonic state in its fermionic representation
// finalStateUp = reference on the array where the monomial representation for the spin up has to be stored
// finalStateDown = reference on the array where the monomial representation for the spin up has to be stored

inline void BosonOnSphereWithSU2Spin::ConvertToMonomial(unsigned long initialStateUp, unsigned long initialStateDown,
							unsigned long*& finalStateUp, unsigned long*& finalStateDown) 
{
  int InitialStateLzMax = this->NUpLzMax;
  int Index = 0;
  int TmpLz = this->NUpLzMax - this->NbrBosonsUp + 1;
  while (InitialStateLzMax >= 0)
    {
      while ((InitialStateLzMax >= 0) && (((initialStateUp >> InitialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalStateUp[Index++] = TmpLz;
	  --InitialStateLzMax;
	}
      while ((InitialStateLzMax >= 0) && (((initialStateUp >> InitialStateLzMax) & 0x1ul) == 0x0ul))
	{
	  --TmpLz;
	  --InitialStateLzMax;
	}
    }
  InitialStateLzMax = this->NDownLzMax;
  TmpLz = this->NDownLzMax - this->NbrBosonsDown + 1;
  Index = 0;
  while (InitialStateLzMax >= 0)
    {
      while ((InitialStateLzMax >= 0) && (((initialStateDown >> InitialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalStateDown[Index++] = TmpLz;
	  --InitialStateLzMax;
	}
      while ((InitialStateLzMax >= 0) && (((initialStateDown >> InitialStateLzMax) & 0x1ul) == 0x0ul))
	{
	  --TmpLz;
	  --InitialStateLzMax;
	}
    }
}

// convert a bosonic state from its monomial representation for single componentw
//
// initialStateUp = array where the monomial representation for the up spin is stored
// initialStateDown = array where the monomial representation for the down spin is stored
// finalStateUp = reference on the bosonic up state in its fermionic representation
// finalStateDown = reference on the bosonic down state in its fermionic representation

inline void BosonOnSphereWithSU2Spin::ConvertFromMonomial(unsigned long* initialStateUp, unsigned long* initialStateDown, 
							  unsigned long& finalStateUp, unsigned long& finalStateDown)
{
  finalStateUp = 0x0ul;
  for (int i = 0; i < this->NbrBosonsUp; ++i)
    finalStateUp |= 0x1ul << (initialStateUp[i] + ((unsigned long) (this->NbrBosonsUp - i)) - 1ul);
  finalStateDown = 0x0ul;
  for (int i = 0; i < this->NbrBosonsDown; ++i)
    finalStateDown |= 0x1ul << (initialStateDown[i] + ((unsigned long) (this->NbrBosonsDown - i)) - 1ul);
}

// convert a bosonic state from its monomial representation for single component
//
// initialStateUp = array where the monomial representation for the up spin is stored
// initialStateDown = array where the monomial representation for the down spin is stored
// finalStateUp = reference on the bosonic up state
// finalStateDown = reference on the bosonic down state

inline void BosonOnSphereWithSU2Spin::ConvertFromMonomial(unsigned long* initialStateUp, unsigned long* initialStateDown, 
							  unsigned long*& finalStateUp, unsigned long*& finalStateDown)
{
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateUp[i] = 0ul;
      finalStateDown[i] = 0ul;
    }
  for (int i = 0; i < this->NbrBosonsUp; ++i)
    {
      ++finalStateUp[initialStateUp[i]]; 
    }
  for (int i = 0; i < this->NbrBosonsDown; ++i)
    {
      ++finalStateDown[initialStateDown[i]]; 
    }
}

// convert a bosonic state from its monomial representation for single component
//
// initialState = array where the monomial representation for a single spin component is stored
// finalState = reference on the bosonic single spin component state
// nbrParticles = numebr of particles in the single spin component

inline void BosonOnSphereWithSU2Spin::ConvertFromMonomial(unsigned long* initialState, unsigned long*& finalState, int nbrParticles)
{
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalState[i] = 0ul;
    }
  for (int i = 0; i < nbrParticles; ++i)
    {
      ++finalState[initialState[i]]; 
    }
}

// apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// temporaryStatei= reference on the temporary array for the type i particles
// temporaryStatej= reference on the temporary array for the type j particles
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int BosonOnSphereWithSU2Spin::AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej,
					     double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryStateUp[i] = this->ProdATemporaryStateUp[i];
      this->TemporaryStateDown[i] = this->ProdATemporaryStateDown[i];
    }
  ++temporaryStatej[m2];
  coefficient = temporaryStatej[m2];
  ++temporaryStatei[m1];
  coefficient *= temporaryStatei[m1];
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(2) index for the first annihilation operator
// sigma2 = SU(2) index for the second annihilation operator
// return value =  multiplicative factor 

inline double BosonOnSphereWithSU2Spin::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateSigma[sigma1][n1] == 0) || (this->ProdATemporaryStateSigma[sigma2][n2] == 0) || 
      ((n1 == n2) && (sigma1 == sigma2) && (this->ProdATemporaryStateSigma[sigma1][n1] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryStateSigma[sigma2][n2];
  --this->ProdATemporaryStateSigma[sigma2][n2];
  Coefficient *= this->ProdATemporaryStateSigma[sigma1][n1];
  --this->ProdATemporaryStateSigma[sigma1][n1];
  return sqrt(Coefficient);
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(3) index for the first creation operator
// sigma2 = SU(3) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int BosonOnSphereWithSU2Spin::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryStateUp[i] = this->ProdATemporaryStateUp[i];
      this->TemporaryStateDown[i] = this->ProdATemporaryStateDown[i];
    }
  ++this->TemporaryStateSigma[sigma2][m2];
  coefficient = this->TemporaryStateSigma[sigma2][m2];
  ++this->TemporaryStateSigma[sigma1][m1];
  coefficient *= this->TemporaryStateSigma[sigma1][m1];
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
}

// find state index
//
// stateDescriptionUp = array describing the bosonic state for type up particles
// stateDescriptionDown = array describing the bosonic state for type down particles
// return value = corresponding index

inline int BosonOnSphereWithSU2Spin::FindStateIndex(unsigned long*& stateDescriptionUp, unsigned long*& stateDescriptionDown)
{
  unsigned long Tmp1;
  unsigned long Tmp2;
  this->BosonToFermion(stateDescriptionUp, stateDescriptionDown, Tmp1, Tmp2);
  return this->FindStateIndex(Tmp1, Tmp2);
}

// compute the number of particles in a given state
//
// stateDescription = unsigned integer describing the state
// return value = number of particles

inline int BosonOnSphereWithSU2Spin::ComputeNbrParticles(unsigned long stateDescription)
{
  unsigned long TmpNbrParticle = 0l;
  while(stateDescription != 0x0ul) 
    {
      TmpNbrParticle += (stateDescription & 0x1ul);
      stateDescription >>= 1;
    }
  return (int) TmpNbrParticle;
}

#endif


