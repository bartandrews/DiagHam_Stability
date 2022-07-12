////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//     class for bosons on sphere with SU(2) spin and Sz<->-Sz symmetry       //
//                                                                            //
//                        last modification : 23/09/2016                      //
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


#ifndef BOSONONSPHEREWITHSU2SPINSZSYMMETRY_H
#define BOSONONSPHEREWITHSU2SPINSZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"

#include <iostream>

using std::cout;
using std::endl;



class BosonOnSphereWithSU2SpinSzSymmetry :  public BosonOnSphereWithSU2Spin
{


 protected:


  // sign of the parity sector for the Sz<->-Sz symmetry
  double SzParitySign;
  
  // a temporary variable to store the normalization factor coming from the orbit size
  int ProdATemporaryNbrStateInOrbit;

  // pointer to the Hilbert space where the result of any operator lies
  BosonOnSphereWithSU2SpinSzSymmetry* TargetSpace;

  // maximum number of states in a given orbit
  int MaxOrbitSize;
  // number of state in each orbit
  int* NbrStateInOrbit;
  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;

 public:

  // default constructor
  // 
  BosonOnSphereWithSU2SpinSzSymmetry ();

  // basic constructor without any constraint on Sz
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  BosonOnSphereWithSU2SpinSzSymmetry (int nbrBosons, int totalLz, int lzMax, bool minusSzParity);

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // totalSpin = twice the total spin value (not taken into account)
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2SpinSzSymmetry (int nbrBosons, int totalLz, int lzMax, int totalSpin, bool minusSzParity, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  // 
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2SpinSzSymmetry (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSU2SpinSzSymmetry(const BosonOnSphereWithSU2SpinSzSymmetry& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSU2SpinSzSymmetry ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSU2SpinSzSymmetry& operator = (const BosonOnSphereWithSU2SpinSzSymmetry& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphereWithSpin* targetSpace);

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

  // apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // sigma1 = SU(2) index for the first annihilation operator
  // sigma2 = SU(2) index for the second annihilation operator
  // return value =  multiplicative factor 
  virtual double AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2);

  // apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // sigma1 = SU(2) index for the first creation operator
  // sigma2 = SU(2) index for the second creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient);

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
   
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
  // The entanglement matrix is only evaluated in a given Lz,Sz=0, Sz parity sectors.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated
  // szSector = Sz sector in which the density matrix has to be evaluated. It should be equal to zero
  // szParitySector = parity sector for the discrete symmetry Sz<->-Sz
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, int szParity, RealVector& groundState, 
									 bool removeBinomialCoefficient = false, AbstractArchitecture* architecture = 0);
   
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
													 int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, 
													 RealMatrix& entanglementMatrix);

  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
  // 
  // nbrParticleSector = number of particles that belong to the subsystem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // szParitySector = parity sector for the discrete symmetry Sz<->-Sz. Can be either -1 or +1
  // nbrOrbitalA = number of orbitals that have to be kept for the A part
  // weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
  // weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
  // nbrOrbitalB = number of orbitals that have to be kept for the B part
  // weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
  // weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
  // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
  virtual RealMatrix& EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector, int szParity,
													 int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
													 int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, 
													 RealMatrix& entanglementMatrix);

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

  protected:

  // read Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description is stored
  // return value = true if no error occured
  virtual bool ReadHilbertSpace (char* fileName);

  // generate the Hilbert space with the discrete symmetry constraint
  //
  virtual void GenerateStatesWithDiscreteSymmetry();

  // find state index
  //
  // stateDescriptionUp = unsigned integer describing the fermionic state for type up particles
  // stateDescriptionDown = unsigned integer describing the fermionic state for type down particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown);

  // compute the rescaling factors
  //
  virtual void ComputeRescalingFactors();

  // find state index
  //
  // stateDescriptionUp = array describing the bosonic state for type up particles
  // stateDescriptionDown = array describing the bosonic state for type down particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long*& stateDescriptionUp, unsigned long*& stateDescriptionDown);

  // apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // temporaryStatei= reference on the temporary array for the type i particles
  // temporaryStatej= reference on the temporary array for the type j particles
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state
  virtual int AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej, double& coefficient);

  // factorized code that is used to symmetrize the result of any operator action
  //
  // stateDescriptionUp = reference on the state up part that has been produced with the operator action
  // stateDescriptionDown = reference on the state down part that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown, double& coefficient);

  // factorized code that is used to symmetrize the result of any operator action
  //
  // stateDescriptionUp = reference on the state up part that has been produced with the operator action
  // stateDescriptionDown = reference on the state down part that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(unsigned long*& stateDescriptionUp, unsigned long*& stateDescriptionDown, double& coefficient);

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
   

};

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(2) index for the first annihilation operator
// sigma2 = SU(2) index for the second annihilation operator
// return value =  multiplicative factor 

inline double BosonOnSphereWithSU2SpinSzSymmetry::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateSigma[sigma1][n1] == 0) || (this->ProdATemporaryStateSigma[sigma2][n2] == 0) || 
      ((n1 == n2) && (sigma1 == sigma2) && (this->ProdATemporaryStateSigma[sigma1][n1] == 1)))
    {
      return 0.0;
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = this->ProdATemporaryStateSigma[sigma2][n2];
  --this->ProdATemporaryStateSigma[sigma2][n2];
  Coefficient *= this->ProdATemporaryStateSigma[sigma1][n1];
  --this->ProdATemporaryStateSigma[sigma1][n1];
  return sqrt(Coefficient);
}

// apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// temporaryStatei= reference on the temporary array for the type i particles
// temporaryStatej= reference on the temporary array for the type j particles
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int BosonOnSphereWithSU2SpinSzSymmetry::AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej,
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
  return this->SymmetrizeAdAdResult(this->TemporaryStateUp, this->TemporaryStateDown, coefficient);
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(3) index for the first creation operator
// sigma2 = SU(3) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int BosonOnSphereWithSU2SpinSzSymmetry::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
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
  return this->SymmetrizeAdAdResult(this->TemporaryStateUp, this->TemporaryStateDown, coefficient);
}

// find state index
//
// stateDescriptionUp = array describing the bosonic state for type up particles
// stateDescriptionDown = array describing the bosonic state for type down particles
// return value = corresponding index

inline int BosonOnSphereWithSU2SpinSzSymmetry::FindStateIndex(unsigned long*& stateDescriptionUp, unsigned long*& stateDescriptionDown)
{
  unsigned long Tmp1;
  unsigned long Tmp2;
  this->BosonToFermion(stateDescriptionUp, stateDescriptionDown, Tmp1, Tmp2);
  return this->FindStateIndex(Tmp1, Tmp2);
}

// factorized code that is used to symmetrize the result of any operator action
//
// stateDescriptionUp = reference on the state up part that has been produced with the operator action
// stateDescriptionDown = reference on the state down part that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state  

inline int BosonOnSphereWithSU2SpinSzSymmetry::SymmetrizeAdAdResult(unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown, double& coefficient)
{
  if (stateDescriptionUp < stateDescriptionDown)
    {
      unsigned long Tmp = stateDescriptionUp;
      stateDescriptionUp = stateDescriptionDown;
      stateDescriptionDown = Tmp;
      coefficient *= this->SzParitySign * this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][2];
    }
  else
    {
      if (stateDescriptionUp != stateDescriptionDown)
	{
	  coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][2];
	}
      else
	{
	  coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][1];
	}
    }
  return this->TargetSpace->FindStateIndex(stateDescriptionUp, stateDescriptionDown);
}

// factorized code that is used to symmetrize the result of any operator action
//
// stateDescriptionUp = reference on the state up part that has been produced with the operator action
// stateDescriptionDown = reference on the state down part that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state  

inline int BosonOnSphereWithSU2SpinSzSymmetry::SymmetrizeAdAdResult(unsigned long*& stateDescriptionUp, unsigned long*& stateDescriptionDown, double& coefficient)
{
  unsigned long Tmp1;
  unsigned long Tmp2;
  this->BosonToFermion(stateDescriptionUp, stateDescriptionDown, Tmp1, Tmp2);
  return this->SymmetrizeAdAdResult(Tmp1, Tmp2, coefficient);
}

#endif


