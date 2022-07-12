////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of particle on sphere                        //
//                                                                            //
//                        last modification : 15/07/2002                      //
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


#ifndef PARTICLEONSPHERE_H
#define PARTICLEONSPHERE_H


#include "config.h"
#include "MathTools/Complex.h"
#include "HilbertSpace/AbstractQHEParticle.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/HermitianMatrix.h"


class LongRationalVector;
class LongRational;
class LongRationalMatrix;
class AbstractArchitecture;
class SparseRealMatrix;
class SparseComplexMatrix;


class ParticleOnSphere :  public AbstractQHEParticle 
{

 public:

  enum 
    {
      NoSymmetry = 0x0,
      LzMinusLzSymmetry = 0x1
    };


  // virtual destructor
  //
  virtual ~ParticleOnSphere ();

  // get the number of orbitals
  //
  // return value = number of orbitals
  virtual int GetNbrOrbitals();

  // get the number of particles
  //
  // return value = number of particles
  virtual int GetNbrParticles();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic() = 0;

  // get information about any additional symmetry of the Hilbert space
  //
  // return value = symmetry id
  virtual int GetHilbertSpaceAdditionalSymmetry();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphere* targetSpace);

  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

  // get the number of particles in the target space
  //
  // return value = number of particles in the target space
  virtual int GetTargetNbrParticles();

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2), safe version i.e. works with any numbers of particles
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAASafe (int index, int m1, int m2, int n1, int n2, double& coefficient);

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

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual long AdAdAA (long index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m1 a^+_m2 a^+_m3 a_n1 a_n2 a_n3 operator to a given state (with m1+m2+m3=n1+n2+n3)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // m3 = third index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // n3 = third index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAdAAA (int index, int m1, int m2, int m3, int n1, int n2, int n3, double& coefficient);

  // apply a^+_m1 a^+_m2 a^+_m3 a^+_m4 a_n1 a_n2 a_n3 a_n4 operator to a given state (with m1+m2+m3+m4=n1+n2+n3+n4)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // m3 = third index for creation operator
  // m4 = fourth index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // n3 = third index for annihilation operator
  // n4 = fourth index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAdAdAAAA (int index, int m1, int m2, int m3, int m4, int n1, int n2, int n3, int n4, double& coefficient);

  // apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
  //
  // index = index of the state on which the operator has to be applied
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient);
  
  // apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2);

  // apply a_n1 a_n2 operator to a given state without keeping it in cache
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AA (int index, int n1, int n2, double& coefficient);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next ProdA call
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

  // apply Prod_i a_ni operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdA (int index, int* n, int nbrIndices, double& coefficient);

  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient);

  // apply a^+_m1 a^+_m2 operator to the state 
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int index, int m1, int m2, double& coefficient);

  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to apply to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int nbrIndices, double& coefficient);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to apply to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int nbrIndices, double& coefficient, int& nbrTranslation);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int nbrIndices, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  // use double when calculating normalization factors to avoid overflow
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAdL (int* m, int nbrIndices, double& coefficient);

  // apply Prod_i a^+_ni operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int index, int* n, int nbrIndices, double& coefficient);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (int index, int m);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient);

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
  virtual long AdA (long index, int m, int n, double& coefficient);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual long AdA (long index, int m, int n, Complex& coefficient);

  // apply creation operator to a word, using the conventions
  // for state-coding and quantum numbers of this space
  // state = word to be acted upon
  // m = Lz value of particle to be added
  // coefficient = reference on the double where the multiplicative factor has to be stored
  virtual unsigned long Ad (unsigned long state, int m, double& coefficient);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient, int& nbrTranslation);
    
  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // apply a_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A (int index, int n);

  // apply a^+_m  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
  //
  // index = index of the state on which the operator has to be applied
  // m = index for annihilation operator
  // return value =  multiplicative factor 
  virtual double Ad (int index, int m);

  // apply a_n operator to the state produced using the A or Ad method (without destroying it)
  //
  // n = first index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int A (int n, double& coefficient);

  // apply a^+_m operator to the state produced using the A or Ad method (without destroying it)
  //
  // m = first index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad (int m, double& coefficient);

  // apply a^+_m operator to the state produced using A method (without destroying it)
  //
  // m = first index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int Ad (int m, double& coefficient, int& nbrTranslationX);

  // apply a^+_m operator to the state produced using A method (without destroying it)
  //
  // m = first index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int Ad (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
    
  // apply a_n  operator to a given state. 
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value =  index of the resulting state 
  virtual int A (int index, int n, double& coefficient);

  // apply a^+_n  operator to a given state. 
  //
  // index = index of the state on which the operator has to be applied
  // n = index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value =  index of the resulting state 
  virtual int Ad (int index, int n, double& coefficient);

  // check whether HilbertSpace implements ordering of operators
  //
  virtual bool HaveOrder ();
  
  // check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
  virtual int CheckOrder (int* m, int* n, int nbrIndices);

  // evaluate wave function in real space using a given basis
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis);

  // evaluate wave function in real space using a given basis, using time coherence
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // nextCoordinates = index of the coordinate that will be changed during the next time iteration
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
							 AbstractFunctionBasis& basis, int nextCoordinates);
  
  // evaluate wave function in real space using a given basis and only for agiven range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
					int firstComponent, int nbrComponent);                                
  
  // evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // nextCoordinates = index of the coordinate that will be changed during the next time iteration
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
							 AbstractFunctionBasis& basis, 
							 int nextCoordinates, int firstComponent, int nbrComponent);

  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  virtual void InitializeWaveFunctionEvaluation (bool timeCoherence = false);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, ComplexVector& groundState);
 
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
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particle. The geometrical cut is a stripe.
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax+shitedCut to -Lzmax+shitedCut+subsytemSize-1)
  // shiftedCut = first orbital belonging to the subsystem (with angular momentum -Lzmax+shitedCut)
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluateShiftedPartialDensityMatrix (int subsytemSize, int nbrShiftedOrbitals, int nbrBosonSector, int lzSector, RealVector& groundState);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particle. The geometrical cut is a stripe.
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax+shitedCut to -Lzmax+shitedCut+subsytemSize-1)
  // shiftedCut = first orbital belonging to the subsystem (with angular momentum -Lzmax+shitedCut)
  // nbrBosonSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluateShiftedPartialDensityMatrix (int subsytemSize, int nbrShiftedOrbitals, int nbrBosonSector, int lzSector, ComplexVector& groundState);

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

  // core part of the evaluation density matrix particle partition calculation involving a sum of projetors 
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
								  int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix);

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
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixRealSpacePartition (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using real space partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // perimeter = cylinder perimeter
  // height = height of a cylinder (from -H/2 to H/2) 
  // xcut = x-coordinate of a cylinder cut
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixRealSpacePartitionCylinder (int nbrBosonSector, int lzSector, double perimeter, double height, double xcut, RealVector& groundState, AbstractArchitecture* architecture = 0);

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
  virtual long EvaluatePartialDensityMatrixRealSpacePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
								   RealVector& groundState,  RealSymmetricMatrix* densityMatrix, double* incompleteBetaThetaBottom, double* incompleteBetaThetaTop, double phiRange);

  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, bool removeBinomialCoefficient = false);

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
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = pointer to an array containing the total system ground states
  // nbrEntanglementMatrices = number of vectors whose entanglement matrix has to be computed
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix* EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector* groundState, int nbrEntanglementMatrices, bool removeBinomialCoefficient = false);

  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and both A and B are resticted to a given number of orbitals
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated
  // nbrOrbitalA = number of orbitals that have to be kept for the A part (starting from the leftmost orbital)
  // nbrOrbitalA = number of orbitals that have to be kept for the B part (starting from the rightmost orbital)
  // groundState = pointer to an array containing the total system ground states
  // nbrEntanglementMatrices = number of vectors whose entanglement matrix has to be computed
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix* EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int lzSector, int nbrOrbitalA, int nbrOrbitalB, 
									RealVector* groundState, int nbrEntanglementMatrices, bool removeBinomialCoefficient);
  
  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, ComplexVector& groundState, bool removeBinomialCoefficient = false);

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
  virtual ComplexMatrix EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int lzSector, int nbrOrbitalA, int nbrOrbitalB, 
									ComplexVector& groundState, bool removeBinomialCoefficient);
  
    // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = pointer to an array containing the total system ground states
  // nbrEntanglementMatrices = number of vectors whose entanglement matrix has to be computed
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix* EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, ComplexVector* groundState, int nbrEntanglementMatrices, bool removeBinomialCoefficient = false);

  // evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and both A and B are resticted to a given number of orbitals
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated
  // nbrOrbitalA = number of orbitals that have to be kept for the A part (starting from the leftmost orbital)
  // nbrOrbitalA = number of orbitals that have to be kept for the B part (starting from the rightmost orbital)
  // groundState = pointer to an array containing the total system ground states
  // nbrEntanglementMatrices = number of vectors whose entanglement matrix has to be computed
  // removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
  // return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix* EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int lzSector, int nbrOrbitalA, int nbrOrbitalB, 
									ComplexVector* groundState, int nbrEntanglementMatrices, bool removeBinomialCoefficient);

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
  virtual RealMatrix& EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, RealMatrix& entanglementMatrix);
  
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
  virtual ComplexMatrix& EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, ComplexMatrix& entanglementMatrix);

  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition on a cylinder. The entanglement matrix is only evaluated in a given Lz sector.
  // and computed from precalculated particle entanglement matrix
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // perimeter = cylinder perimeter
  // height = height of a cylinder (from -H/2 to H/2) 
  // xcut = x-coordinate of the cut   // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
  virtual RealMatrix& EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrixCylinder (int nbrBosonSector, int lzSector, double perimeter, double height, double xcut, RealMatrix& entanglementMatrix);

  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrOrbitalA = number of orbitals that have to be kept for the A part
  // weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
  // nbrOrbitalB = number of orbitals that have to be kept for the B part
  // weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
  // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
  virtual RealMatrix& EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, 
													 int nbrOrbitalA, double* weightOrbitalA, 
													 int nbrOrbitalB, double* weightOrbitalB, RealMatrix& entanglementMatrix);
  
   // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrOrbitalA = number of orbitals that have to be kept for the A part
  // weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
  // nbrOrbitalB = number of orbitals that have to be kept for the B part
  // weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
  // entanglementMatrix = pointer to array of entanglement matrix (will be overwritten)
  // nbrEntanglementMatrices = number of entanglement matrices to be computed
  // return value = reference on the entanglement matrix
  virtual RealMatrix* EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, 
													 int nbrOrbitalA, double* weightOrbitalA, 
													 int nbrOrbitalB, double* weightOrbitalB, RealMatrix* entanglementMatrix, int nbrEntanglementMatrices);
  
  // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrOrbitalA = number of orbitals that have to be kept for the A part
  // weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
  // nbrOrbitalB = number of orbitals that have to be kept for the B part
  // weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
  // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
  // return value = reference on the entanglement matrix
  virtual ComplexMatrix& EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, 
													 int nbrOrbitalA, double* weightOrbitalA, 
													 int nbrOrbitalB, double* weightOrbitalB, ComplexMatrix& entanglementMatrix);
  
   // evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
  // The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
  // 
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // nbrOrbitalA = number of orbitals that have to be kept for the A part
  // weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
  // nbrOrbitalB = number of orbitals that have to be kept for the B part
  // weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
  // entanglementMatrix = pointer to array of entanglement matrix (will be overwritten)
  // nbrEntanglementMatrices = number of entanglement matrices to be computed
  // return value = reference on the entanglement matrix
  virtual ComplexMatrix* EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, 
													 int nbrOrbitalA, double* weightOrbitalA, 
													 int nbrOrbitalB, double* weightOrbitalB, ComplexMatrix* entanglementMatrix, int nbrEntanglementMatrices);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector
  // and computed from precalculated entanglement matrix
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // entanglementMatrix = reference on the entanglement matrix
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartitionFromEntanglementMatrix (int nbrBosonSector, int lzSector, RealMatrix& entanglementMatrix);

  // compute part of the Schmidt decomposition, allowing cut in the reduced denisty matrix eigenvalue space
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // eigenvalueCut = discard all contribution from the reduced density matrix whose eigenvalues is lower than eigenvalueCut
  // rebuiltSchmidtGroundState = reference on the state to whose current sector contribution to the Schmidt decomposition will be added 
  // diagonalizedDensityMatrix = reference on the diagonalized reduced density matrix associated to the current sector (with down ordered diagonal elements)
  // transformationMatrix =  reference on the transformation matric that diagonalizes the reduced density matrix
  // return value = reference on rebuiltSchmidtGroundState
  virtual RealVector& EvaluatePartialSchmidtDecomposition(int subsytemSize, int nbrFermionSector, int lzSector, double eigenvalueCut,
							  RealVector& groundState, RealVector& rebuiltSchmidtGroundState,
							  RealDiagonalMatrix& diagonalizedDensityMatrix, RealMatrix& transformationMatrix);

  // compute part of the Schmidt decomposition for the particle partition, allowing cut in the reduced denisty matrix eigenvalue space
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // eigenvalueCut = discard all contribution from the reduced density matrix whose eigenvalues is lower than eigenvalueCut
  // rebuiltSchmidtGroundState = reference on the state to whose current sector contribution to the Schmidt decomposition will be added 
  // diagonalizedDensityMatrix = reference on the diagonalized reduced density matrix associated to the current sector (with down ordered diagonal elements)
  // transformationMatrix =  reference on the transformation matric that diagonalizes the reduced density matrix
  // return value = reference on rebuiltSchmidtGroundState
  virtual RealVector& EvaluatePartialSchmidtDecompositionParticlePartition(int nbrParticleSector, int lzSector, double eigenvalueCut,
									   RealVector& groundState, RealVector& rebuiltSchmidtGroundState,
									   RealDiagonalMatrix& diagonalizedDensityMatrix, RealMatrix& transformationMatrix);

  // rebuild a state from its Schmidt decomposition for the particle partition
  // 
  // nbrParticleSector = number of particles that belong to the subsytem (i.e. part A)
  // lzSector = Lz sector in which the density matrix has to be evaluated  (i.e. part A)
  // schmidtDecomposedState = reference on the vector to which the rebuild state will be added
  // nbrSingularValues = number of singular values (can be lower than the actual number of ingular values to perform a truncation)
  // singularValues = array containing the singular values
  // aVectors = matrix than contains the singular vectors of the part A
  // bVectors = transposed matrix than contains the singular vectors of the part B
  virtual void RebuildStateFromSchmidtDecompositionParticlePartition(int nbrParticleSector, int lzSector, RealVector& schmidtDecomposedState, 
									     int nbrSingularValues, double* singularValues, RealMatrix& aVectors, RealMatrix& bVectors);

  // convert a state to its occupation number representation
  //
  // index = index of the state
  // finalState = reference on the array where the occupation number representation has to be stored
  virtual void GetOccupationNumber(long index, unsigned long*& finalState);

  // get the list of occupied orbitals in a given state
  //
  // state = ID of the state
  // orbitals = list of orbitals to be filled
  virtual void GetOccupied(int state, int* orbitals);

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
 
  // convert a state such that its components are now expressed in the normalized basis, shifting all orbitals
  //
  // state = reference to the state to convert
  // shift = shift to apply to each orbitals
  // reference = set which component has been normalized to 1
  // return value = converted state
  virtual RealVector& ShiftedConvertFromUnnormalizedMonomial(RealVector& state, int shift, long reference = 0);

  // convert a state such that its components, given in the conformal limit,  are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // return value = converted state
  virtual RealVector& ConvertFromConformalLimit(RealVector& state, long reference);

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
					  ParticleOnSphere** inputSpaces, bool symmetrizedFlag = false, double coefficient = 1.0);

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
  virtual double JackSqrNormalization (RealVector& state, long minIndex = 0l, long nbrComponents = 0l);

  // compute part of the Jack polynomial square normalization in a given range of indices
  //
  // state = reference on the unnormalized Jack polynomial
  // minIndex = first index to compute 
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization 
  virtual LongRational JackSqrNormalization (LongRationalVector& state, long minIndex = 0l, long nbrComponents = 0l);

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
  
  // remove part of each Fock state, discarding component if the Fock state does not a given pattern
  //
  // inputVector = state to truncate
  // reducedSpace = Hilbert space where the truncated state will lie
  // pattern = array describing the pattern 
  // patternSize = pattern size
  // patternShift = indicate where the pattern has to be applied
  // return value = trucated state
  virtual RealVector TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnSphere* reducedSpace, int* pattern, int patternSize, int patternShift = 0);

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual bool HasPauliExclusions(int index, int pauliK, int pauliR);
  
  // get Lz component of a component
  //
  // j = index of the component in Hilbert space
  // return value = twice the Lz component
  virtual int GetLzValue(int j=0);

  // get Sz component of the spin
  //
  // j = index of the vector in Hilbert space
  // return value = Sz component
  virtual int GetSzValue(int j);

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

  // evaluate coeffecicents requested to compute the real space partition
  //
  // lzMax = twice the maximum angular momentum
  // thetaTop = inclination angle defining one edge of the cut in radians
  // thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in radians
  // incompleteBetaThetaTop = reference on the pointer to the array (allocation done by the method) where the top part coefficients will be stored
  // incompleteBetaThetaBotton = reference on the pointer to the array (allocation done by the method) where the bottom part coefficients will be stored
  virtual void EvaluatePartialDensityMatrixRealSpacePartitionCoefficient(int lzMax, double thetaTop, double thetaBottom, double*& incompleteBetaThetaTop, double*& incompleteBetaThetaBottom);

  // evaluate coeffecicents requested to compute the real space partition (cylinder geometry)
  //
  // lzMax = twice the maximum angular momentum
  // perimeter = cylinder perimeter
  // xcut = x-coordinate of the cut
  // incompleteBetaThetaTop = reference on the pointer to the array (allocation done by the method) where the top part coefficients will be stored
  virtual void EvaluatePartialDensityMatrixRealSpacePartitionCoefficientCylinder(int lzMax, double perimeter, double xcut, double*& incompleteBetaThetaTop);

  // compute the number of particles in each Landau level
  //
  // state = ID of the state to handle
  // lLOccupationConfiguration = array where the decomposition will be store
  virtual void LandauLevelOccupationNumber(int state, int* lLOccupationConfiguration);

  virtual void EvaluatePartialDensityMatrixMultipartiteParticlePartition(ParticleOnSphere * spaceA, ParticleOnSphere * spaceB, ParticleOnSphere * spaceC,  RealVector groundstate, RealSymmetricMatrix* densityMatrix, AbstractArchitecture* architecture = 0);

  // convert the vector with a given Lz to the full space (all Lz components)
  // inputState = input vector
  // inputSpace = input Hilbert space with given Lz
  // return value = vector in the full Hilbert space
  void ConvertToAllLz (ComplexVector& inputState, ParticleOnSphere* inputSpace, ComplexVector& outputState);

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

  // create a state from its MPS description
  //
  // bMatrices = array that gives the B matrices 
  // state = reference to vector that will contain the state description
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // memory = amount of memory that can be use to precompute matrix multiplications  
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, RealVector& state, int mPSRowIndex, int mPSColumnIndex, 
					      long memory = 0l, long initialIndex = 0l, long nbrComponents = 0l);

  // create a state from its MPS description, inclusing additional quasihole matrices
  //
  // bMatrices = array that gives the B matrices 
  // state = reference to vector that will contain the state description
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // memory = amount of memory that can be use to precompute matrix multiplications  
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, SparseComplexMatrix* quasiholeBMatrices, int nbrQuasiholeBMatrices,
					      ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
					      long memory = 0l, long initialIndex = 0l, long nbrComponents = 0l);

  // create a state from its site-dependent MPS description
  //
  // bMatrices = array that gives the site-dependent MPS
  // state = reference to vector that will contain the state description
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromSiteDependentMPSDescription (SparseRealMatrix** bMatrices, RealVector& state, int mPSRowIndex, int mPSColumnIndex, 
							   long initialIndex = 0l, long nbrComponents = 0l);

  //get the dimension of the subspace generated by application of all symmetry to a state stored as a bosonic state in TemporaryState 
  //
  //state = ID in the symmetrized basis of the space whose symmetry has to be determined
  //return value = dimension of the corresponding subspace  
  virtual int GetSymmetryDimension(int state);
  
  // convert a state defined in the Ky basis into a state in the (Kx,Ky) basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertToKxKyBasis(ComplexVector& state, ParticleOnSphere* space);

  // convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space);

  // apply a Gutzwiller projection (in the orbital space) to a given state
  //
  // state = reference on the state to project
  // space = pointer to the Hilbert space where state is defined
  // return value = Gutzwiller projected state
  virtual ComplexVector GutzwillerProjection(ComplexVector& state, ParticleOnSphere* space);

  // Compute the overlap of two states made from different one-body wavefunction
  //
  //  firstVector = reference on the first vector 
  //  secondVector = reference on the second vector 
  // overlapMatrix = pointer to the table with the overlap between the one-body states
  virtual Complex ComputeOverlapWaveFunctionsWithDifferentGamma (ComplexVector& firstVector, ComplexVector& secondVector, Complex * overlapMatrix);

  // compute sum of positions in the x and y direction for lattice class
  //
  // index = index of the state in the basis whose position sums are to be computed
  // positionX = reference on the sum of positions in the x direction
  // positionY = reference on the sum of positions in the y direction
  virtual void GetPositionSum(int index,int & positionX, int & positionY);
  virtual void GetPositionSum(unsigned long * monomial, int & positionX, int & positionY);
  
};

// convert the Lz value from the sphere geometry to the disk geometry
// 
// lzValue = twice the Lz value on the sphere geometry
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = Lz on the disk geometry

inline int ConvertLzFromSphereToDisk (int lzValue, int nbrParticles, int nbrFluxQuanta)
{
  return ((lzValue + (nbrParticles * nbrFluxQuanta)) >> 1);
}

// convert the Lz value from the disk geometry to the sphere geometry
// 
// lzValue = Lz value on the disk geometry
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// return value = twice the Lz on the sphere geometry

inline int ConvertLzFromDiskToSphere (int lzValue, int nbrParticles, int nbrFluxQuanta)
{
  return ((2 * lzValue) - (nbrParticles * nbrFluxQuanta));
}

// get the number of orbitals
//
// return value = number of orbitals

inline int ParticleOnSphere::GetNbrOrbitals()
{
  return -1;
}

// get the number of particles
//
// return value = number of particles

inline int ParticleOnSphere::GetNbrParticles()
{
  return -1;
}

#endif
