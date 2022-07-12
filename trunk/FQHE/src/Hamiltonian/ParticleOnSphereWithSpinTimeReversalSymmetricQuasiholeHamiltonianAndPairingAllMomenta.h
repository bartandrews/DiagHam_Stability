////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//       class of hamiltonian associated to Laughlin qh on a sphere with      //
//        SU(2) spin with opposite magnetic field for each species,           //
//                      pairing and no momentum conservation                  //
//                                                                            //
//                        last modification : 26/08/2016                      //
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


#ifndef PARTICLEONSPHEREWITHSPINTIMEREVERSALSYMMETRYQUASIHOLEHAMILTONIANANDPAIRINGALLMOMENTA_H
#define PARTICLEONSPHEREWITHSPINTIMEREVERSALSYMMETRYQUASIHOLEHAMILTONIANANDPAIRINGALLMOMENTA_H


#include "config.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "Hamiltonian/ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian.h"

#include <iostream>


using std::ostream;


class ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta : public ParticleOnLatticeTimeReversalBreakingSingleBandHamiltonian
{

   friend class QHEParticlePrecalculationOperation;

 protected:

   // maxixum monentum transfer that can appear in a one body operator
   int MaximumMomentumTransfer;

  // array that contains all one-body interaction factors for particles with spin up
  double* OneBodyInteractionFactorsupup;
  // array that contains all one-body interaction factors for particles with spin down
  double* OneBodyInteractionFactorsdowndown;
  // off-diagonal contribution of the one-body potential for particles with spin up, the first entry is the annihilation index, the second entry is the momentum tranfer
  Complex** OneBodyOffDiagonalInteractionFactorsupup;
  // off-diagonal contribution of the one-body potential for particles with spin down, the first entry is the annihilation index, the second entry is the momentum tranfer
  Complex** OneBodyOffDiagonalInteractionFactorsdowndown;


  // array that contains all one-body interaction factors for the pairing term
  Complex* OneBodyInteractionFactorsPairing;
  // off diagonal contribution of the one-body pairing term, the first entry is the index of the rightmost creation operator, the second entry is the momentum tranfer
  Complex** OneBodyOffDiagonalInteractionFactorsPairing; 

  // factor in front of the charging energy (i.e 1/(2C))
  double ChargingEnergy;
  // avearge number of particles in the system
  double AverageNumberParticles;
  // global shift to apply to the diagonal matrix elements
  //  double HamiltonianShift;
  



 public:

  ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // lzmax = maximum Lz value reached by a particle in the state
  // maxMomentumTransfer = maxixum monentum transfer that can appear in a one body operator
  // onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
  // onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
  // onebodyOffDiagonalPotentialUpUp = off-diagonal contribution of the one-body potential for particles with spin up, the first entry is the annihilation index, 
  //                                   the second entry is the momentum tranfer
  // onebodyOffDiagonalPotentialDownDown = off-diagonal contribution of the one-body potential for particles with spin down, the first entry is the annihilation index, 
  //                                       the second entry is the momentum tranfer
  // onebodyPotentialPairing = one-body pairing term (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
  // onebodyOffDiagonalPotentialPairing = off diagonal contribution of the one-body pairing term, the first entry is the index of the rightmost creation operator, 
  //	    				  the second entry is the momentum tranfer
  // chargingEnergy = factor in front of the charging energy (i.e 1/(2C))
  // averageNumberParticles = average number of particles in the system
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta(QuasiholeOnSphereWithSpinAndPairing* particles, int lzmax, int maxMomentumTransfer,
											double* onebodyPotentialUpUp, double* onebodyPotentialDownDown,
											Complex** onebodyOffDiagonalPotentialUpUp, Complex** onebodyOffDiagonalPotentialDownDown,
											Complex* onebodyPotentialPairing, Complex** onebodyOffDiagonalPotentialPairing, 
											double chargingEnergy, double averageNumberParticles,
											AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);
  
  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  virtual AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  virtual int GetHilbertSpaceDimension ();
  
   // ask if Hamiltonian implements hermitian symmetry operations
  //
  virtual bool IsHermitian();

 protected:
 
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						     int firstComponent, int nbrComponent);
  
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						      int firstComponent, int nbrComponent);
 
  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
							      int firstComponent, int nbrComponent);
  
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								int firstComponent, int nbrComponent);
  
  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent);

  
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& HermitianLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									 int firstComponent, int nbrComponent);
  
  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										 int firstComponent, int nbrComponent);
  
  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = number of components that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int nbrComponent);

  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = number of components that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int nbrComponent);
  
};


// ask if Hamiltonian implements hermitian symmetry operations
//
// return value = true if the Hamiltonian implements hermitian symmetry operations

inline bool ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::IsHermitian()
{
  return this->HermitianSymmetryFlag;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

inline int ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairingAllMomenta::GetHilbertSpaceDimension ()
{
  return (this->Particles->GetHilbertSpaceDimension());
}
  
  
#endif
