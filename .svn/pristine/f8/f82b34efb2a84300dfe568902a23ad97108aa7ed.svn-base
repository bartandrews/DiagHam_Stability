////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of hamiltonian associated to particles on a             //
//        4D space  torus x sphere with a generic two-body interaction        //
//                          and magnetic translations                         //
//                                                                            //
//                        last modification : 24/02/2017                      //
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


#ifndef PARTICLEONT2XS2WITHMAGNETICTRANSLATIONSGENERICTWOBODYHAMILTONIAN_H
#define PARTICLEONT2XS2WITHMAGNETICTRANSLATIONSGENERICTWOBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnTorusWithMagneticTranslationsHamiltonian.h"
#include "Polynomial/Polynomial.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;


class ParticleOnT2xS2WithMagneticTranslationsGenericTwoBodyHamiltonian : public AbstractQHEOnTorusWithMagneticTranslationsHamiltonian
{

 protected:

  // number of flux quanta for the torus
  int NbrFluxQuantumTorus;
  // number of flux quanta for the sphere
  int NbrFluxQuantumSphere;
  // number of orbitals on the sphere
  int NbrLzValues;

  // number of pseudo-potentials
  int NbrPseudoPotentials;
  // pseudo-potential torus relative momentum
  int* PseudoPotentialMomentumTorus;
  // pseudo-potential sphere relative angular momenta
  int* PseudoPotentialMomentumSphere;
  //pseudo-potential coefficients
  double* PseudoPotentials;

  // Laguerre polynomial for the pseudopotentials
  Polynomial* LaguerrePolynomials;

 public:

  // default constructor
  //
  ParticleOnT2xS2WithMagneticTranslationsGenericTwoBodyHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuantumTorus = number of flux quanta piercing the torus
  // kxMomentum = momentum in the x direction for the torus
  // maxMomentum = maximum Lz value reached by a particle in the state
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // nbrFluxQuantumSphere = number of flux quanta piercing the sphere
  // nbrPseudoPotentials = number of pseudo-potentials
  // pseudoPotentialMomentumTorus = torus pseudo-potential relative momenta
  // pseudoPotentialMomentumSphere = sphere pseudo-potential relative angular momenta
  // pseudoPotentials = pseudo-potential coefficients
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnT2xS2WithMagneticTranslationsGenericTwoBodyHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int nbrFluxQuantumTorus, int xMomentum,
								   double ratio, int nbrFluxQuantumSphere, int nbrPseudoPotentials, int* pseudoPotentialMomentumTorus, int* pseudoPotentialMomentumSphere, 
								   double* pseudoPotentials, AbstractArchitecture* architecture, long memory = -1, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnT2xS2WithMagneticTranslationsGenericTwoBodyHamiltonian();

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

 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // evaluate the numerical coefficient in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term for the torus part
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // nbrFluxQuanta = number of flux quanta
  // ratio = ratio between the width in the x direction and the width in the y direction
  // pseudopotentialMomentum = pseudo-potential relative momentum
  // return value = numerical coefficient
  virtual double EvaluateTorusInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrFluxQuanta, 
						     double ratio, int pseudopotentialMomentum);

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term for the sphere [art
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // clebsch = reference on the Clebsch-Gordan coefficients for the sphere geometry
  // nbrFluxQuanta = number of flux quanta
  // pseudopotentialMomentum = pseudo-potential relative momentum
  // return value = numerical coefficient
  virtual double EvaluateSphereInteractionCoefficient(int m1, int m2, int m3, int m4, ClebschGordanCoefficients& clebsch,
						      int nbrFluxQuanta,int pseudopotentialMomentum);

  // evaluate the numerical coefficient  in front of the a+_(lz1,kz1) a+_(lz2,kz2) a_(lz3,kz3) a_(lz4,kz4) coupling term
  //
  // ky1 = first ky index
  // ky2 = second ky index
  // ky3 = third ky index
  // ky4 = fourth ky index
  // lz1 = first lz index
  // lz2 = second lz index
  // lz3 = third lz index
  // lz4 = fourth lz index
  // clebsch = reference on the Clebsch-Gordan coefficients for the sphere geometry
  // pseudopotentialIndex1 = pseudopotential index for the interaction on the torus 
  // pseudopotentialIndex2 = pseudopotential index for the interaction on the sphere
  // pseudoPotential = pseudopotential amplitude
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int ky1, int ky2, int ky3, int ky4, 
						int lz1, int lz2, int lz3, int lz4, 
						ClebschGordanCoefficients& clebsch,
						int pseudopotentialIndex1, int pseudopotentialIndex2, double pseudoPotential);

  // get a linearized index from the two momenta
  //
  // ky = momentum along the y direction for the torus
  // lz = z projection of the angular momentum for the sphere
  // return value = linearized index 
  virtual int GetLinearizedIndex(int ky, int lz);

  // get the two momenta associated to a given linearized index
  //
  // index = linearized index 
  // ky = reference on the momentum along the y direction for the torus
  // lz = reference on the z projection of the angular momentum for the sphere
  virtual void GetLinearizedIndex(int index, int& ky, int& lz);

  // get all the indices that should appear in the annihilation/creation operators
  //
  virtual void GetIndices();

};

// get a linearized index from the two momenta
//
// ky = momentum along the y direction for the torus
// lz = z projection of the angular momentum for the sphere
// return value = linearized index 

inline int ParticleOnT2xS2WithMagneticTranslationsGenericTwoBodyHamiltonian::GetLinearizedIndex(int ky, int lz)
{
  return ((ky * this->NbrLzValues) + lz);
}

// get the two momenta associated to a given linearized index
//
// index = linearized index 
// ky = reference on the momentum along the y direction for the torus
// lz = reference on the z projection of the angular momentum for the sphere

inline void ParticleOnT2xS2WithMagneticTranslationsGenericTwoBodyHamiltonian::GetLinearizedIndex(int index, int& ky, int& lz)
{
  lz = index % this->NbrLzValues;
  ky = index / this->NbrLzValues;
}


#endif
