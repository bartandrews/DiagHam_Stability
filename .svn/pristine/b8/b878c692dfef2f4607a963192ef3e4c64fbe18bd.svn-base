////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                         hollowcore n-body interaction                      //
//                                                                            //
//                        last modification : 17/01/2015                      //
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


#ifndef PARTICLEONTORUSNBODYHOLLOWCOREWITHMAGNETICTRANSLATIONSHAMILTONIAN_H
#define PARTICLEONTORUSNBODYHOLLOWCOREWITHMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian.h"
#include "MathTools/FactorialCoefficient.h" 

#include <iostream>
#include <algorithm>


using std::ostream;



class ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian : public ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian
{

 protected:

 public:

  // default constructor
  //
  ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // nbrNBody = value of the n (i.e. the n-body interaction)
  // nBodyStrength = strength of the n-body interaction
  // regenerateElementFlag = regenerate th interaction matrix elements instead of reading them from the harddrive
  // nbrPseudopotentials = number of two-body pseudopotentials (can be zero for a pure n-body interaction)
  // pseudopotentials = pseudopotential coefficients for the two-body interaction 
  // oneBodyPotentials = array of additional one-body potentials
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, 
								    int nbrNBody, double nBodyStrength, bool regenerateElementFlag, int nbrPseudopotentials,
								    double* pseudopotentials, double* oneBodyPotentials,
								    AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
								    char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian();
  
 protected:
  
  // evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term
  //
  // creationCoefficients = array that contains the creation coefficients
  // annihilationCoefficients = array that contains the annihilation coefficients
  // return value = numerical coefficient  
  virtual double EvaluateInteractionCoefficient(double* creationCoefficients, double* annihilationCoefficients);
  
  // evaluate all the numerical coefficients in front of the \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation operators only) for each integer modulo the NBodyValue (without using symmetries)
  //
  void EvaluateInteractionCoefficientCreation();
  
  // evaluate all the numerical coefficients in front of the \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation operators only) for each integer modulo the NBodyValue (using symmetries)
  //
  void EvaluateInteractionCoefficientCreationUsingSymmetries();
  
  // evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation operators only) for each integer modulo the NBodyValue
  //
  // mIndices = array that contains the creation indices 
  // g = modulo of the indices that are being summed on
  // momentumTransfer = momentum transfer operated by the \prod_i a+_mi \prod_j a_nj operator, in units of the number of flux quanta
  // return value = array of numerical coefficients 
  double EvaluateIndividualInteractionCoefficientCreation(int* mIndices, int g, int momentumTransfer);
  
  // evaluate the two nested Gaussian sum for a three body interaction (for test purposes)
  //
  // momFactor = array of indices that contains the information of the creation (or annihilation) indices
  // TmpIndices = array of indices that gives the initial indices that will be incremented in the sum
  // return value = value of the sum
  double EvaluateGaussianSum(int* momFactor, int* TmpIndices);
  
  
  // evaluate the N nested infinite sums of EvaluateInteractionCoefficientCreation
  //
  // nBodyValue = current index being incremented 
  //TmpIndices = array containing the current value of all indices
  // Sum = current value of the sum
  // countIter = array of integers giving, for each TmpIndices[i], the number of times it has been incremented
  // momFactor = array of indices that contains the information of the creation (or annihilation) indices
  //return value = value of the coefficient
  double EvaluateGaussianSum(int nBodyValue, int* TmpIndices, double Sum, int* countIter, int* momFactor);
  
};

#endif
