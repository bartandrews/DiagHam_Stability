////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a twisted torus      //
//                       with generic n-body interaction                      //
//                                                                            //
//                        last modification : 31/01/2015                      //
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


#ifndef PARTICLEONTWISTEDTORUSGENERICNBODYBODYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H
#define PARTICLEONTWISTEDTORUSGENERICNBODYBODYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian.h"

#include <iostream>


using std::ostream;



class ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian : public AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian
{

  friend class FQHETorusComputeMatrixElementOperation;

 protected:

  // angle (in radian) between the two fundamental cycles of the torus, along (L1 sin, L1 cos) and (0, L2)
  double Theta;
  // cosine of theta
  double CosTheta;
  // sine of theta
  double SinTheta;

  // number of monomials in the Fourier transformed interaction
  int NbrMonomials;
  // coefficients in front of each monomial in the Fourier transformed interaction
  double* MonomialCoefficients;
  // description of each monomial in the Fourier transformed interaction
  int** MonomialDescription;
  
 public:

  // default constructor
  //
  ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = number of flux quanta
  // xMomentum = relative angular momentum along x 
  // ratio = torus aspect ratio (Lx/Ly)
  // theta =  angle (in pi units) between the two fundamental cycles of the torus, along (Lx sin theta, Lx cos theta) and (0, Ly)
  // nbrNBody = type of interaction i.e. the number of density operators that are involved in the interaction
  // interactionName = name of the interaction, will be use to generate the interaction matrix element output file name
  // nbrMonomials = number of monomials in the Fourier transformed interaction
  // monomialCoefficients = coefficients in front of each monomial in the Fourier transformed interaction
  // monomialDescription = description of each monomial in the Fourier transformed interaction
  // regenerateElementFlag = regenerate th interaction matrix elements instead of reading them from the harddrive
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, double theta,
									int nbrNBody, char* interactionName, 
									int nbrMonomials, double* monomialCoefficients, int** monomialDescription, 
									bool regenerateElementFlag, AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
									char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian();
  
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
  virtual void EvaluateInteractionFactors();
	
  // evaluate the numerical coefficient  in front of the Prod a^+_mi Prod a+_n coupling term
  //
  // mIndices = array containing the creation operator indices
  // nIndices = array containing the annihilation operator indices
  // return value = numerical coefficient  
  //  virtual Complex EvaluateInteractionCoefficient(int* mIndices, int* nIndices, long& nbrOperations);
  virtual Complex EvaluateInteractionCoefficient(int* mIndices, int* nIndices, double* qxValues, double* qyValues, double* q2Values, double* cosineCoefficients);
  
  //  virtual Complex RecursiveEvaluateInteractionCoefficient(double currentSumQx, double currentSumQy, double currentSumQ2, double currentSumPhase, double& currentPrecision, double* qxValues, double* qyValues, double* q2Values, double* cosineCoefficients);

  virtual Complex RecursiveEvaluateInteractionCoefficient(int xPosition, double currentSumQx, double currentSumQy, double currentSumQ2, double currentSumPhase, double& currentPrecision, double* qxValues, double* qyValues, double* q2Values, double* cosineCoefficients);

  virtual double VFactor(double* q2Values);


  // read the interaction matrix elements from disk
  //
  // fileName = name of the file where the interaction matrix elements are stored
  // return value = true if no error occured
  virtual bool ReadInteractionFactors(char* fileName);

  // write the interaction matrix elements from disk
  //
  // fileName = name of the file where the interaction matrix elements are stored
  // return value = true if no error occured
  virtual bool WriteInteractionFactors(char* fileName);
  
  // test if creation/annihilation are in canonical form
  //
  // mIndices = array that contains the creation indices (will be modified)
  // nIndices = array that contains the annihilation indices (will be modified)
  // linearizedMIndices = linearized creation index
  // linearizedNIndices = linearized annihilation index
  // totalMomentum = momentum sector of the creation/annihilation indices
  // canonicalMIndices = reference on the linearized creation index of the canonical form
  // canonicalMIndices = reference on the linearized annihilation index of the canonical form
  // canonicalTotalMomentum = reference on the momentum sector of the creation/annihilation indices of the canonical form
  // canonicalPermutationCoefficient = additional permutation coefficient to get the canonical form
  // conjugationFlag = true if an additional conjugation should be applied to get the canonical form
  // return value = true if the indices are in canonical form
  virtual bool FindCanonicalIndices(int* mIndices, int* nIndices, int linearizedMIndices, int linearizedNIndices, 
				    int totalMomentum, int& canonicalMIndices, int& canonicalNIndices, 
				    int& canonicalTotalMomentum, double& canonicalPermutationCoefficient,
				    bool& conjugationFlag);

};


#endif
