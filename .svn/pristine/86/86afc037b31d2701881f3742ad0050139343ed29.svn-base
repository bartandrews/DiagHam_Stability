////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          hardcore 3-body interaction                       //
//                                                                            //
//                        last modification : 30/04/2014                      //
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


#ifndef PARTICLEONTORUSGENERICTHREEBODYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H
#define PARTICLEONTORUSGENERICTHREEBODYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian.h"

#include <iostream>


using std::ostream;



class ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian : public AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian
{

 protected:

  // array where the interaction coefficients are stored
  double****** PrecalculatedInteractionCoefficients;

 public:

  // default constructor
  //
  ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // twoBodyDeltaStrength = strength of the additional two body delta interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio,
								      AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
								      char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian();
  
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
	
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a+_m3 a_n1 a_n2 a_n3 coupling term
  //
  // m1 = first creation operator index
  // m2 = second creation operator index
  // m3 = third creation operator index
  // n1 = first annihilation operator index
  // n2 = second annihilation operator index
  // n3 = thrid annihilation operator index
  //
  // return value = numerical coefficient  
  virtual double EvaluateInteractionCoefficient(int m1, int m2, int m3, int n1, int n2, int n3);
  
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  virtual double EvaluateTwoBodyInteractionCoefficient(int m1, int m2, int m3, int m4);

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a+_m3 a_n1 a_n2 a_n3 coupling term, symmetrized on the n indices
  //
  // m1 = first creation operator index
  // m2 = second creation operator index
  // m3 = third creation operator index
  // n1 = first annihilation operator index
  // n2 = second annihilation operator index
  // n3 = thrid annihilation operator index
  // return value = numerical coefficient
  inline double EvaluateInteractionNIndexSymmetrizedCoefficient(int m1, int m2, int m3, int n1,int n2, int n3);

};

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a+_m3 a_n1 a_n2 a_n3 coupling term, symmetrized on the n indices
//
// m1 = first creation operator index
// m2 = second creation operator index
// m3 = third creation operator index
// n1 = first annihilation operator index
// n2 = second annihilation operator index
// n3 = thrid annihilation operator index
// return value = numerical coefficient

inline double ParticleOnTorusThreeBodyHardcoreWithMagneticTranslationsHamiltonian::EvaluateInteractionNIndexSymmetrizedCoefficient(int m1, int m2, int m3, int n1,int n2, int n3)
{  
  double TmpInteraction;
  if (n1 != n2)
    {
      if (n2 != n3)
	{	  
	  TmpInteraction  = this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n2][n3];
	  TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n3][n2];
	  TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][m2][m3][n2][n1][n3];
	  TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][m2][m3][n2][n3][n1];
	  TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][m2][m3][n3][n1][n2];
	  TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][m2][m3][n3][n2][n1];
	}
      else
	{	  
	  TmpInteraction  = this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n2][n3];
	  TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][m2][m3][n2][n1][n3];
	  TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][m2][m3][n3][n2][n1];
	}
    }
  else
    {
      if (n2 != n3)
	{	  
	  TmpInteraction  = this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n2][n3];
	  TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n3][n2];
	  TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][m2][m3][n2][n1][n3];
	}
      else
	{	  
	  TmpInteraction  = this->PrecalculatedInteractionCoefficients[m1][m2][m3][n1][n2][n3];
	}
    }
  return TmpInteraction;
}

#endif
