////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of SkyrmionOnSphere state wave function                    //
//                                                                            //
//                        last modification : 20/04/2005                      //
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


#ifndef SKYRMIONONSPHEREWAVEFUNCTION_H
#define SKYRMIONONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "QHEWaveFunctionManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"
#include "Vector/RealVector.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "GeneralTools/GarbageFlag.h"

#ifdef USE_HILBERT_SPACE

class AbstractQHEParticle;
class AbstractFunctionBasis;


class SkyrmionOnSphereWaveFunction: public Abstract1DComplexTrialFunction
{

 protected:

  // Pointer to architecture
  AbstractArchitecture* Architecture;

  // Number of particles;
  int NbrParticles;

  // Nbr of flux in polarized part
  int PolarizedLzMax;

  // Angular momentum in polarized part
  int PolarizedLz;

  // data for bosonic hilbert-space
  int BosonLzMax;
  int BosonLz;
  int BosonSz;
  
  // vector that describes polarized state components in PolarizedSpace basis
  RealVector PolarizedState;

  // vector that describes bosonic spin texture components in BosonicSpace basis
  RealVector BosonicState;

  // Hilbert space associated to the polarized SkyrmionOnSphere state
  ParticleOnSphere* PolarizedSpace;

  // Hilbert space associated to the polarized SkyrmionOnSphere state
  ParticleOnSphereWithSpin *BosonicSpace;

  // one body real space basis to use for polarized state
  AbstractFunctionBasis* OneBodyBasisPol;
  // one body real space basis to use for bosonic state
  AbstractFunctionBasis* OneBodyBasisBos;

  // last value of the bosonic state
  Complex LastBosonicValue;
  
  // flag indicating whether we are using an exact polarized wavefunction
  bool UseExact;

  // Analytic polarized wavefunction
  Abstract1DComplexFunction* AnalyticPolarizedWaveFunction;

  // type cast of function pointer, if suitable
  Abstract1DComplexTrialFunction* AnalyticPolarizedTrialWaveFunction;

  // minimum dimension of bosonic space before parallelisation is used
  int MinParallel;

  // Garbage flag
  GarbageFlag Flag;

 public:

  // constructor
  //
  // create object to be initialized from system options
  //
  SkyrmionOnSphereWaveFunction(AbstractArchitecture* architecture, OptionManager &manager, int nbrParticles,
			       int totalLzMax, int totalLz, int totalSz, QHEWaveFunctionManager *wfManager,
			       int basisType=ParticleOnSphereFunctionBasis::LeftHanded);


  // copy constructor
  //
  // function = reference on the wave function to copy
  SkyrmionOnSphereWaveFunction(const SkyrmionOnSphereWaveFunction& function);

  // destructor
  //
   ~SkyrmionOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

  // get a value of the wavefunction for the last set of coordinates, but with different variational parameters
  // parameters =  alternative set of parameters
  virtual Complex GetForOtherParameters( double *parameters);

  // do many evaluations of the function, storing the result in the vector results given in the call
  // result = vector of leading dimension of the array coefficients for returning values
  // x = positions to evaluate the wavefuntion in
  // format for passing parameters as [nbrSet][nbrParameter],
  virtual void GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients);

  // access internal values of parameters
  virtual double *GetTrialParameters();

  // get number of parameters
  virtual int GetNbrParameters();
  
  // set new values of the trial coefficients (keeping the initial number of parameters)
  virtual void SetTrialParameters(double * coefficients);

  // test parity under reversal of all spins and rotation
  void TestSymmetries(ParticleOnSphereCollection *particles);

  // get function properties, and possible extensions of interface 
  // 
  virtual unsigned GetProperties();

  // add an option group containing all options related to the skyrmion wave functions
  //
  // manager = pointer to the option manager
  static void AddSkyrmionOptionGroup(OptionManager &manager, QHEWaveFunctionManager *wfManager);

};


#endif

#endif

