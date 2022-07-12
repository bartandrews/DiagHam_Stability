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


#ifndef POLARIZEDPRODUCTWAVEFUNCTION_H
#define POLARIZEDPRODUCTWAVEFUNCTION_H


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


class PolarizedProductWavefunction: public Abstract1DComplexTrialFunction
{

 protected:

  // Pointer to architecture
  AbstractArchitecture* Architecture;

  // Number of particles;
  int NbrParticles;

  // Nbr of flux in first factor
  int FirstLzMax;

  // Angular moment(a) in first factor
  int FirstLz;

  // statistics of first state
  bool FirstBosonicStatistic;

  // data for second hilbert-space
  int SecondLzMax;
  // angular momentum of second state(s)
  int SecondLz;
  
  // total angular momentum projection
  int TotalLz;

  // statistics of second state
  bool SecondBosonicStatistic;
  
  // arrays of temporary results to be used internally in operator()
  // for first state
  Complex LastValueFirst;
  // for second state
  Complex LastValueSecond;
  
  // vector that describes first state components in PolarizedSpace basis
  RealVector *FirstState;

  // vector that describes second state components in BosonicSpace basis
  RealVector *SecondState;
  
  // Hilbert space associated to the first and second state
  ParticleOnSphere* FirstSpace;
  ParticleOnSphere* SecondSpace;
  
  // one body real space basis to use for polarized state
  AbstractFunctionBasis* FirstOneBodyBasis;
  AbstractFunctionBasis* SecondOneBodyBasis;
  
  // flag indicating whether we are using an exact wavefunction as a second factor
  bool UseExactSecond;

  // Analytic polarized wavefunction
  Abstract1DComplexFunction* AnalyticSecondWaveFunction;

  // type cast of function pointer, if suitable
  Abstract1DComplexTrialFunction* AnalyticSecondTrialWaveFunction;

  // minimum dimension of bosonic space before parallelisation is used
  int MinParallel;

  // power of optional jastrow factor
  int JastrowPower;

  // velue of Jastrow factor
  Complex Jastrow;

  // temporary array used to store u spinor coordinates
  Complex* SpinorUCoordinates;
  // temporary array used to store v spinor coordinates
  Complex* SpinorVCoordinates;

  // Garbage flag
  GarbageFlag Flag;

 public:

  // constructor
  //
  // create object to be initialized from system options
  //
  PolarizedProductWavefunction(AbstractArchitecture* architecture, OptionManager &manager, int nbrParticles,
			       int totalLzMax, int totalLz, QHEWaveFunctionManager *wfManager,
			       int basisType=ParticleOnSphereFunctionBasis::LeftHanded);


  // copy constructor
  //
  // function = reference on the wave function to copy
  PolarizedProductWavefunction(const PolarizedProductWavefunction& function);

  // destructor
  //
   ~PolarizedProductWavefunction();

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
  static void AddPolarizedProductStateOptionGroup(OptionManager &manager);

};


#endif

#endif

