////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  abstract class of particle in QHE systems                 //
//                                                                            //
//                        last modification : 10/10/2004                      //
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


#ifndef ABSTRACTQHEPARTICLE_H
#define ABSTRACTQHEPARTICLE_H


#include "config.h"
#include "MathTools/Complex.h"
#include "HilbertSpace/AbstractHilbertSpace.h"


class RealVector;
class ComplexVector;
class AbstractFunctionBasis;


class AbstractQHEParticle :  public AbstractHilbertSpace
{

 public:

  enum 
    {
      BosonicStatistic = 0x1,
      FermionicStatistic = 0x2,
    };
  
  // virtual destructor
  //
  virtual ~AbstractQHEParticle ();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic() = 0;

  // forge an eigenstate from a description given by a file
  //
  // filename = name of the file that contains the state description
  // state = reference on the vector where the state has to be stored
  // return value = true if no error occured
  virtual bool ForgeEigenstate(char* filename, RealVector& state);

  // forge an eigenstate from a description given by a file
  //
  // filename = name of the file that contains the state description
  // state = reference on the vector where the state has to be stored
  // return value = true if no error occured
  virtual bool ForgeEigenstate(char* filename, ComplexVector& state);

  // get the variance of the state
  // index = index of state to consider
  virtual int StateVariance (int index);

  // evaluate wave function in real space using a given basis
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis);

  // evaluate wave functions in real space using a given basis
  //
  // states = array of vector corresponding to the state in the Fock basis
  // nbrStates = number of states in the states array
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // waveFuntions = array where the  wave function values at the given location will be stored
  virtual void EvaluateWaveFunctions (RealVector* states, int nbrStates, RealVector& position, AbstractFunctionBasis& basis, Complex* waveFuntions);

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
  
  // evaluate wave functions in real space using a given basis and only for agiven range of components
  //
  // states = array of vector corresponding to the state in the Fock basis
  // nbrStates = number of states in the states array
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // waveFuntions = array where the  wave function values at the given location will be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  virtual void EvaluateWaveFunctions (RealVector* states, int nbrStates, RealVector& position, AbstractFunctionBasis& basis,
				      Complex* waveFuntions, int firstComponent, int nbrComponent);                                
  
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
  

 protected:

  // forge an eigenstate from a description given by a file
  //
  // filename = name of the file that contains the state description
  // nbrElements = number of elements per basis state 
  // realFlag = flag to indicates if coefficients are real numbers (true for real, false for complex)
  // nbrComponents = reference om the integer where the number of components will be stored
  // coefficients = reference on the array where coefficients (memory allocation is done by the method itself)
  // componentDescription = reference on the array where components description (memory allocation is done by the method itself)
  // return value = true if no error occured
  virtual bool LoadEigenstateDescrition(char* filename, int nbrElements, bool realFlag, int& nbrComponents, Complex*& coefficients, int**& componentDescription);

  
};

#endif


