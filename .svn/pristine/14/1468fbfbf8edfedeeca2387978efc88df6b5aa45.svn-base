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


#include "config.h"
#include "HilbertSpace/AbstractQHEParticle.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Vector/RealVector.h"

#include <iostream>


using std::cout;
using std::endl;


// virtual destructor
//

AbstractQHEParticle::~AbstractQHEParticle ()
{
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured
bool AbstractQHEParticle::WriteHilbertSpace (char* fileName)
{
  cout << "Warning : WriteHilbertSpace not implemented" <<endl;
  return false;
}

// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex AbstractQHEParticle::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}

// evaluate wave functions in real space using a given basis
//
// states = array of vector corresponding to the state in the Fock basis
// nbrStates = number of states in the states array
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// waveFuntions = array where the  wave function values at the given location will be stored

void AbstractQHEParticle::EvaluateWaveFunctions (RealVector* states, int nbrStates, RealVector& position, AbstractFunctionBasis& basis, 
						 Complex* waveFuntions)
{
  this->EvaluateWaveFunctions(states, nbrStates, position, basis, waveFuntions, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// return value = wave function evaluated at the given location

Complex AbstractQHEParticle::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								    AbstractFunctionBasis& basis, int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, 
						     this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location
Complex AbstractQHEParticle::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
						int firstComponent, int nbrComponent)
{
  return Complex(0.0, 0.0);
}

// evaluate wave functions in real space using a given basis and only for agiven range of components
//
// states = array of vector corresponding to the state in the Fock basis
// nbrStates = number of states in the states array
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// waveFuntions = array where the  wave function values at the given location will be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void AbstractQHEParticle::EvaluateWaveFunctions (RealVector* states, int nbrStates, RealVector& position, AbstractFunctionBasis& basis,
						 Complex* waveFuntions, int firstComponent, int nbrComponent)
{
  for (int i = 0; i < nbrStates; ++i)
    waveFuntions[i] = this->EvaluateWaveFunction(states[i], position, basis, firstComponent, nbrComponent);
}
  

// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex AbstractQHEParticle::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								    AbstractFunctionBasis& basis, 
								    int nextCoordinates, int firstComponent, 
								    int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void AbstractQHEParticle::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}

                                    
// forge an eigenstate from a description given by a file
//
// filename = name of the file that contains the state description
// state = reference on the vector where the state has to be stored
// return value = true if no error occured

bool AbstractQHEParticle::ForgeEigenstate(char* filename, RealVector& state)
{
  return false;
}

// forge an eigenstate from a description given by a file
//
// filename = name of the file that contains the state description
// state = reference on the vector where the state has to be stored
// return value = true if no error occured

bool AbstractQHEParticle::ForgeEigenstate(char* filename, ComplexVector& state)
{
  return false;
}

// get the variance of the state
// index = index of state to consider
int AbstractQHEParticle::StateVariance (int index)
{
  cout << "Variance requested - needs to be defined in derived class"<<endl;
  return 0.0;
}


// forge an eigenstate from a description given by a file
//
// filename = name of the file that contains the state description
// nbrElements = number of elements per basis state 
// realFlag = flag to indicates if coefficients are real numbers (true for real, false for complex)
// nbrComponents = reference om the integer where the number of components will be stored
// coefficients = reference on the array where coefficients (memory allocation is done by the method itself)
// componentDescription = reference on the array where components description (memory allocation is done by the method itself)
// return value = true if no error occured

bool AbstractQHEParticle::LoadEigenstateDescrition(char* filename, int nbrElements, bool realFlag, int& nbrComponents, Complex*& coefficients, int**& componentDescription)
{
  ConfigurationParser Configuration;
  if (Configuration.Parse(filename) == false)
    {
      Configuration.DumpErrors(cout);
      return false;
    }

  double* TmpCoefficients;
  int TmpNbrCoefficients;
  if (Configuration.GetAsDoubleArray("coefficients", ' ', TmpCoefficients, TmpNbrCoefficients) == false)
    {
      Configuration.DumpErrors(cout);
      return false;
    }
  if (realFlag == true)
    {
      nbrComponents = TmpNbrCoefficients;
      coefficients = new Complex[nbrComponents];
      for (int i = 0; i < nbrComponents; ++i)
	coefficients[i] = TmpCoefficients[i];
    }
  else
    {
      if ((TmpNbrCoefficients & 1) != 0)
	{
	  cout << "wrong number of elements for complex definition" << endl;
	  delete[] TmpCoefficients;
	  return false;
	}
      nbrComponents = TmpNbrCoefficients >> 1;
      coefficients = new Complex[nbrComponents];
      for (int i = 0; i < TmpNbrCoefficients; ++i)
	{
	  coefficients[i >> 1].Re = TmpCoefficients[i];
	  ++i;
	  coefficients[i >> 1].Im = TmpCoefficients[i];
	}
    }
  delete[] TmpCoefficients;

  int* TmpElements;
  int TmpNbrElements;
  if (Configuration.GetAsIntegerArray("descriptions", ' ', TmpElements, TmpNbrElements) == false)
    {
      delete[] coefficients;
      Configuration.DumpErrors(cout);
      return false;
    }
  if (TmpNbrElements != (nbrElements * nbrComponents))
    {
      cout << "wrong number of elements for basis state definition" << endl;
      delete[] coefficients;
      delete[] TmpElements;
      return false;
    }
  int Pos = 0;
  componentDescription = new int*[nbrComponents];
  for (int i = 0; i < nbrComponents; ++i)
    {
      componentDescription[i] = new int[nbrElements];      
      for (int j = 0; j < nbrElements; ++j)
	{
	  componentDescription[i][j] = TmpElements[Pos];
	  ++Pos;
	}
    }
  delete[] TmpElements;
  return true;
}

  // filter a Hilbert space according to some exclusion rules
  //
  //return value = dimension of the filtered Hilbert space
  
  long AbstractQHEParticle::FilterHilbertSpace()
  {
    return this->HilbertSpaceDimension; 
  }