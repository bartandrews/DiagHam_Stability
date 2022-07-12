////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of FQHE Jack generator parallelization operation         //
//                                                                            //
//                        last modification : 02/02/2010                      //
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


#ifndef FQHESPHEREJACKGENERATORSUMRATIONALPOLYNOMIALOPERATION_H
#define FQHESPHEREJACKGENERATORSUMRATIONALPOLYNOMIALOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Polynomial/LongRationalPolynomial.h"
#include "GeneralTools/GarbageFlag.h"


class FQHESphereJackGeneratorSumRationalPolynomialOperation: public AbstractArchitectureOperation
{

 protected:

  // array of polynomials describing the rational polynomials 
  LongRationalPolynomial* Numerators;
  LongRationalPolynomial* Denominators;

  //  indices to the rational polynomials to consider in the polynomial arrays
  long* ConnectedIndices;
  // linear coefficients for the rational polynomial sum
  long* ConnectedCoefficients;
  // local number of rational polynomials to handle
  int LocalNbrRationalPolynomials;

  // local storage of polynomials
  LongRationalPolynomial** LocalNumerators;
  LongRationalPolynomial** LocalDenominators;
  // garbage collector flag for the local arrays
  GarbageFlag Flag;
 
  // local shift for global array indices
  int LocalShift;
  // local index to access local arrays
  int LocalIndex;
  // local step when summing up rational polynomials
  int LocalStep;
  // maximum number of jobs to handle
  int LocalNbrJobs;

  // a temporary rational used to evaluate polynomials
  LongRational TmpRational;
  
  // indicates that this is the first pass in the rational polynmial summation
  bool FirstPassFlag;

 public:
  
  // constructor 
  //
  // space = pointer to the Hilbert space to use
  // invAlpha = inverse of the Jack polynomial alpha coefficient
  // rootPartition = root partition (in fermionic binary representation)
  // indexArray = array where state indices are stored
  // stateArray = array use to store computed state description
  // componentArray = array where computed component numerical factors are stored
  // rhoArray = rho factor associated to each state
  // nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
  // fermionicFlag = true if we are dealing with fermions
  // symmetricFlag = true if the state is Lz<->-Lz symmetric
  FQHESphereJackGeneratorSumRationalPolynomialOperation(long index, LongRationalPolynomial* numerators, LongRationalPolynomial* denominators, long* connectedIndices, long* connectedCoefficients, int nbrConnected);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereJackGeneratorSumRationalPolynomialOperation(const FQHESphereJackGeneratorSumRationalPolynomialOperation& operation);
  
  // destructor
  //
  ~FQHESphereJackGeneratorSumRationalPolynomialOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
};

#endif
