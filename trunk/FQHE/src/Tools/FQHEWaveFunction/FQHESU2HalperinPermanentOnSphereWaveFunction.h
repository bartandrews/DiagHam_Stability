////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of SU(2) Halperin wave function on sphere             //
//                            times the permament                             //
//                                                                            //
//                        last modification : 09/04/2008                      //
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


#ifndef FQHESU2HALPERINPERMANENTONSPHEREWAVEFUNCTION_H
#define FQHESU2HALPERINPERMANENTONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "Matrix/ComplexMatrix.h"


class FQHESU2HalperinPermanentOnSphereWaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // number of particles with spin up
  int NbrSpinUpParticles;
  // number of particles with spin down
  int NbrSpinDownParticles;
  // total number of particles
  int TotalNbrParticles;

  // m1 index ( (m1, m2, n) Halperin wave function)
  int M1Index;
  // m2 index ( (m1, m2, n) Halperin wave function)
  int M2Index;
  // n index ( (m1, m2, n) Halperin wave function)
  int NIndex;

  // temporary arrays used during wave function evaluation
  Complex* UCoordinates;
  Complex* VCoordinates;

  // indicate if the elements of  GeneralizedPermanentMatrix are the invert of the correlations
  bool InvertFlag;

  // temporary matrix that are used to compute the permanent
  ComplexMatrix Permanent12;

 public:

  // constructor
  //
  // nbrSpinUpParticles = number of particles with spin up
  // nbrSpinDownParticles = number of particles with spin down
  // m1Index = m1 index ( (m1, m2, n) Halperin wave function)
  // m2Index = m2 index ( (m1, m2, n) Halperin wave function)
  // nIndex = n index ( (m1, m2, n) Halperin wave function)
  // invertFlag = if true, use the invert of the matrix elements for the permanent
  FQHESU2HalperinPermanentOnSphereWaveFunction(int nbrSpinUpParticles, int nbrSpinDownParticles, int m1Index, int m2Index, int nIndex,
					       bool invertFlag = false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  FQHESU2HalperinPermanentOnSphereWaveFunction(const FQHESU2HalperinPermanentOnSphereWaveFunction& function);

  // destructor
  //
   ~FQHESU2HalperinPermanentOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  // return value = function value at (uv)
  Complex CalculateFromSpinorVariables(ComplexVector& uv);

};

#endif
