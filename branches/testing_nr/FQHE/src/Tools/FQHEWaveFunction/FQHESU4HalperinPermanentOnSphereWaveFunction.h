////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of SU(4) generalized Halperine wave function on sphere       //
//                        times the product of permaments                     //
//                                                                            //
//                        last modification : 17/04/2008                      //
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


#ifndef FQHESU4HALPERINPERMANENTONSPHEREWAVEFUNCTION_H
#define FQHESU4HALPERINPERMANENTONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"
#include "Matrix/ComplexMatrix.h"


class ConfigurationParser;


class FQHESU4HalperinPermanentOnSphereWaveFunction: public Abstract1DComplexFunctionOnSphere
{

 protected:

  // total number of particles
  int TotalNbrParticles;

  // number of type 1 particles
  int NbrN1;
  // number of type 2 particles
  int NbrN2;
  // number of type 3 particles
  int NbrN3;
  // number of type 4 particles
  int NbrN4;  

  // coefficient of the intra-component correlations in the 1 - 1 sector
  int M11;
  // coefficient of the intra-component correlations in the 2 - 2 sector
  int M22;
  // coefficient of the intra-component correlations in the 3 - 3 sector
  int M33;
  // coefficient of the intra-component correlations in the 4 - 4 sector
  int M44;
  // coefficient of the inter-component correlations in the 1 - 2 sector
  int M12;
  // coefficient of the inter-component correlations in the 1 - 3 sector
  int M13;
  // coefficient of the inter-component correlations in the 1 - 4 sector
  int M14;
  // coefficient of the inter-component correlations in the 2 - 3 sector
  int M23;
  // coefficient of the inter-component correlations in the 2 - 4 sector
  int M24;
  // coefficient of the inter-component correlations in the 3 - 4 sector
  int M34;

  // indicate if the elements of  GeneralizedPermanentMatrix are the invert of the correlations
  bool InvertFlag;

  // temporary matrices that are used to compute the permanents
  ComplexMatrix Permanent12;
  ComplexMatrix Permanent13;
  ComplexMatrix Permanent14;
  ComplexMatrix Permanent23;
  ComplexMatrix Permanent24;
  ComplexMatrix Permanent34;

  
 public:

  // constructor
  //
  // nbrN1 = number of type 1 particles
  // nbrN2 = number of type 2 particles
  // nbrN3 = number of type 3 particles
  // nbrN4 = number of type 4 particles
  // m11 = coefficient of the intra-component correlations in the 1 - 1 sector
  // m22 = coefficient of the intra-component correlations in the 2 - 2 sector
  // m33 = coefficient of the intra-component correlations in the 3 - 3 sector
  // m44 = coefficient of the intra-component correlations in the 4 - 4 sector
  // m12 = coefficient of the inter-component correlations in the 1 - 2 sector
  // m13 = coefficient of the inter-component correlations in the 1 - 3 sector
  // m14 = coefficient of the inter-component correlations in the 1 - 4 sector
  // m23 = coefficient of the inter-component correlations in the 2 - 3 sector
  // m24 = coefficient of the inter-component correlations in the 2 - 4 sector
  // m34 = coefficient of the inter-component correlations in the 3 - 4 sector
  // invertFlag = if true, use the invert of the matrix elements for the permanents
  FQHESU4HalperinPermanentOnSphereWaveFunction(int nbrN1, int nbrN2, int nbrN3, int nbrN4,
					       int m11, int m22, int m33, int m44,
					       int m12, int m13, int m14, int m23,
					       int m24, int m34,
					       bool invertFlag = false);

  // copy constructor
  //
  // function = reference on the wave function to copy
  FQHESU4HalperinPermanentOnSphereWaveFunction(const FQHESU4HalperinPermanentOnSphereWaveFunction& function);

  // constructor from configuration file
  //
  // configuration = reference on the configuration file parser
  // errorFlag = reference on the error flag that is set to false if an error occured while retriving the configuration
  // nbrParticles = reference on the total number of particles (computed from the configuration file datas)
  // lzMax = twice the maximum angular momentum for a particle (computed from the configuration file datas)
  FQHESU4HalperinPermanentOnSphereWaveFunction(ConfigurationParser& configuration, bool& errorFlag, int& nbrParticles, int& lzMax);

  // destructor
  //
   ~FQHESU4HalperinPermanentOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point (the first 2*N1 coordinates correspond to the position of the type 1 particles, 
  //                                     the following 2*N2 coordinates correspond to the position of the type 2 particles,
  //                                     the following 2*N3 coordinates correspond to the position of the type 3 particles,
  //                                     the last 2*N4 coordinates correspond to the position of the type 4 particles)
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  // return value = function value at (uv)
  Complex CalculateFromSpinorVariables(ComplexVector& uv);

};

#endif
