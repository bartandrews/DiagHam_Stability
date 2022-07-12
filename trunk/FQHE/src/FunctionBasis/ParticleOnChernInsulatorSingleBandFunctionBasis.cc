////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of function basis for particle on Chern insulator         //
//                        within single band approximation                    //
//                                                                            //
//                        last modification : 02/03/2011                      //
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
#include "FunctionBasis/ParticleOnChernInsulatorSingleBandFunctionBasis.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"



// constructor
//

ParticleOnChernInsulatorSingleBandFunctionBasis::ParticleOnChernInsulatorSingleBandFunctionBasis(int nbrSiteX, int nbrSiteY, double mass)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->BandParameter = mass;
  this->InvNbrSiteX = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->InvNbrSiteY = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NormalizationCoefficients = new Complex[this->NbrSiteX * this->NbrSiteY];
  double Factor = 1.0 / sqrt ((double) (this->NbrSiteX * this->NbrSiteY));
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
      {
	int Index = (kx1 * this->NbrSiteY) + ky1;
	double d1 = sin (this->InvNbrSiteX * ((double) kx1));
	double d2 = sin (this->InvNbrSiteY * ((double) ky1));
	double d3 = (this->BandParameter - cos (this->InvNbrSiteX * ((double) ky1))
		     - cos (this->InvNbrSiteY * ((double) kx1)));
	HermitianMatrix TmpOneBobyHamiltonian(2, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d3);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, Complex(d1, -d2));
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, -d3);
	ComplexMatrix TmpMatrix(2, 2, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBobyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	this->NormalizationCoefficients[Index] = Factor;//(TmpMatrix[0][0] + TmpMatrix[0][1]) * Factor;	
      }
}

// destructor
//

ParticleOnChernInsulatorSingleBandFunctionBasis::~ParticleOnChernInsulatorSingleBandFunctionBasis ()
{
  delete[] this->NormalizationCoefficients;
}

// get value of the index-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// index = linearized momentum index i.e. kx * Ny + ky
// result = reference on the value where the result has to be stored
// index = function index 

void ParticleOnChernInsulatorSingleBandFunctionBasis::GetFunctionValue(RealVector& value, Complex& result, int index)
{
  int Kx = index / this->NbrSiteY;
  int Ky = index % this->NbrSiteY;
  double Arg = (this->InvNbrSiteX * ((double) Kx) * value[0]) + (this->InvNbrSiteY * ((double) Ky) * value[1]);
  result = (this->NormalizationCoefficients[index] * Complex(cos(Arg), sin(Arg)));
}



