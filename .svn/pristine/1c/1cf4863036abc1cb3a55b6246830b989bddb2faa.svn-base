////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of function basis for particle on checkerboard lattice      //
//                                                                            //
//                        last modification : 29/06/2011                      //
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
#include "FunctionBasis/ParticleOnCheckerboardLatticeFunctionBasis.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"



// constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// t2p = hoping amplitude between second next neareast neighbor sites
// mus = sublattice staggered chemical potential 

ParticleOnCheckerboardLatticeFunctionBasis::ParticleOnCheckerboardLatticeFunctionBasis(int nbrSiteX, int nbrSiteY, double t1, double t2, double t2p, double mus)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->InvNbrSiteX = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->InvNbrSiteY = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->SecondNextNNHoping = t2p;
  this->MuS = mus;
  this->OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  double GammaX = 0.0;
  double GammaY = 0.0;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = (kx * this->NbrSiteY) + ky;
	Complex B1 = 4.0 * this->NNHoping * Complex (cos (1.0 * M_PI * (((double) kx) + GammaX) / ((double) this->NbrSiteX)) * cos (1.0 * M_PI * (((double) ky) + GammaY) / ((double) this->NbrSiteY)) * cos(M_PI * 0.25), 
						     sin (1.0 * M_PI * (((double) kx) + GammaX) / ((double) this->NbrSiteX)) * sin (1.0 * M_PI * (((double) ky) + GammaY) / ((double) this->NbrSiteY)) * sin(M_PI * 0.25));
	double d1 = 4.0 * this->SecondNextNNHoping * cos (2.0 * M_PI * (((double) kx) + GammaX) / ((double) this->NbrSiteX)) * cos (2.0 * M_PI * (((double) ky) + GammaY) / ((double) this->NbrSiteY));
	double d3 =  this->MuS + (2.0 * this->NextNNHoping * (cos (2.0 * M_PI * (((double) kx) + GammaX) / ((double) this->NbrSiteX))
							      - cos (2.0 * M_PI * (((double) ky) + GammaY) / ((double) this->NbrSiteY))));
	HermitianMatrix TmpOneBobyHamiltonian(2, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, B1);
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	ComplexMatrix TmpMatrix(2, 2, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBobyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	this->OneBodyBasis[Index] = TmpMatrix;	
      }
}

// destructor
//

ParticleOnCheckerboardLatticeFunctionBasis::~ParticleOnCheckerboardLatticeFunctionBasis ()
{
  delete[] this->OneBodyBasis;
}

// get value of the index-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// index = linearized momentum index i.e. kx * Ny + ky
// result = reference on the value where the result has to be stored
// index = function index 

void ParticleOnCheckerboardLatticeFunctionBasis::GetFunctionValue(RealVector& value, Complex& result, int index)
{
  int Kx = index / this->NbrSiteY;
  int Ky = index % this->NbrSiteY;
  double Arg = (this->InvNbrSiteX * ((double) Kx) * value[0]) + (this->InvNbrSiteY * ((double) Ky) * value[1]);
  //  result = (this->NormalizationCoefficients[index] * Complex(cos(Arg), sin(Arg)));
  result = Phase(Arg);
}



