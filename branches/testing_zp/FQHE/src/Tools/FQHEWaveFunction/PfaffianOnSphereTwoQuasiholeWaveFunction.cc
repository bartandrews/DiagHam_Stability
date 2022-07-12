////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of Pfaffian wave function with two quasiholes on sphere         //
//                                                                            //
//                        last modification : 18/07/2006                      //
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
#include "Tools/FQHEWaveFunction/PfaffianOnSphereTwoQuasiholeWaveFunction.h"
#include "Vector/RealVector.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"


// constructor
//
// nbrParticles = number of particles
// theta1 = position of the first quasihole (spherical coordinates, theta angle)
// phi1 = position of the first quasihole (spherical coordinates, phi angle)
// theta2 = position of the second quasihole (spherical coordinates, theta angle)
// phi2 = position of the second quasihole (spherical coordinates, phi angle)
// fermions = flag indicating whether to calculate bosonic or fermionic pfaffian

PfaffianOnSphereTwoQuasiholeWaveFunction::PfaffianOnSphereTwoQuasiholeWaveFunction(int nbrParticles, double theta1, double phi1, double theta2, double phi2, bool fermions)
{
  this->NbrParticles = nbrParticles;
  this->Theta1=theta1;
  this->Phi1=phi1;
  this->Theta2=theta2;
  this->Phi2=phi2;
  
  this->U1.Re=cos(0.5*phi1);
  this->U1.Im=-sin(0.5*phi1);
  this->U1*=cos(0.5*theta1);

  this->V1.Re=cos(0.5*phi1);
  this->V1.Im=sin(0.5*phi1);
  this->V1*=sin(0.5*theta1);
  
  this->U2.Re=cos(0.5*phi2);
  this->U2.Im=-sin(0.5*phi2);
  this->U2*=cos(0.5*theta2);

  this->V2.Re=cos(0.5*phi2);
  this->V2.Im=sin(0.5*phi2);
  this->V2*=sin(0.5*theta2);

  this->TmpPfaffian = ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->FermionFlag = fermions;
}



// copy constructor
//
// function = reference on the wave function to copy

PfaffianOnSphereTwoQuasiholeWaveFunction::PfaffianOnSphereTwoQuasiholeWaveFunction(const PfaffianOnSphereTwoQuasiholeWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->Theta1=function.Theta1;
  this->Phi1=function.Phi1;
  this->Theta2=function.Theta2;
  this->Phi2=function.Phi2;
  this->U1=function.U1;
  this->V1=function.V1;
  this->U2=function.U2;
  this->V2=function.V2;
  this->TmpPfaffian = ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->FermionFlag=function.FermionFlag;
}

// destructor
//

PfaffianOnSphereTwoQuasiholeWaveFunction::~PfaffianOnSphereTwoQuasiholeWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PfaffianOnSphereTwoQuasiholeWaveFunction::Clone ()
{
  return new PfaffianOnSphereTwoQuasiholeWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PfaffianOnSphereTwoQuasiholeWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp, T1, T2, T3, T4;
  Complex WaveFunction(1.0);
  double Theta;
  double Phi;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      Theta = x[i << 1];
      Phi = x[1 + (i << 1)];
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = sin(0.5 * (x[j << 1] - Theta)) * cos(0.5 * (Phi - x[1 + (j << 1)]));
	  Tmp.Im = sin(0.5 * (Theta + x[j << 1])) * sin(0.5 * (Phi - x[1 + (j << 1)]));
	  WaveFunction *= Tmp;
	  if (FermionFlag)
	    WaveFunction *= Tmp;
	  T1.Re = sin(0.5 * (x[i << 1] - Theta1)) * cos(0.5 * (Phi1 - x[1 + (i << 1)]));
	  T1.Im = sin(0.5 * (Theta1 + x[i << 1])) * sin(0.5 * (Phi1 - x[1 + (i << 1)]));

	  T2.Re = sin(0.5 * (x[j << 1] - Theta2)) * cos(0.5 * (Phi2 - x[1 + (j << 1)]));
	  T2.Im = sin(0.5 * (Theta1 + x[j << 1])) * sin(0.5 * (Phi2 - x[1 + (j << 1)]));

	  T3.Re = sin(0.5 * (x[i << 1] - Theta2)) * cos(0.5 * (Phi2 - x[1 + (i << 1)]));
	  T3.Im = sin(0.5 * (Theta2 + x[i << 1])) * sin(0.5 * (Phi2 - x[1 + (i << 1)]));

	  T4.Re = sin(0.5 * (x[j << 1] - Theta1)) * cos(0.5 * (Phi1 - x[1 + (j << 1)]));
	  T4.Im = sin(0.5 * (Theta1 + x[j << 1])) * sin(0.5 * (Phi1 - x[1 + (j << 1)]));
	  
	  Tmp = ( T1*T2-T3*T4 ) / Tmp;
	  this->TmpPfaffian.SetMatrixElement (i , j, Tmp);
	}
    }
  WaveFunction *= TmpPfaffian.Pfaffian();
  return WaveFunction;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex PfaffianOnSphereTwoQuasiholeWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  Complex Tmp;
  Complex TmpU;
  Complex TmpV;
  Complex TmpHole1;
  Complex TmpHole2;
  Complex WaveFunction(1.0);
  double Factor = M_PI * 0.5;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      TmpU = uv[2*i];
      TmpV = uv[2*i+1];
      TmpHole1 = (TmpU * this->V1 - TmpV * this->U1);
      TmpHole2 = (TmpU * this->V2 - TmpV * this->U2);
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{	  
	  Tmp = Factor * (TmpU * uv[2*j+1] - TmpV * uv[2*j] );
	  WaveFunction *= Tmp;
	  Tmp = ((TmpHole1 * (uv[2*j] * this->V2 - uv[2*j+1] * this->U2)) + ((uv[2*j] * this->V1 - uv[2*j+1] * this->U1) * TmpHole2)) / Tmp;
	  this->TmpPfaffian.SetMatrixElement (i , j, Tmp);
	}
    }
  if (this->FermionFlag == true)
    WaveFunction *= WaveFunction;
  WaveFunction *= this->TmpPfaffian.Pfaffian();
  return WaveFunction;
}
