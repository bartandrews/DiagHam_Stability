////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2008 Gunnar Moeller                   //
//                                                                            //
//                                                                            //
//           class implementing composite fermion state with partially        //
//             filled highest CF shell for a wave function on sphere          //
//                                                                            //
//                        last modification : 16/01/2008                      //
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


#include "HundRuleBilayerSinglet.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
HundRuleBilayerSinglet::HundRuleBilayerSinglet()
{
  this->NbrParticlesPerLayer = 0;
}

// standard constructor
HundRuleBilayerSinglet::HundRuleBilayerSinglet(int nbrParticlesPerLayer, int nbrEffectiveFlux, int jastrowP )
{
  this->NbrParticlesPerLayer = nbrParticlesPerLayer;
  this->Part.Resize(2*NbrParticlesPerLayer);
  this->PartC.Resize(2*NbrParticlesPerLayer);  
  this->CFStates=new HundRuleCFStates(nbrParticlesPerLayer, nbrEffectiveFlux, jastrowP);
  this->LPerLayer=CFStates->GetTotalL();
  this->ResultsLayer1 = new Complex[LPerLayer+1];
  this->ResultsLayer2 = new Complex[LPerLayer+1];
  ClebschGordanCoefficients TmpCoefficients(LPerLayer,LPerLayer);
  this->NbrCouplings = 0;
  for (int m=-LPerLayer; m<=LPerLayer; m+=2)
    if (TmpCoefficients.GetCoefficient (m, -m, 0)!=0.0)
      ++this->NbrCouplings;
  this->MPositions = new int[this->NbrCouplings];
  this->Couplings = new double[this->NbrCouplings];
  int i=0;
  for (int m=-LPerLayer; m<=LPerLayer; m+=2)
    if (TmpCoefficients.GetCoefficient (m, -m, 0)!=0.0)
      {
	this->MPositions[i] = (m+LPerLayer)/2;
	this->Couplings[i] = TmpCoefficients.GetCoefficient (m, -m, 0);
	++i;
      }  
  this->Flag.Initialize();
}


// copy constructor
HundRuleBilayerSinglet::HundRuleBilayerSinglet(HundRuleBilayerSinglet &toCopy)
{
  this->NbrParticlesPerLayer = toCopy.NbrParticlesPerLayer;
  this->Part.Resize(2*NbrParticlesPerLayer);
  this->PartC.Resize(2*NbrParticlesPerLayer);  
  this->CFStates = toCopy.CFStates;
  this->LPerLayer = toCopy.LPerLayer;
  this->ResultsLayer1 = toCopy.ResultsLayer1;
  this->ResultsLayer2 = toCopy.ResultsLayer2;
  this->NbrCouplings = toCopy.NbrCouplings;
  this->MPositions = toCopy.MPositions;
  this->Couplings = toCopy.Couplings;
  this->Flag = toCopy.Flag;
}


// destructor
HundRuleBilayerSinglet::~HundRuleBilayerSinglet()
{
  if ((this->NbrParticlesPerLayer!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete CFStates;
      delete [] MPositions;
      delete [] Couplings;
      delete [] ResultsLayer1;
      delete [] ResultsLayer2;
    }
}

// assignment operator
HundRuleBilayerSinglet& HundRuleBilayerSinglet::operator = (HundRuleBilayerSinglet &toCopy)
{
  if ((this->NbrParticlesPerLayer!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete CFStates;
      delete [] MPositions;
      delete [] Couplings;
      delete [] ResultsLayer1;
      delete [] ResultsLayer2;
    }
  this->NbrParticlesPerLayer = toCopy.NbrParticlesPerLayer;
  this->Part.Resize(2*NbrParticlesPerLayer);
  this->PartC.Resize(2*NbrParticlesPerLayer);  
  this->CFStates = toCopy.CFStates;
  this->LPerLayer = toCopy.LPerLayer;
  this->ResultsLayer1 = toCopy.ResultsLayer1;
  this->ResultsLayer2 = toCopy.ResultsLayer2;
  this->NbrCouplings = toCopy.NbrCouplings;
  this->MPositions = toCopy.MPositions;
  this->Couplings = toCopy.Couplings;
  this->Flag = toCopy.Flag;
  return *this;
}

// clone function 
//
// return value = clone of the function 
Abstract1DComplexFunction* HundRuleBilayerSinglet::Clone ()
{
  return (Abstract1DComplexFunction*)( new HundRuleBilayerSinglet(*this));
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  
Complex HundRuleBilayerSinglet::operator() (RealVector& x)
{
  this->Part=x.Extract(0, 2*this->NbrParticlesPerLayer-1);
  this->CFStates->GetValues(Part, ResultsLayer1);

  this->Part=x.Extract(2*this->NbrParticlesPerLayer, 4*this->NbrParticlesPerLayer-1);
  this->CFStates->GetValues(Part, ResultsLayer2);

  Complex Result=0.0;
  Complex Tmp;

  for (int i=0; i<NbrCouplings; ++i)
    {
      Tmp = ResultsLayer1[this->MPositions[i]] * ResultsLayer2[this->LPerLayer-this->MPositions[i]];
      Tmp*=this->Couplings[i];
      Result += Tmp;
    }
  
  return Result;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex HundRuleBilayerSinglet::CalculateFromSpinorVariables(ComplexVector& uv)
{
  this->PartC=uv.Extract(0, 2*this->NbrParticlesPerLayer-1);
  this->CFStates->GetValuesFromSpinorVariables(PartC, ResultsLayer1);

  this->PartC=uv.Extract(2*this->NbrParticlesPerLayer, 4*this->NbrParticlesPerLayer-1);
  this->CFStates->GetValuesFromSpinorVariables(PartC, ResultsLayer2);

  Complex Result=0.0;
  Complex Tmp;

  for (int i=0; i<NbrCouplings; ++i)
    {
      Tmp = ResultsLayer1[this->MPositions[i]] * ResultsLayer2[this->LPerLayer-this->MPositions[i]];
      Tmp*=this->Couplings[i];
      Result += Tmp;
    }
  
  return Result;
}  



void HundRuleBilayerSinglet::AdaptNorm(RealVector& x)
{
  this->Part=x.Extract(0, 2*this->NbrParticlesPerLayer-1);
  this->CFStates->AdaptNorm(Part);
}


// utility function to set the right dynamic interval for Monte-Carlo
void HundRuleBilayerSinglet::AdaptAverageMCNorm(int thermalize, int average)
{
  this->CFStates->AdaptAverageMCNorm(thermalize, average);
}


