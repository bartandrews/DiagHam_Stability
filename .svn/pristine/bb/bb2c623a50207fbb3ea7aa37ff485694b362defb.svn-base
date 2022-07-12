////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of Moore Read state wave function on sphere             //
//                                                                            //
//                        last modification : 19/09/2004                      //
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
#include "Tools/FQHEWaveFunction/AdvancedReadRezayiOnSphereWaveFunction.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/GroupedPermutations.h"
#include "Vector/RealVector.h"

#include "GeneralTools/OrderedList.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;


// constructor
//
// nbrParticlesPerCluster = number of particles per cluster (=N/k)
// nbrClusters = number of clusters
// fermionicStatistics = flag indicating whether the pfaffian should be multiplied by a squared Jastrow Factor
AdvancedReadRezayiOnSphereWaveFunction::AdvancedReadRezayiOnSphereWaveFunction(int nbrParticlesPerCluster, int nbrClusters, bool fermionicStatistics)
{
  this->NbrParticles = 2*nbrParticlesPerCluster;
  this->ClusterSize = nbrParticlesPerCluster;
  this->NbrClusters = nbrClusters;
  this->FermionicStatistics = fermionicStatistics;
  this->BlockJk = new Complex[NbrClusters];
  this->SpinorUCoordinates = new Complex[this->NbrParticles];
  this->SpinorVCoordinates = new Complex[this->NbrParticles];
  this->JastrowFactorElements = new Complex*[this->NbrParticles];
  this->JastrowFactorSquares = new Complex*[this->NbrParticles];  
  for (int i=0; i<NbrParticles; ++i)
    {
      this->JastrowFactorElements[i] = new Complex[this->NbrParticles];
      this->JastrowFactorSquares[i] = new Complex[this->NbrParticles];
    }
  this->EvaluatePermutations();
  this->Flag.Initialize();
}

// copy constructor
//
// function = reference on the wave function to copy

AdvancedReadRezayiOnSphereWaveFunction::AdvancedReadRezayiOnSphereWaveFunction(const AdvancedReadRezayiOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->ClusterSize = function.ClusterSize;
  this->NbrClusters = function.NbrClusters;
  this->FermionicStatistics = function.FermionicStatistics;
  this->BlockJk = function.BlockJk;
  this->SpinorUCoordinates = function.SpinorUCoordinates;
  this->SpinorVCoordinates = function.SpinorVCoordinates;
  this->Permutations = function.Permutations;
  this->NbrPermutations = function.NbrPermutations;
  this->WeightOfPermutations = function.WeightOfPermutations;
  this->JastrowFactorElements = function.JastrowFactorElements;
  this->JastrowFactorSquares = function.JastrowFactorSquares;
  this->Flag = function.Flag;
}

// destructor
//

AdvancedReadRezayiOnSphereWaveFunction::~AdvancedReadRezayiOnSphereWaveFunction()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      for (unsigned i = 0; i < this->NbrPermutations; ++i)
	delete[] this->Permutations[i];
      delete[] this->Permutations;
      delete[] this->BlockJk;
      delete[] this->SpinorUCoordinates;
      delete[] this->SpinorVCoordinates;
      for (int i=0; i<NbrParticles; ++i)
	{
	  delete [] this->JastrowFactorSquares[i];
	  delete [] this->JastrowFactorElements[i];
	}
      delete [] this->JastrowFactorSquares;
      delete [] this->JastrowFactorElements;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* AdvancedReadRezayiOnSphereWaveFunction::Clone ()
{
  return new AdvancedReadRezayiOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex AdvancedReadRezayiOnSphereWaveFunction::operator ()(RealVector& x)
{  
  // CalculateSpinors
  double s,c;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
      this->SpinorUCoordinates[i].Im = this->SpinorUCoordinates[i].Re;
      this->SpinorUCoordinates[i].Re *= (c=cos(0.5 * x[1 + (i << 1)]));
      this->SpinorUCoordinates[i].Im *= -(s=sin(0.5 * x[1 + (i << 1)]));
      this->SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
      this->SpinorVCoordinates[i].Im = this->SpinorVCoordinates[i].Re;
      this->SpinorVCoordinates[i].Re *= c;
      this->SpinorVCoordinates[i].Im *= s;
      //cout << "U["<<i<<"]="<<SpinorUCoordinates[i]<<", "<< "V["<<i<<"]="<<SpinorVCoordinates[i]<<endl;
    }

  return this->ComplexEvaluations();
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex AdvancedReadRezayiOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{  
  // Import from spinors
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = uv.Re(2*i);
      this->SpinorUCoordinates[i].Im = uv.Im(2*i);
      this->SpinorVCoordinates[i].Re = uv.Re(2*i+1);
      this->SpinorVCoordinates[i].Im = uv.Im(2*i+1);
    }

  return this->ComplexEvaluations();
}  


// evaluate permutations required for the Moore-Read state evaluation
// using: two symmetric blocs, only permutations changing particles
//        between blocks are required

void AdvancedReadRezayiOnSphereWaveFunction::EvaluatePermutations()
{  
  GroupedPermutations Generator(this->NbrClusters, this->ClusterSize);
  FactorialCoefficient NbrP;  
  this->NbrPermutations = Generator.GetNbrPermutations();
  this->Permutations = new unsigned*[NbrPermutations];
  this->WeightOfPermutations = new double[NbrPermutations];
  unsigned long *Multiplicities = Generator.GetMultiplicities();
  SmallIntegerArray *TmpPermutations = Generator.GetPermutations();
  for (unsigned i=0; i<NbrPermutations; ++i)
    {
      this->Permutations[i] = new unsigned[NbrParticles];
      TmpPermutations[i].GetElements(this->Permutations[i]);
      cout << "Permutation["<<i<<"]= ["<<Permutations[i][0];
      for (int j=1; j<NbrParticles; ++j) cout << " " << Permutations[i][j];
      cout << "]"<<endl;
      NbrP.SetToOne();
      NbrP.FactorialDivide(NbrParticles);
      NbrP *= Multiplicities[i];
      this->WeightOfPermutations[i] = NbrP.GetNumericalValue();
    }

  return;
}

// perform complex part of calculations
// uses internal spinor coordinates as input
//
Complex AdvancedReadRezayiOnSphereWaveFunction::ComplexEvaluations()
{
  Complex J, Tmp;

  for (int i = 0; i < this->NbrParticles; ++i)
    for (int j = 0; j < i; ++j)
      {
	Tmp = ((this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]));
	JastrowFactorElements[i][j] = Tmp;
	JastrowFactorElements[j][i] = -Tmp;
      }

  unsigned *TmpP;
  Complex Value(0.0,0.0);
  // cout << "Evaluating function"<<endl;
  for (unsigned i=0; i<NbrPermutations; ++i)
    {
      TmpP=this->Permutations[i];
      J=1.0;
      for (int m=0; m<NbrClusters; ++m)
	{
	  BlockJk[m]=1.0;
	  for (int k=m*ClusterSize+1; k<(m+1)*ClusterSize; ++k)
	    for (int j=m*ClusterSize; j<k; ++j)
	      BlockJk[m]*=JastrowFactorElements[TmpP[k]][TmpP[j]];
	  J*=BlockJk[m];
	}
      Value += WeightOfPermutations[i]*J*J;
    }
  // cout << "Symmetric Part: "<<Value<<endl;
  if (this->FermionicStatistics)
    {
      for (int i = 1; i < this->NbrParticles; ++i)
	for (int j = 0; j < i; ++j)
	  Value *=  JastrowFactorElements[i][j];
    }
  //  cout << "Value ="<<Value<<endl;
  
  return Value;
}



/* code that might be useful later on for evaluation of overlaps (wrong for individual function values)

// evaluate permutations required for the Moore-Read state evaluation
// using: two symmetric blocs, only permutations changing particles
//        between blocks are required

void AdvancedReadRezayiOnSphereWaveFunction::EvaluatePermutations()
{
  double totalWeight=0.0;
  BinomialCoefficients bico(NbrParticles);
  this->NbrPermutations = this->ClusterSize+1;
  this->WeightOfPermutations = new double[NbrPermutations];
  this->Permutations = new unsigned*[NbrPermutations];
  for (unsigned i=0; i<NbrPermutations; ++i)
    this->Permutations[i] = new unsigned[NbrParticles];
  for (int i=0; i<NbrParticles; ++i)
    this->Permutations[0][i] = i;
  this->WeightOfPermutations[0] = 1.0/(double)bico(NbrParticles,ClusterSize);
  totalWeight+=WeightOfPermutations[0];
  for (int i=0; i<ClusterSize; ++i)
    {
      for (int j=0; j<NbrParticles; ++j)
	this->Permutations[i+1][j] = this->Permutations[i][j];
      this->Permutations[i+1][i] = i+ClusterSize;
      this->Permutations[i+1][i+ClusterSize] = i;
      this->WeightOfPermutations[i+1] = (double)bico(ClusterSize,i+1)*(double)bico(ClusterSize,i+1)/(double)bico(NbrParticles,ClusterSize);
      totalWeight+=WeightOfPermutations[i+1];
      cout << "Permutation "<<i+1<<": [ "<<Permutations[i+1][0];
      for (int k=1; k<NbrParticles; ++k) cout<<", "<<Permutations[i+1][k];
      cout << "] has weight "<<WeightOfPermutations[i+1]<<endl;
    }
  cout << "TotalWeight="<<totalWeight<<endl;
  return;
}
*/
