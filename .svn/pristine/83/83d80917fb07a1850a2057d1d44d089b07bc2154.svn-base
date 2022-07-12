////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                       Copyright (C) 2009 Gunnar Möller                     //
//                                                                            //
//                                                                            //
//           class for calculation of a Gutzwiller state on the lattice       //
//                                                                            //
//                        last modification : 27/10/2009                      //
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


#include "GutzwillerOnLatticeWaveFunction.h"

#include "Hamiltonian/AbstractQHEOnLatticeHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "Tools/NewUnconstrainedOptimizsation.h"
#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;

// switch for controlling mode of operation for optimizer (keep first parameter constant, or not)
// #define HAVE_CONSTANT__1ST_PARAMETER
// switch for number of parameters controlling empty sites (1 parameter globally, or 1 parameter per site)
// #define SINGLE_EMPTY_PARAMETER

// constructor
// nbrParticles = particles in the condensate (should match space)
// hardCore = flag indicating whether double occupations may occur or whether hardcore particles are present
// space = target-space of many-body state
// variationalParameters = initial set of trial parameters
// symmetryType = assume some symmetry among parameters 0 = none, 1 = 2xsublattice density
GutzwillerOnLatticeWaveFunction::GutzwillerOnLatticeWaveFunction(int nbrParticles, bool hardCore, ParticleOnLattice *space, RealVector *variationalParameters, int symmetryType, AbstractRandomNumberGenerator *randomGenerator)
{
  this->NbrParticles = nbrParticles;
  this->Space = space;
  this->NbrSites = Space->GetNbrSites();
  this->HardCoreFlag = hardCore;
#ifdef SINGLE_EMPTY_PARAMETER
  this->NbrEmptyParameters = 1;
#else
  this->NbrEmptyParameters = this->NbrSites;
#endif
  this->SymmetryType = symmetryType;
  if (SymmetryType!=0)
    {
      if (SymmetryType==1)
	{
	  
	}
    }
  else
    if (this->HardCoreFlag)
      {
	this->MaxOccupation = 1;
	this->NbrVariationalParameters = 2*this->NbrSites+this->NbrEmptyParameters;
      }
    else
      {
	this->MaxOccupation = NbrParticles;
	this->NbrVariationalParameters = 2*this->NbrParticles*this->NbrSites+this->NbrEmptyParameters;
      }

  if (variationalParameters!=NULL)
    {
      this->VariationalParameters = RealVector(*variationalParameters,true);
      if (variationalParameters->GetVectorDimension()!=NbrVariationalParameters)
	{
	  cout << "Attention, inconsistent number of variational parameters"<<endl;
	  this->VariationalParameters.Resize(NbrVariationalParameters);
	}
    }
  else
    this->VariationalParameters = RealVector(NbrVariationalParameters);
  if (randomGenerator!=NULL)
    {
      this->RandomNumbers = randomGenerator;
      ExternalGenerator=true;
    }
  else
    {
      timeval RandomTime;
      gettimeofday (&(RandomTime), 0);
      this->RandomNumbers = new NumRecRandomGenerator(RandomTime.tv_sec);
      ExternalGenerator=false;
    }
    
  this->Dim = Space->GetHilbertSpaceDimension();
  this->Hamiltonian=NULL;
  this->Architecture=NULL;
}

// destructor
GutzwillerOnLatticeWaveFunction::~GutzwillerOnLatticeWaveFunction()
{
  delete this->RandomNumbers;
}

// get Many-Body state
// return = resultingState
ComplexVector & GutzwillerOnLatticeWaveFunction::GetGutzwillerWaveFunction()
{
  this->TargetVector.Resize(Dim);
  this->TargetVector.ClearVector();

  // call main recursion
  this->Product(NbrSites-1, 0, 0x0ul, 1.0);
  if (TargetVector.Norm()==0.0)
    cout << "Attention, obtained state with zero norm"<<endl;
  this->TargetVector/=this->TargetVector.Norm();
  //   cout <<"Test norm: "<<TargetVector.Norm()<<endl;
  return this->TargetVector;
}


void GutzwillerOnLatticeWaveFunction::SetVariationalParameters(RealVector &variationalParameters)
{
  if (variationalParameters.GetVectorDimension()!=NbrVariationalParameters)
    {
      cout << "Attention, inconsistent number of variational parameters"<<endl;
      variationalParameters.Resize(NbrVariationalParameters);
    }
  this->VariationalParameters.Copy(variationalParameters);
}

// set parameters to a random initial distribution (random phase)
void GutzwillerOnLatticeWaveFunction::SetToRandomPhase()
{
  Complex TmpC;
  
  for (int i=0; i<NbrEmptyParameters; ++i)
    this->VariationalParameters[i]=1.0;
  for (int i=0; i<NbrSites; ++i)
    {
      TmpC=Polar((double)NbrParticles/NbrSites,2.0*M_PI*RandomNumbers->GetRealRandomNumber());
      this->VariationalParameters[NbrEmptyParameters+2*i]=TmpC.Re;
      this->VariationalParameters[NbrEmptyParameters+2*i+1]=TmpC.Im;
      for (int n=2; n<MaxOccupation; ++n)
	for (int i=0; i<NbrSites; ++i)
	  {
	    this->VariationalParameters[NbrEmptyParameters+2*(n-1)*NbrSites+2*i]=TmpC.Re/n;
	    this->VariationalParameters[NbrEmptyParameters+2*(n-1)*NbrSites+2*i+1]=TmpC.Im/n;
	  }
    }
}

// import parameters in format of FQHELatticeCondensateState
void GutzwillerOnLatticeWaveFunction::ImportCondensate(ComplexVector &condensateState, RealVector &emptyNorms)
{
  if (condensateState.GetVectorDimension()!=this->NbrSites)
    {
      cout << "Wrong number of sites in condensate"<<endl;
      exit(1);
    }
  bool HaveNorms=true;
  if (emptyNorms.Norm()<1e-6)
    HaveNorms=false;
#ifdef SINGLE_EMPTY_PARAMETER
  if (HaveNorms)
    {
      cout << "Attention, ignoring constant parameters for empty sites"<<endl;
    }
  VariationalParameters[0]=1.0;
#else
  if (HaveNorms)
    for (int i=0; i<NbrSites; ++i)
      VariationalParameters[i]=emptyNorms[i];
  else
    for (int i=0; i<NbrSites; ++i)
      VariationalParameters[i]=1.0;
#endif
  if (HardCoreFlag)
    {
      for (int i=0; i<NbrSites; ++i)
	{
	  this->VariationalParameters[NbrEmptyParameters+2*i]  =condensateState[NbrSites-1-i].Re;
	  this->VariationalParameters[NbrEmptyParameters+2*i+1]=condensateState[NbrSites-1-i].Im;
	}
    }
  else
    {
      for (int i=0; i<NbrSites; ++i)
	{
	  Complex TmpC=1.0;
	  for (int n=0; n<MaxOccupation; ++n)
	    {
	      TmpC*=condensateState[NbrSites-1-i];
	      VariationalParameters[NbrEmptyParameters+2*n*NbrSites+2*i]  =TmpC.Re;// /(n+1);
	      VariationalParameters[NbrEmptyParameters+2*n*NbrSites+2*i+1]=TmpC.Im;// /(n+1);
	    }
	}
    }
}

// define a Hamiltonian to enable immediate evaluation of the energy
void GutzwillerOnLatticeWaveFunction::SetHamiltonian(AbstractQHEOnLatticeHamiltonian *hamiltonian)
{
  this->Hamiltonian=hamiltonian;
}
  
// define an architecture to enable multi-processor operations
void GutzwillerOnLatticeWaveFunction::SetArchitecture(AbstractArchitecture *architecture)
{
  this->Architecture=architecture;
}

// get expectation value of the energy
double GutzwillerOnLatticeWaveFunction::GetEnergy()
{
  if (this->Hamiltonian==NULL)
    {
      cout << "Please define a Hamiltonian first" << endl;
      exit(1);
    }
  if (this->Architecture==NULL)
    {
      cout << "Please define an Architecture first" << endl;
      exit(1);
    }
  this->GetGutzwillerWaveFunction();
  ComplexVector TmpState(this->Space->GetHilbertSpaceDimension());
  VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->TargetVector), &TmpState);
  Operation.ApplyOperation(this->Architecture);
  return Real(TargetVector*TmpState);
}




// main recursion to calculate State \prod_i (\sum \chi_i + \psi_i a^\dagger_i) |state>
// nextQ = value quantum number in next operator to be applied
// nbrBosons = number of bosons already in state
// state = state to be acted upon
// prefactor = previous coefficients applied to state
// in last stage of recursion, writes to this->TargetVector
void GutzwillerOnLatticeWaveFunction::Product (int nextQ, int nbrBosons, unsigned long state, Complex prefactor)
{
  // cout << "Calling: Product ("<<nextQ<<", "<< nbrBosons<<", "<<state<<", "<< prefactor<<")"<<endl;
  int Index;
  unsigned long ResultingState;
  double AdFactor;
  if (nextQ>0)
    {
      int NbrPlaced = 0;
      double TotalAdFactor=1.0;
      ResultingState = state;
      while ((nbrBosons+NbrPlaced<this->NbrParticles)&&((!HardCoreFlag)||(NbrPlaced<1)))
	{
	  ResultingState = Space->Ad(ResultingState,nextQ,AdFactor);
	  ++NbrPlaced;
	  TotalAdFactor/=AdFactor;
	  Product(nextQ-1, nbrBosons+NbrPlaced, ResultingState, prefactor*TotalAdFactor*
		  Complex(VariationalParameters[NbrEmptyParameters+2*nextQ+  2*(NbrPlaced-1)*NbrSites],
			  VariationalParameters[NbrEmptyParameters+2*nextQ+1+2*(NbrPlaced-1)*NbrSites]));
	}
#ifdef SINGLE_EMPTY_PARAMETER
      Product(nextQ-1, nbrBosons, state, prefactor*VariationalParameters[0]);
#else
      Product(nextQ-1, nbrBosons, state, prefactor*VariationalParameters[nextQ]);
#endif
    }
  else
    {
      int NbrPlaced = 0;
      double TotalAdFactor=1.0;
      ResultingState = state;
      while ((nbrBosons+NbrPlaced<this->NbrParticles)&&((!HardCoreFlag)||(NbrPlaced<1)))
	{
	  ResultingState = Space->Ad(ResultingState,nextQ,AdFactor);
	  ++NbrPlaced;
	  TotalAdFactor/=AdFactor;
	}
      if ((Index=Space->CarefulFindStateIndex(ResultingState,-1))<Dim)
	{
	  if (NbrPlaced>0)
	    TargetVector[Index]+= prefactor*TotalAdFactor*Complex(VariationalParameters[NbrSites+2*nextQ+2*(NbrPlaced-1)*NbrSites],
								  VariationalParameters[NbrSites+2*nextQ+1+2*(NbrPlaced-1)*NbrSites]);
	  else
#ifdef SINGLE_EMPTY_PARAMETER
	    TargetVector[Index]+= prefactor*TotalAdFactor*VariationalParameters[0];
#else
	    TargetVector[Index]+= prefactor*TotalAdFactor*VariationalParameters[nextQ];
#endif
	}
    }
}


#ifdef HAVE_CONSTANT__1ST_PARAMETER

// target function for optimizer routine:
// version for 1st parameter fixed to constant
double GutzwillerOnLatticeWaveFunction::EvaluateEnergy(int nbrParameters, double *x)
{
  for (int i=1; i<this->NbrVariationalParameters; ++i)
    if (this->VariationalParameters[i]!=x[i])
      this->VariationalParameters[i]=x[i];
#ifdef SINGLE_EMPTY_PARAMETER
  if (this->VariationalParameters[0]==0.0)
    this->VariationalParameters[0]=1e-6;
#endif
  ++this->NbrEvaluations;
  if ((this->NbrEvaluations<20)
      || ((this->NbrEvaluations<1000)&&(this->NbrEvaluations%10==0))
      || (this->NbrEvaluations%100==0))
    cout << ".";
  cout.flush();
  //cout << "new parameters:" << endl << this->VariationalParameters;
  return this->GetEnergy();
}

// optimize wavefunction starting from present settings of VariationalParameters
// tolerance = final tolerance on the variational parameters
// maxIter = maximal number of function evaluations
//
double GutzwillerOnLatticeWaveFunction::Optimize(double tolerance, int maxIter)
{
  double InitialStepSize=1.0;
  int EffectiveNbrVariationalParameters = NbrVariationalParameters-1;
  cout << "Starting Optimization ";
  this->NbrEvaluations=0;
  int NbrPoints = 2 * EffectiveNbrVariationalParameters + 1, rnf;
  double Result;
  double *Work = new double[(NbrPoints+13)*(NbrPoints+EffectiveNbrVariationalParameters)
			    + 3*EffectiveNbrVariationalParameters*(EffectiveNbrVariationalParameters+3)/2 + 12];
  // passing parameter vector to optimizer as vector indexed from 1, not 0:
  double *x = &this->VariationalParameters[0];
  double (GutzwillerOnLatticeWaveFunction::*TargetFunction)(int, double*)=&GutzwillerOnLatticeWaveFunction::EvaluateEnergy;
  GutzwillerOnLatticeWaveFunction *TargetObject=this;
  Result = NewUOA::newuoa(EffectiveNbrVariationalParameters, NbrPoints, x, InitialStepSize,
			  tolerance, &rnf, maxIter, Work, TargetObject, TargetFunction);
  cout << endl << "total: "<<NbrEvaluations<< " evaluations"<<endl;
  delete [] Work;
  return Result;
  
}


#else

// target function for optimizer routine:
// version with all parameters varying (norm not enforced)
double GutzwillerOnLatticeWaveFunction::EvaluateEnergy(int nbrParameters, double *x)
{
  if (nbrParameters!=this->NbrVariationalParameters)
    {
      cout << "Unexpected error"<<endl;
    }
  for (int i=0; i<this->NbrVariationalParameters; ++i)
    this->VariationalParameters[i]=x[i];
#ifdef SINGLE_EMPTY_PARAMETER
  if (fabs(this->VariationalParameters[0])<1e-15)
    this->VariationalParameters[0]=1e-6;
#endif
  ++this->NbrEvaluations;
  if ((this->NbrEvaluations<20)
      || ((this->NbrEvaluations<1000)&&(this->NbrEvaluations%10==0))
      || (this->NbrEvaluations%100==0))
    cout << ".";
  cout.flush();
  return this->GetEnergy();
}

// optimize wavefunction starting from present settings of VariationalParameters
// tolerance = final tolerance on the variational parameters
// maxIter = maximal number of function evaluations
//
double GutzwillerOnLatticeWaveFunction::Optimize(double tolerance, int maxIter)
{
  double InitialStepSize=1.0;
  int EffectiveNbrVariationalParameters = NbrVariationalParameters;
  cout << "Starting Optimization ";
  this->NbrEvaluations=0;
  int NbrPoints = 2 * EffectiveNbrVariationalParameters + 1, rnf;
  double Result;
  double *Work = new double[(NbrPoints+13)*(NbrPoints+EffectiveNbrVariationalParameters)
			    + 3*EffectiveNbrVariationalParameters*(EffectiveNbrVariationalParameters+3)/2 + 12];
  // passing parameter vector to optimizer as vector indexed from 1, not 0:
  double *x = new double[EffectiveNbrVariationalParameters];
  for (int i=0; i<EffectiveNbrVariationalParameters; ++i)
    x[i]=this->VariationalParameters[i];
  double (GutzwillerOnLatticeWaveFunction::*TargetFunction)(int, double*)=&GutzwillerOnLatticeWaveFunction::EvaluateEnergy;
  GutzwillerOnLatticeWaveFunction *TargetObject=this;
  Result = NewUOA::newuoa(EffectiveNbrVariationalParameters, NbrPoints, x, InitialStepSize,
			  tolerance, &rnf, maxIter, Work, TargetObject, TargetFunction);
  cout << endl << "total: "<<NbrEvaluations<< " evaluations"<<endl;
  delete [] Work;
  return Result;
  
}


#endif


