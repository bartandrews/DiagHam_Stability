////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of SkyrmionOnSphere state wave function                    //
//                                                                            //
//                        last modification : 20/04/2005                      //
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
#include "Tools/FQHEWaveFunction/SkyrmionOnSphereWaveFunction.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "Options/Options.h"

#include <iostream>
using std::cout;
using std::endl;

// constructor
//
// create object to be initialized from system options
//
SkyrmionOnSphereWaveFunction::SkyrmionOnSphereWaveFunction(AbstractArchitecture* architecture, OptionManager &manager, int nbrParticles, int totalLzMax, int totalLz, int totalSz, QHEWaveFunctionManager *waveFunctionManager, int basisType)
{
  this->Architecture=architecture;
  this->UseExact=true;
  this->NbrParticles=nbrParticles;
  int PolarizedParticles = manager.GetInteger("nbr-particles");
  this->PolarizedLzMax = manager.GetInteger("polarized-lzmax"); 
  this->PolarizedLz = manager.GetInteger("polarized-lz");
  this->MinParallel = manager.GetInteger("min-parallel");
  this->LastBosonicValue = 0.0;
  this->AnalyticPolarizedWaveFunction=NULL;
  this->AnalyticPolarizedTrialWaveFunction=NULL;

  this->Flag.Initialize();

  if (manager.GetString("polarized-state")==0)
    {
      if (waveFunctionManager!=NULL)
	{
	  UseExact=false;
	  PolarizedParticles=nbrParticles;
	  AnalyticPolarizedWaveFunction = waveFunctionManager->GetWaveFunction();	  
	  if (AnalyticPolarizedWaveFunction==0)
	    {
	      cout << "Failed to option valid analytic wavefunction. Please check your wavefunction options!"<<endl;
	      exit(-1);
	    }
	  if (AnalyticPolarizedWaveFunction->GetProperties() & Abstract1DComplexFunction::Trial)
	    AnalyticPolarizedTrialWaveFunction = (Abstract1DComplexTrialFunction*) AnalyticPolarizedWaveFunction;
	}
      else
	{
	  cout << "An exact state vector is required for Skyrmion::\"polarized state\""<<endl;
	  exit(-1);
	}
    }

  bool Statistics = true;
  
  if (UseExact)
    {
      PolarizedParticles = manager.GetInteger("nbr-particles"); 
      PolarizedLzMax = manager.GetInteger("polarized-lzmax"); 
      PolarizedLz = manager.GetInteger("polarized-lz");
      
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(manager.GetString("polarized-state"),
						       PolarizedParticles, PolarizedLzMax, PolarizedLz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters for polarized state " <<
	    manager.GetString("polarized-state") << endl;
	  exit(-1);
	}
      
      if (PolarizedParticles!=nbrParticles)
	{
	  cout << "Error: polarized and exact state have to have the same number of particles";
	  exit(-1);
	}
      
      if (PolarizedLzMax>totalLzMax)
	{
	  cout << "Error: polarized state has to be at a lower flux than the exact state";
	  exit(-1);
	}
      
      if (PolarizedState.ReadVector (manager.GetString("polarized-state")) == false)
	{
	  cout << "can't open vector file " << manager.GetString("polarized-state") << endl;
	  exit(-1);      
	}
    }
  else
    {
      if (PolarizedLzMax>totalLzMax)
	{
	  cout << "Error: polarized state has to be at a lower flux than the exact state";
	  exit(-1);
	}
    }


  int NbrBosons = manager.GetInteger("nbr-particles"); 
  this->BosonLzMax = 0;
  this->BosonLz = 0;
  this->BosonSz = 0;
  bool SzSymmetrizedBasis = false;
  bool SzMinusParity = false;
  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  Statistics = false;  

  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(manager.GetString("bosonic-state"), NbrBosons,
							   BosonLzMax, BosonLz, BosonSz,
							   SzSymmetrizedBasis, SzMinusParity, 
							   LzSymmetrizedBasis, LzMinusParity, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << manager.GetString("bosonic-state") << endl;
      exit(-1);
    }

  if (NbrBosons!=nbrParticles)
    {
      cout << "Error: bosonic and exact state have to have the same number of particles";
      exit(-1);
    }
  
  if (BosonLzMax==0)
    BosonLzMax = totalLzMax-PolarizedLzMax;
  else if (BosonLzMax!=totalLzMax-PolarizedLzMax)
    {
      cout << "Error: total flux has to match: BosonLzMax == TotalLzMax-PolarizedLzMax";
      exit(-1);
    }

  if (BosonLz==0)
    BosonLz = totalLz;
  else if (BosonLz!=totalLz)
    {
      cout << "Error: total angular momentum has to match: BosonLz == TotalLz";
      exit(-1);
    }

  if (BosonSz==0)
    BosonSz = totalSz;
  else if (BosonSz!=totalSz)
    {
      cout << "Error: total spin has to match: BosonSz == TotalSz";
      exit(-1);
    }

  if (BosonicState.ReadVector (manager.GetString("bosonic-state")) == false)
    {
      cout << "can't open vector file " << manager.GetString("bosonic-state") << endl;
      exit(-1);      
    }
  
  this->PolarizedSpace=NULL;
  if (UseExact)
    {
#ifdef __64_BITS__
      if (PolarizedLzMax <= 63)
#else
	if (PolarizedLzMax <= 31)
#endif
	  {	
	    PolarizedSpace = new FermionOnSphere(NbrParticles, PolarizedLz, PolarizedLzMax);
	  }
	else
#ifdef __128_BIT_LONGLONG__
	  if (PolarizedLzMax <= 126)
#else
	    if (PolarizedLzMax <= 62)
#endif
	      {	    
		PolarizedSpace = new FermionOnSphereLong(NbrParticles, PolarizedLz, PolarizedLzMax);
	      }
	    else
	      {
		cout << "States of this polarized Hilbert space cannot be represented in a single word." << endl;
		exit(-1);
	      }
    }
  
  this->BosonicSpace=new BosonOnSphereWithSpin(NbrParticles, BosonLz, BosonLzMax, BosonSz);
  
  this->OneBodyBasisPol = new ParticleOnSphereFunctionBasis(PolarizedLzMax, basisType);
  this->OneBodyBasisBos = new ParticleOnSphereFunctionBasis(BosonLzMax, basisType);
}

// copy constructor
//
// function = reference on the wave function to copy

SkyrmionOnSphereWaveFunction::SkyrmionOnSphereWaveFunction(const SkyrmionOnSphereWaveFunction& function)
{
  this->NbrParticles=function.NbrParticles;
  this->PolarizedLzMax=function.PolarizedLzMax;
  this->PolarizedLz=function.PolarizedLz;
  this->PolarizedState=function.PolarizedState;
  this->BosonicState=function.BosonicState;
  this->PolarizedSpace=function.PolarizedSpace;
  this->BosonicSpace=function.BosonicSpace;
  this->OneBodyBasisBos=function.OneBodyBasisBos;
  this->OneBodyBasisPol=function.OneBodyBasisPol;
  this->UseExact=function.UseExact;  
  this->AnalyticPolarizedWaveFunction=function.AnalyticPolarizedWaveFunction;
  this->Flag=function.Flag;
}

// destructor
//
SkyrmionOnSphereWaveFunction::~SkyrmionOnSphereWaveFunction()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (UseExact)
	delete this->PolarizedSpace;
      delete this->BosonicSpace;
      delete this->OneBodyBasisPol;
      delete this->OneBodyBasisBos;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* SkyrmionOnSphereWaveFunction::Clone ()
{
  return new SkyrmionOnSphereWaveFunction(*this);
}

void SkyrmionOnSphereWaveFunction::TestSymmetries(ParticleOnSphereCollection *particles)
{
  Complex ValueSpin, ValuePolarized, ValueSpin2, ValuePolarized2, ValueSpin3, ValuePolarized3;
  if (UseExact)
    {
      QHEParticleWaveFunctionOperation Operation(PolarizedSpace, &PolarizedState, &(particles->GetPositions()),
						 this->OneBodyBasisPol, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      ValuePolarized = Operation.GetScalar();
    }
  else
    {
      ValuePolarized = (*(this->AnalyticPolarizedWaveFunction))(particles->GetPositions());
    }
  if (BosonicSpace->GetHilbertSpaceDimension()>this->MinParallel)
    {
      QHEParticleWaveFunctionOperation Operation(BosonicSpace, &BosonicState, &(particles->GetPositions()),
						 this->OneBodyBasisBos, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      ValueSpin = Operation.GetScalar();
      }
  else ValueSpin = BosonicSpace->EvaluateWaveFunction(BosonicState, particles->GetPositions(), *(this->OneBodyBasisBos));
  
  particles->ToggleHalfHalf();

  // recalculate:
  if (UseExact)
    {
      QHEParticleWaveFunctionOperation Operation(PolarizedSpace, &PolarizedState, &(particles->GetPositions()),
						 this->OneBodyBasisPol, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      ValuePolarized2 = Operation.GetScalar();
    }
  else
    {
      ValuePolarized2 = (*(this->AnalyticPolarizedWaveFunction))(particles->GetPositions());
    }
  if (BosonicSpace->GetHilbertSpaceDimension()>this->MinParallel)
    {
      QHEParticleWaveFunctionOperation Operation2(BosonicSpace, &BosonicState, &(particles->GetPositions()),
						  this->OneBodyBasisBos, /* TimeCoherence */ -1);
      Operation2.ApplyOperation(this->Architecture);      
      ValueSpin2 = Operation2.GetScalar();
    }
  else ValueSpin2 = BosonicSpace->EvaluateWaveFunction(BosonicState, particles->GetPositions(), (*this->OneBodyBasisBos));
  
  // rotate all particles      
  particles->RotateAll(0.781723465, 2.13428571);

  // recalculate:
  if (UseExact)
    {
      QHEParticleWaveFunctionOperation Operation(PolarizedSpace, &PolarizedState, &(particles->GetPositions()),
						 this->OneBodyBasisPol, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      ValuePolarized3 = Operation.GetScalar();
    }
  else
    {
      ValuePolarized3 = (*(this->AnalyticPolarizedWaveFunction))(particles->GetPositions());
    }
  if (BosonicSpace->GetHilbertSpaceDimension()>this->MinParallel)
    {
      QHEParticleWaveFunctionOperation Operation3(BosonicSpace, &BosonicState, &(particles->GetPositions()),
						  this->OneBodyBasisBos, /* TimeCoherence */ -1);
      Operation3.ApplyOperation(this->Architecture);
      ValueSpin3 = Operation3.GetScalar();
    }
  else ValueSpin3 = BosonicSpace->EvaluateWaveFunction(BosonicState, particles->GetPositions(), *(this->OneBodyBasisBos));
  
  cout << "Pol Before exchange: "<< ValuePolarized << endl << "After exchange:  " << ValuePolarized2 << endl;
  cout << "Pol Parity: " << ValuePolarized/ValuePolarized2 << endl;
  cout << "Pol After rotation: " << ValuePolarized3 << " ratio: "<< Norm(ValuePolarized/ValuePolarized3) << endl;
  cout << "Bos Before exchange: "<< ValueSpin  << endl << "After exchange:  " << ValueSpin2 << endl;
  cout << "Bos Parity: " << ValueSpin / ValueSpin2 << endl;
  cout << "Bos After rotation: " << ValueSpin3 << " ratio: "<< Norm(ValueSpin/ValueSpin3) << endl;
}


// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex SkyrmionOnSphereWaveFunction::operator ()(RealVector& x)
{
  Complex ValuePolarized;
  if (UseExact)
    {
      QHEParticleWaveFunctionOperation Operation(PolarizedSpace, &PolarizedState, &x,
						 this->OneBodyBasisPol, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      ValuePolarized = Operation.GetScalar();
    }
  else
    {
      ValuePolarized = (*(this->AnalyticPolarizedWaveFunction))(x);
    }
  if (BosonicSpace->GetHilbertSpaceDimension()>this->MinParallel)
    {  
      QHEParticleWaveFunctionOperation Operation(BosonicSpace, &BosonicState, &x,
						 this->OneBodyBasisBos, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      this->LastBosonicValue = Operation.GetScalar();
    }
  else this->LastBosonicValue = BosonicSpace->EvaluateWaveFunction(BosonicState, x, *(this->OneBodyBasisBos));
  
  return ValuePolarized*this->LastBosonicValue;
}

// get a value of the wavefunction for the last set of coordinates, but with different variational parameters
// parameters =  alternative set of parameters
Complex SkyrmionOnSphereWaveFunction::GetForOtherParameters(double *parameters)
{
  if ((UseExact)|(AnalyticPolarizedTrialWaveFunction==NULL))
    {
      cout << "Cannot change any parameters in SkyrmionOnSphereWaveFunction when calculating from exact state"
	   << ", or using non variational wavefunction."<<endl;
      exit(1);
    }
  Complex Result=this->AnalyticPolarizedTrialWaveFunction->GetForOtherParameters(parameters);
  
  return Result*this->LastBosonicValue;
}

// do many evaluations of the function, storing the result in the vector results given in the call
// result = vector of leading dimension of the array coefficients for returning values
// x = positions to evaluate the wavefuntion in
// format for passing parameters as [nbrSet][nbrParameter],
void SkyrmionOnSphereWaveFunction::GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients)
{
  if ((UseExact)|(AnalyticPolarizedTrialWaveFunction==NULL))
    {
      cout << "Cannot change any parameters in SkyrmionOnSphereWaveFunction when calculating from exact state"
	   << ", or using non variational wavefunction."<<endl;
      exit(1);
    }
  this->AnalyticPolarizedTrialWaveFunction->GetForManyParameters(results,x,coefficients);

  Complex BosonicPart;
  if (BosonicSpace->GetHilbertSpaceDimension()>this->MinParallel)
    {  
      QHEParticleWaveFunctionOperation Operation(BosonicSpace, &BosonicState, &x,
						 this->OneBodyBasisBos, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      BosonicPart = Operation.GetScalar();
    }
  else BosonicPart = BosonicSpace->EvaluateWaveFunction(BosonicState, x, *(this->OneBodyBasisBos));
  for (int i=0; i<results.GetVectorDimension(); ++i)
    results[i]*=BosonicPart;
}

// access internal values of parameters
double *SkyrmionOnSphereWaveFunction::GetTrialParameters()
{
  if ((UseExact)|(AnalyticPolarizedTrialWaveFunction==NULL))
    {
      cout << "Cannot get any parameters in SkyrmionOnSphereWaveFunction when calculating from exact state"
	   << ", or using non variational wavefunction."<<endl;
      exit(1);
    }
  return this->AnalyticPolarizedTrialWaveFunction->GetTrialParameters();
}

// get number of parameters
int SkyrmionOnSphereWaveFunction::GetNbrParameters() 
{
  if ((UseExact)|(AnalyticPolarizedTrialWaveFunction==NULL))
    {
      cout << "No parameters available in SkyrmionOnSphereWaveFunction when calculating from exact state"
	   << ", or using non variational wavefunction."<<endl;
      exit(1);
    }
  return this->AnalyticPolarizedTrialWaveFunction->GetNbrParameters();
}

  
// set new values of the trial coefficients (keeping the initial number of parameters)
void SkyrmionOnSphereWaveFunction::SetTrialParameters(double * coefficients)
{
  if ((UseExact)|(AnalyticPolarizedTrialWaveFunction==NULL))
    {
      cout << "Cannot change any parameters in SkyrmionOnSphereWaveFunction when calculating from exact state"
	   << ", or using non variational wavefunction."<<endl;
      exit(1);
    }
  this->AnalyticPolarizedTrialWaveFunction->SetTrialParameters(coefficients);
}


// get function properties, and possible extensions of interface 
// 
unsigned SkyrmionOnSphereWaveFunction::GetProperties()
{
  if (UseExact)
    return Abstract1DComplexFunction::Basic;
  else
    return this->AnalyticPolarizedWaveFunction->GetProperties();
}




// add an option group containing all options related to the skyrmion wave functions
//
// manager = pointer to the option manager
void SkyrmionOnSphereWaveFunction::AddSkyrmionOptionGroup(OptionManager &manager, QHEWaveFunctionManager *wfManager)
{  
  OptionGroup* SkyrmionGroup = new OptionGroup ("skyrmion options");
  manager+=SkyrmionGroup;
  
  (*SkyrmionGroup) += new SingleStringOption  ('\n', "polarized-state", "file name of polarized fermionic reference wave function (if omitted using analytic function)",0);
  (*SkyrmionGroup) += new SingleIntegerOption  ('l', "polarized-lzmax", "total number of flux quanta (0 if it has to be guessed from input file name)", 0);  
  (*SkyrmionGroup) += new SingleIntegerOption  ('z', "polarized-lz", "twice the total lz value of the system (0 if it has to be guessed from input file name)", 0);
  (*SkyrmionGroup) += new SingleStringOption  ('\n', "bosonic-state", "file name of spinful bosonic part of wave function",0);
  (*SkyrmionGroup) += new SingleIntegerOption  ('\n', "min-parallel", "minimum dimension of bosonic space before parallel evaluation is applied", 5000);

  if (wfManager!=NULL)
    wfManager->AddOptionGroup(&manager);
}
