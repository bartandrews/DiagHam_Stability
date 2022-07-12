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
#include "Tools/FQHEWaveFunction/PolarizedProductWavefunction.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "Options/Options.h"

#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;

// constructor
//
// create object to be initialized from system options
//
PolarizedProductWavefunction::PolarizedProductWavefunction(AbstractArchitecture* architecture, OptionManager &manager, int nbrParticles, int totalLzMax, int totalLz, QHEWaveFunctionManager *waveFunctionManager, int basisType)
{
  this->Architecture=architecture;
  this->UseExactSecond=true;
  this->NbrParticles=nbrParticles;
  int FirstParticles = manager.GetInteger("nbr-particles");

  if (NbrParticles!=FirstParticles)
    {
      cout << "Inconsistent numbers of particles"<<endl;
      exit(-1);
    }
  this->FirstLzMax = manager.GetInteger("first-lzmax");  
  this->FirstLz = manager.GetInteger("first-lz");
  this->FirstBosonicStatistic = manager.GetBoolean("first-boson");
  
  this->SecondLzMax = manager.GetInteger("second-lzmax");  
  this->SecondLz = manager.GetInteger("second-lz");
  this->SecondBosonicStatistic = manager.GetBoolean("second-boson");

  
  this->MinParallel = manager.GetInteger("min-parallel");  
  this->AnalyticSecondWaveFunction=NULL;
  this->AnalyticSecondTrialWaveFunction=NULL;

  this->TotalLz=totalLz;
  
  this->Flag.Initialize();

  this->JastrowPower= manager.GetInteger("jastrow-power");
  cout << "Adding Jastrow factor (z_i-z_j)^"<<this->JastrowPower<<endl;
  if (JastrowPower!=0)
    {
      this->SpinorUCoordinates = new Complex[NbrParticles];
      this->SpinorVCoordinates = new Complex[NbrParticles];
    }

  bool Statistics = true;
    
  char *TmpState = manager.GetString("first-state");
  this->FirstState = new RealVector();
  FirstParticles=0;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(TmpState,FirstParticles, FirstLzMax, FirstLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters for first state " <<
	TmpState << endl;
      exit(-1);
    }
  if (FirstParticles!=NbrParticles)
    {
      cout << "Error: first and first state have to have the same number of particles"<<endl;
      exit(-1);
    }
      
  if (FirstLzMax>totalLzMax)
    {
      cout << "Error: first state has to be at a lower flux than the exact state"<<endl;
      exit(-1);
    }
      
  if (FirstState->ReadVector (TmpState) == false)
    {
      cout << "can't open vector file " << manager.GetString("first-state") << endl;
      exit(-1);      
    }
  else
    {
      cout << "Read first state "<<TmpState<<" as vector with angular momentum Lz="<<FirstLz<<endl;
    }
  if (Statistics==FirstBosonicStatistic) // inverse conventions, here!
    {
      cout << "Statistics of first vector file do not match your settings"<<endl;
      exit(-1);
    }


  if (manager.GetString("second-state")==NULL)
    {
      if (waveFunctionManager!=NULL)
	{
	  if (this->SecondLz!=0)
	    {
	      cout << "Analytic product wavefunction can only be generated for cases with Lz_second = 0"<<endl;
	      exit(-1);
	    }
	  UseExactSecond=false;
	  FirstParticles=nbrParticles;
	  cout << "Using analytic second wavefunction for product"<<endl;
	  AnalyticSecondWaveFunction = waveFunctionManager->GetWaveFunction();
	  if (AnalyticSecondWaveFunction==0)
	    {
	      cout << "Failed to option valid analytic wavefunction. Please check your wavefunction options!"<<endl;
	      exit(-1);
	    }
	  
	  ((Abstract1DComplexFunction*)AnalyticSecondWaveFunction)->AdaptAverageMCNorm(500,1000);
	}
      else
	{
	  cout << "An exact state vector is required for the second state"<<endl;
	  exit(-1);
	}
    }

  Statistics = true;
  this->SecondState = new RealVector();
      
  if (UseExactSecond)
    {
      char *TmpState = manager.GetString("second-state");
      int SecondParticles;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(TmpState,SecondParticles, SecondLzMax, SecondLz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters for second state " <<
	    TmpState << endl;
	  exit(-1);
	}
      if (SecondParticles!=NbrParticles)
	{
	  cout << "Error: first and second state have to have the same number of particles"<<endl;
	  exit(-1);
	}
      
      if (SecondLzMax>totalLzMax)
	{
	  cout << "Error: second state has to be at a lower flux than the exact state"<<endl;
	  exit(-1);
	}
      
      if (SecondState->ReadVector (TmpState) == false)
	{
	  cout << "can't open vector file " << manager.GetString("first-state") << endl;
	  exit(-1);      
	}
      else
	{
	  cout << "Read second state "<<TmpState<<" as vector with angular momentum Lz="<<SecondLz<<endl;
	}
      if (Statistics==SecondBosonicStatistic) // inverse conventions, here!
	{
	  cout << "Statistics of second vector file do not match your settings"<<endl;
	  exit(-1);
	}

      delete [] TmpState;
    }
  else
    {
      if (SecondLzMax>totalLzMax)
	{
	  cout << "Error: second state has to be at a lower flux than the exact state"<<endl;
	  exit(-1);
	}
    }


  if (FirstBosonicStatistic)
    {
      cout << "Creating bosonic space for first state with Lz="<<FirstLz<<", LzMax="<<FirstLzMax<<": ";
      FirstSpace = new BosonOnSphereShort(NbrParticles, FirstLz, FirstLzMax);
    }
  else
    {
      cout << "Creating fermionic space for first state with Lz="<<FirstLz<<", LzMax="<<FirstLzMax<<": ";
#ifdef __64_BITS__
  if (FirstLzMax <= 63)
#else
    if (FirstLzMax <= 31)
#endif
      {	
	FirstSpace = new FermionOnSphere(NbrParticles, FirstLz, FirstLzMax);
      }
    else
#ifdef __128_BIT_LONGLONG__
      if (FirstLzMax <= 126)
#else
	if (FirstLzMax <= 62)
#endif
	  {	    
	    FirstSpace = new FermionOnSphereLong(NbrParticles, FirstLz, FirstLzMax);
	  }
	else
	  {
	    cout << "States of this polarized Hilbert space cannot be represented in a single word." << endl;
	    exit(-1);
	  }
    }

  
  
  this->SecondSpace=NULL;  
  if (UseExactSecond)
    {
      if (SecondBosonicStatistic)
	{
	  cout << "Creating bosonic space for first state with Lz="<<SecondLz<<", LzMax="<<SecondLzMax<<": ";
	  SecondSpace = new BosonOnSphereShort(NbrParticles, SecondLz, SecondLzMax);
	}
      else
	{

	  cout << "Creating fermionic space for second space with Lz="<<SecondLz<<", LzMax="<<SecondLzMax<<": ";
#ifdef __64_BITS__
	  if (SecondLzMax <= 63)
#else
	    if (SecondLzMax <= 31)
#endif
	      {	
		SecondSpace = new FermionOnSphere(NbrParticles, SecondLz, SecondLzMax);
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	      if (SecondLzMax <= 126)
#else
		if (SecondLzMax <= 62)
#endif
		  {	    
		    SecondSpace = new FermionOnSphereLong(NbrParticles, SecondLz, SecondLzMax);
		  }
		else
		  {
		    cout << "States of this polarized Hilbert space cannot be represented in a single word." << endl;
		    exit(-1);
		  }
	}
    }
  
  
  
  this->FirstOneBodyBasis = new ParticleOnSphereFunctionBasis(FirstLzMax, basisType);
  this->SecondOneBodyBasis = new ParticleOnSphereFunctionBasis(SecondLzMax, basisType);
  
  LastValueFirst=0.0;
  LastValueSecond=0.0;
  
}

// copy constructor
//
// function = reference on the wave function to copy

PolarizedProductWavefunction::PolarizedProductWavefunction(const PolarizedProductWavefunction& function)
{
  this->Architecture=function.Architecture;
  this->NbrParticles=function.NbrParticles;
  this->FirstLzMax=function.FirstLzMax;
  this->FirstLz=function.FirstLz;
  this->FirstBosonicStatistic=function.FirstBosonicStatistic;
  
  this->SecondLzMax=function.SecondLzMax;
  this->SecondLz=function.SecondLz;
  this->SecondBosonicStatistic=function.SecondBosonicStatistic;
  this->TotalLz=function.TotalLz;

  this->FirstState=function.FirstState;
  this->SecondState=function.SecondState;
  this->FirstSpace=function.FirstSpace;
  this->SecondSpace=function.SecondSpace;
  this->FirstOneBodyBasis=function.FirstOneBodyBasis;
  this->SecondOneBodyBasis=function.SecondOneBodyBasis;
  this->UseExactSecond=function.UseExactSecond;
  this->AnalyticSecondWaveFunction=function.AnalyticSecondWaveFunction;
  this->AnalyticSecondTrialWaveFunction=function.AnalyticSecondTrialWaveFunction;
  this->Flag=function.Flag;
  this->LastValueFirst = function.LastValueFirst;
  this->LastValueSecond = function.LastValueSecond;
  this->JastrowPower = function.JastrowPower;
  if (JastrowPower!=0)
    {
      this->SpinorUCoordinates = new Complex[NbrParticles];
      this->SpinorVCoordinates = new Complex[NbrParticles];
    }
}

// destructor
//
PolarizedProductWavefunction::~PolarizedProductWavefunction()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete this->FirstSpace;
      if (UseExactSecond)
	{
	  delete this->SecondSpace;
	}
      delete this->FirstState;
      delete this->SecondState;
      delete this->FirstOneBodyBasis;
      delete this->SecondOneBodyBasis;
    }
  if (this->JastrowPower!=0)
    {
      delete [] SpinorUCoordinates;
      delete [] SpinorVCoordinates;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PolarizedProductWavefunction::Clone ()
{
  return new PolarizedProductWavefunction(*this);
}


void PolarizedProductWavefunction::TestSymmetries(ParticleOnSphereCollection *particles)
{
  Complex ValueFirst, ValueSecond, ValueFirst2, ValueSecond2, ValueFirst3, ValueSecond3;
  
  
  if (UseExactSecond)
    {
      QHEParticleWaveFunctionOperation Operation(SecondSpace, &(SecondState[0]), &(particles->GetPositions()),
						 this->SecondOneBodyBasis, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      ValueSecond = Operation.GetScalar();
    }
  else
    {
      ValueSecond = (*(this->AnalyticSecondWaveFunction))(particles->GetPositions());
    }
  if (FirstSpace->GetHilbertSpaceDimension()>this->MinParallel)
    {
      QHEParticleWaveFunctionOperation Operation(FirstSpace, &(FirstState[0]), &(particles->GetPositions()),
						 this->FirstOneBodyBasis, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);
      ValueFirst = Operation.GetScalar();
    }
  else ValueFirst = FirstSpace->EvaluateWaveFunction(*FirstState, particles->GetPositions(), *(this->FirstOneBodyBasis));
  
  particles->ToggleHalfHalf();
  
  // recalculate:
  if (UseExactSecond)
    {
      QHEParticleWaveFunctionOperation Operation(SecondSpace, &(SecondState[0]), &(particles->GetPositions()),
						 this->SecondOneBodyBasis, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      ValueSecond2 = Operation.GetScalar();
    }
  else
    {
      ValueSecond2 = (*(this->AnalyticSecondWaveFunction))(particles->GetPositions());
    }
      
  if (FirstSpace->GetHilbertSpaceDimension()>this->MinParallel)
    {
      QHEParticleWaveFunctionOperation Operation2(FirstSpace, &(FirstState[0]), &(particles->GetPositions()),
						  this->FirstOneBodyBasis, /* TimeCoherence */ -1);
      Operation2.ApplyOperation(this->Architecture);
      ValueFirst2 = Operation2.GetScalar();
    }
  else ValueFirst2 = FirstSpace->EvaluateWaveFunction(*FirstState, particles->GetPositions(), (*this->FirstOneBodyBasis));
  
  // rotate all particles      
  particles->RotateAll(0.781723465, 2.13428571);
  
  // recalculate:
  if (UseExactSecond)
    {
      QHEParticleWaveFunctionOperation Operation(SecondSpace, &(SecondState[0]), &(particles->GetPositions()),
						 this->SecondOneBodyBasis, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      ValueSecond3 = Operation.GetScalar();
    }
  else
    {
      ValueSecond3 = (*(this->AnalyticSecondWaveFunction))(particles->GetPositions());
    }
  if (FirstSpace->GetHilbertSpaceDimension()>this->MinParallel)
    {
      QHEParticleWaveFunctionOperation Operation3(FirstSpace, &(FirstState[0]), &(particles->GetPositions()),
						  this->FirstOneBodyBasis, /* TimeCoherence */ -1);
      Operation3.ApplyOperation(this->Architecture);
      ValueFirst3 = Operation3.GetScalar();
    }
  else ValueFirst3 = FirstSpace->EvaluateWaveFunction(*FirstState, particles->GetPositions(), *(this->FirstOneBodyBasis));
  
  cout << "First Before exchange: "<< ValueSecond << endl << "After exchange:  " << ValueSecond2 << endl;
  cout << "First Parity: " << ValueSecond/ValueSecond2 << endl;
  cout << "First After rotation: " << ValueSecond3 << " ratio: "<< Norm(ValueSecond/ValueSecond3) << endl;
  cout << "Second Before exchange: "<< ValueFirst  << endl << "After exchange:  " << ValueFirst2 << endl;
  cout << "Second Parity: " << ValueFirst / ValueFirst2 << endl;
  cout << "Second After rotation: " << ValueFirst3 << " ratio: "<< Norm(ValueFirst/ValueFirst3) << endl;
}


// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PolarizedProductWavefunction::operator ()(RealVector& x)
{
  this->Jastrow=1.0;
  if (this->JastrowPower!=0)
    {
      double s,c;
      // Convert particle positions into spinor form:
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
	}
      
      double Factor = M_PI * 0.5;
      for (int i=1;i<this->NbrParticles;++i)
	for(int j=0;j<i;++j)
	  Jastrow *= Factor*((this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]));

      if (this->JastrowPower<0) Jastrow=1.0/Jastrow;
      if (abs(this->JastrowPower)>1)
	{
	  Complex Tmp=Jastrow;
	  for (int j=abs(this->JastrowPower); j>1; --j) Jastrow*=Tmp;
	}
    }

  Complex ValueFirst, ValueSecond;
  if (FirstSpace->GetHilbertSpaceDimension()>this->MinParallel)
    {
      QHEParticleWaveFunctionOperation Operation(FirstSpace, &(FirstState[0]), &x,
						 this->FirstOneBodyBasis, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);
      ValueFirst = Operation.GetScalar();
    }
  else ValueFirst = FirstSpace->EvaluateWaveFunction(*FirstState, x, *(this->FirstOneBodyBasis));

  if (UseExactSecond)
    {
      if (SecondSpace->GetHilbertSpaceDimension()>this->MinParallel)
	{
	  QHEParticleWaveFunctionOperation Operation(SecondSpace, &(SecondState[0]), &x,
						     this->SecondOneBodyBasis, /* TimeCoherence */ -1);
	  Operation.ApplyOperation(this->Architecture);      
	  ValueSecond = Operation.GetScalar();
	}
      else
	ValueSecond = SecondSpace->EvaluateWaveFunction(*SecondState, x, *(this->SecondOneBodyBasis));
    }
  else
    {
      ValueSecond = (*(this->AnalyticSecondWaveFunction))(x);
    }
  
  return Jastrow*ValueFirst*ValueSecond;
}

// get a value of the wavefunction for the last set of coordinates, but with different variational parameters
// parameters =  alternative set of parameters
Complex PolarizedProductWavefunction::GetForOtherParameters(double *parameters)
{
  if ((UseExactSecond)|(AnalyticSecondTrialWaveFunction==NULL))
    {
      cout << "Cannot change any parameters in PolarizedProductWavefunction when calculating from exact state"
	   << ", or using non variational wavefunction."<<endl;
      exit(1);
    }
  Complex Result=this->AnalyticSecondTrialWaveFunction->GetForOtherParameters(parameters);
  
  return Jastrow*Result*this->LastValueFirst;
}

// do many evaluations of the function, storing the result in the vector results given in the call
// result = vector of leading dimension of the array coefficients for returning values
// x = positions to evaluate the wavefuntion in
// format for passing parameters as [nbrSet][nbrParameter],
void PolarizedProductWavefunction::GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients)
{
  if ((UseExactSecond)|(AnalyticSecondTrialWaveFunction==NULL))
    {
      cout << "Cannot change any parameters in PolarizedProductWavefunction when calculating from exact state"
	   << ", or using non variational wavefunction."<<endl;
      exit(1);
    }
  this->AnalyticSecondTrialWaveFunction->GetForManyParameters(results,x,coefficients);

  Complex FirstPart;
  if (FirstSpace->GetHilbertSpaceDimension()>this->MinParallel)
    {  
      QHEParticleWaveFunctionOperation Operation(FirstSpace, &(FirstState[0]), &x,
						 this->FirstOneBodyBasis, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);
      FirstPart = Operation.GetScalar();
    }
  else FirstPart = FirstSpace->EvaluateWaveFunction(*FirstState, x, *(this->FirstOneBodyBasis));
  for (int i=0; i<results.GetVectorDimension(); ++i)
    results[i]*=Jastrow*FirstPart;
}

// access internal values of parameters
double *PolarizedProductWavefunction::GetTrialParameters()
{
  if ((UseExactSecond)|(AnalyticSecondTrialWaveFunction==NULL))
    {
      cout << "Cannot get any parameters in PolarizedProductWavefunction when calculating from exact state"
	   << ", or using non variational wavefunction."<<endl;
      exit(1);
    }
  return this->AnalyticSecondTrialWaveFunction->GetTrialParameters();
}

// get number of parameters
int PolarizedProductWavefunction::GetNbrParameters() 
{
  if ((UseExactSecond)|(AnalyticSecondTrialWaveFunction==NULL))
    {
      cout << "No parameters available in PolarizedProductWavefunction when calculating from exact state"
	   << ", or using non variational wavefunction."<<endl;
      exit(1);
    }
  return this->AnalyticSecondTrialWaveFunction->GetNbrParameters();
}

  
// set new values of the trial coefficients (keeping the initial number of parameters)
void PolarizedProductWavefunction::SetTrialParameters(double * coefficients)
{
  if ((UseExactSecond)|(AnalyticSecondTrialWaveFunction==NULL))
    {
      cout << "Cannot change any parameters in PolarizedProductWavefunction when calculating from exact state"
	   << ", or using non variational wavefunction."<<endl;
      exit(1);
    }
  this->AnalyticSecondTrialWaveFunction->SetTrialParameters(coefficients);
}


// get function properties, and possible extensions of interface 
// 
unsigned PolarizedProductWavefunction::GetProperties()
{
  if (UseExactSecond)
    return Abstract1DComplexFunction::Basic;
  else
    return this->AnalyticSecondWaveFunction->GetProperties();
}




// add an option group containing all options related to the skyrmion wave functions
//
// manager = pointer to the option manager
void PolarizedProductWavefunction::AddPolarizedProductStateOptionGroup(OptionManager &manager)
{  
  OptionGroup* PolarizedGroup = new OptionGroup ("product state options");
  manager+=PolarizedGroup;
  (*PolarizedGroup) += new BooleanOption ('\n', "product-state", "Use the product state interface");
  (*PolarizedGroup) += new SingleStringOption  ('\n', "first-state", "file name of the first wave function");
  (*PolarizedGroup) += new BooleanOption ('\n', "first-boson", "Use a polarized bosonic basis for first reference state");
  
  (*PolarizedGroup) += new SingleIntegerOption  ('\n', "first-lzmax", "total number of flux quanta (0 if it has to be guessed from input file name)", 0);  
  (*PolarizedGroup) += new SingleIntegerOption  ('\n', "first-lz", "twice the total lz value of the polarized vector(s) (0 if it has to be guessed from input file name)",0);

  (*PolarizedGroup) += new SingleStringOption  ('\n', "second-state", "file name of the second wave function (optional)");
  (*PolarizedGroup) += new BooleanOption ('\n', "second-boson", "Use a polarized bosonic basis for first reference state");
  
  (*PolarizedGroup) += new SingleIntegerOption  ('\n', "second-lzmax", "total number of flux quanta (0 if it has to be guessed from input file name)", 0);  
  (*PolarizedGroup) += new SingleIntegerOption  ('\n', "second-lz", "twice the total lz value of the polarized vector(s) (0 if it has to be guessed from input file name)", 0);
  (*PolarizedGroup) += new SingleIntegerOption  ('\n', "min-parallel", "minimum dimension of polarized space before parallel evaluation is applied", 5000);

  (*PolarizedGroup) += new SingleIntegerOption  ('\n', "jastrow-power", "powers of Jastrow factors to be added (can be negative)",0);

}
