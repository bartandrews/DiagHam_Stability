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
#include "MathTools/ClebschGordanCoefficients.h"

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
  this->PolarizedLz = manager.GetIntegers("polarized-lz",this->NbrMultipletPolarized);
  this->PolarizedL = manager.GetInteger("polarizedL");
  this->MinParallel = manager.GetInteger("min-parallel");  
  this->AnalyticPolarizedWaveFunction=NULL;
  this->AnalyticPolarizedTrialWaveFunction=NULL;

  this->TotalL=manager.GetInteger("skyrmionL");
  this->TotalLz=manager.GetInteger("skyrmion-lz");

  if (abs(this->TotalLz)>this->TotalL)
    {
      cout << "Invalid TotalLz: require |Lz| <= L"<<endl;
      exit(-1);
    }  
  
  this->Flag.Initialize();  

  if (manager.GetStrings("polarized-state")==NULL)
    {
      if (waveFunctionManager!=NULL)
	{
	  if ((this->NbrMultipletPolarized!=1)||(this->PolarizedLz[0]!=0))
	    {
	      cout << "Analytic polarized wavefunction can only be generated for cases with L_pol = 0"<<endl;
	      exit(-1);
	    }
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
      int tmpNbrStates;
      char **TmpStates = manager.GetStrings("polarized-state", tmpNbrStates);
      this->PolarizedState = new RealVector[this->NbrMultipletPolarized];
      if (tmpNbrStates != this->NbrMultipletPolarized)
	{
	  cout << "please indicate angular momentum component Lz for all trial states in multiplet"<<endl;
	  exit(-1);
	}
      for (int i=0; i<tmpNbrStates; ++i)
	{
	  if (FQHEOnSphereFindSystemInfoFromVectorFileName(TmpStates[i],
							   PolarizedParticles, PolarizedLzMax, PolarizedLz[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters for polarized state " <<
		TmpStates[i] << endl;
	      exit(-1);
	    }
	  if (PolarizedParticles!=nbrParticles)
	    {
	      cout << "Error: polarized and exact state have to have the same number of particles"<<endl;
	      exit(-1);
	    }
      
	  if (PolarizedLzMax>totalLzMax)
	    {
	      cout << "Error: polarized state has to be at a lower flux than the exact state"<<endl;
	      exit(-1);
	    }
      
	  if (PolarizedState[i].ReadVector (TmpStates[i]) == false)
	    {
	      cout << "can't open vector file " << manager.GetString("polarized-state") << endl;
	      exit(-1);      
	    }
	  else
	    {
	      cout << "Read polarized state "<<TmpStates[i]<<" as vector "<<i<<" with angular momentum Lz="<<PolarizedLz[i]<<endl;
	    }

	}
    }
  else
    {
      if (PolarizedLzMax>totalLzMax)
	{
	  cout << "Error: polarized state has to be at a lower flux than the exact state"<<endl;
	  exit(-1);
	}
    }


  int NbrBosons = manager.GetInteger("nbr-particles"); 
  this->BosonLzMax = 0;
  this->BosonLz = manager.GetIntegers("bosonic-lz",NbrMultipletBosons);
  this->BosonL = manager.GetInteger("bosonicL");
  this->BosonSz = manager.GetInteger("bosonic-sz");
  bool SzSymmetrizedBasis = false;
  bool SzMinusParity = false;
  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  Statistics = false;  

  int tmpNbrStates;
  char **TmpStates = manager.GetStrings("bosonic-state", tmpNbrStates);
  if (tmpNbrStates != this->NbrMultipletBosons)
    {
      cout << "please indicate angular momentum component Lz for all trial states in multiplet"<<endl;
      exit(-1);
    }
  this->BosonicState = new RealVector[this->NbrMultipletBosons];
  for (int i=0; i<tmpNbrStates; ++i)
    {
      if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(TmpStates[i], NbrBosons,
							       BosonLzMax, BosonLz[i], BosonSz,
							       SzSymmetrizedBasis, SzMinusParity, 
							       LzSymmetrizedBasis, LzMinusParity, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << manager.GetString("bosonic-state") << endl;
	  exit(-1);
	}

      if (NbrBosons!=nbrParticles)
	{
	  cout << "Error: bosonic and exact state have to have the same number of particles"<<endl;
	  exit(-1);
	}
  
      if (BosonLzMax==0)
	BosonLzMax = totalLzMax-PolarizedLzMax;
      else if (BosonLzMax!=totalLzMax-PolarizedLzMax)
	{
	  cout << "Error: total flux has to match: BosonLzMax == TotalLzMax-PolarizedLzMax"<<endl;
	  exit(-1);
	}
      
      if (tmpNbrStates==1)
	{
	  if (BosonLz[0]==0)
	    BosonLz[0] = totalLz;
	  else if (BosonLz[0]!=totalLz)
	    {
	      cout << "Error: total angular momentum has to match: BosonLz == TotalLz"<<endl;
	      exit(-1);
	    }
	}
      
      if (BosonSz==0)
	BosonSz = totalSz;
      else if (BosonSz!=totalSz)
	{
	  cout << "Error: total spin has to match: BosonSz == TotalSz ("<<BosonSz<<" vs "<<totalSz<<")"<<endl;
	  exit(-1);
	}
      
      if (BosonicState[i].ReadVector (TmpStates[i]) == false)
	{
	  cout << "can't open vector file " << TmpStates[i] << endl;
	  exit(-1);      
	}
      else
	{
	  cout << "Read bosonic state "<<TmpStates[i]<<" as vector "<<i<<" with angular momentum Lz="<<BosonLz[i]<<endl;
	}
    }
  
  
  if ((this->BosonL==0)&&(this->PolarizedL==0))
    {
      this->NbrCoupling=1;
      this->BosonIndex = new int[1];
      this->BosonIndex[0]=0;
      this->PolarizedIndex = new int[1];
      this->PolarizedIndex[0]=0;
      this->CouplingForIndex = new double[1];
      this->CouplingForIndex[0]=0.0;
    }
  else
    {
      // set up coupling of angular momenta
      ClebschGordanCoefficients Clebsch(this->PolarizedL, this->BosonL);
      // count number of relevant couplings
      this->NbrCoupling=0;
      cout << "Coupling "<<this->PolarizedL << "_P + "<<this->BosonL<<"_B"<<" = "<<this->TotalL<<"_tot at Lz_tot = "
	   << this->TotalLz << " -- relevant couplings: ";
      for (int mP=-this->PolarizedL; mP<=this->PolarizedL; mP+=2)
	for (int mB=-this->BosonL; mB<=this->BosonL; mB+=2)
	  if ((mP+mB==this->TotalLz)&&(Clebsch.CarefulGetCoefficient (mP, mB, this->TotalL)!=0.0))
	    ++NbrCoupling;
      cout << NbrCoupling<<endl;
      this->BosonIndex = new int[NbrCoupling];
      this->PolarizedIndex = new int[NbrCoupling];
      this->CouplingForIndex = new double[NbrCoupling];
      // assign terms in expansion
      this->NbrCoupling=0;
      for (int mP=-this->PolarizedL; mP<=this->PolarizedL; mP+=2)
	for (int mB=-this->BosonL; mB<=this->BosonL; mB+=2)
	  if ((mP+mB==this->TotalLz)&&(Clebsch.CarefulGetCoefficient (mP, mB, this->TotalL)!=0.0))
	    {
	      this->CouplingForIndex[NbrCoupling] = Clebsch.CarefulGetCoefficient (mP, mB, this->TotalL);
	      this->BosonIndex[NbrCoupling] = -1;
	      for (int i=0; i<NbrMultipletBosons; ++i)
		if (BosonLz[i]==mB)
		  this->BosonIndex[NbrCoupling] = i;
	      if (this->BosonIndex[NbrCoupling] == -1)
		{
		  cout << "Require bosonic state with angular momentum lz="<<mB<<endl;
		  exit(-1);
		}
	      this->PolarizedIndex[NbrCoupling] = -1;
	      for (int i=0; i<NbrMultipletPolarized; ++i)
		if (PolarizedLz[i]==mP)
		  this->PolarizedIndex[NbrCoupling] = i;
	      if (this->PolarizedIndex[NbrCoupling] == -1)
		{
		  cout << "Require polarized state with angular momentum lz="<<mB<<endl;
		  exit(-1);
		}
	      cout << "Coupling "<<NbrCoupling<<": "<<CouplingForIndex[NbrCoupling] << " * | ("<<this->BosonL<<", "
		   << BosonLz[BosonIndex[NbrCoupling]]<<")_B, ("<<this->PolarizedL<<", "
		   << PolarizedLz[PolarizedIndex[NbrCoupling]]<<")_P >"<<"  (vectors "<<BosonIndex[NbrCoupling]<< "_B"
		   << " and "<<PolarizedIndex[NbrCoupling] <<"_P)"<< endl;
	      ++NbrCoupling;
	    }
    }

  this->PolarizedSpace=new ParticleOnSphere*[NbrMultipletPolarized];
  for (int i=0; i<NbrMultipletPolarized; ++i)
    this->PolarizedSpace[i]=NULL;  
  if (UseExact)
    {
      for (int i=0; i<NbrCoupling; ++i)
	{
	  cout << "Creating polarized space for Lz="<<PolarizedLz[PolarizedIndex[i]]
	       <<" at position "<<PolarizedIndex[i]<<": ";
#ifdef __64_BITS__
	  if (PolarizedLzMax <= 63)
#else
	    if (PolarizedLzMax <= 31)
#endif
	      {	
		PolarizedSpace[PolarizedIndex[i]] =
		  new FermionOnSphere(NbrParticles, PolarizedLz[PolarizedIndex[i]], PolarizedLzMax);
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	      if (PolarizedLzMax <= 126)
#else
		if (PolarizedLzMax <= 62)
#endif
		  {	    
		    PolarizedSpace[PolarizedIndex[i]] =
		      new FermionOnSphereLong(NbrParticles, PolarizedLz[PolarizedIndex[i]], PolarizedLzMax);
		  }
		else
		  {
		    cout << "States of this polarized Hilbert space cannot be represented in a single word." << endl;
		    exit(-1);
		  }
	}
    }
  // clear up unnecessary vectors
  if (UseExact)
    for (int i=0; i<NbrMultipletPolarized; ++i)
      if (PolarizedSpace[i]==NULL)
	this->PolarizedState[i]=RealVector();  
  
  this->BosonicSpace = new ParticleOnSphereWithSpin*[NbrMultipletBosons];
  for (int i=0; i<NbrMultipletBosons; ++i)
    this->BosonicSpace[i] = NULL;
  for (int i=0; i<NbrCoupling; ++i)
    {
      cout << "Creating bosonic space for Lz="<<BosonLz[i]<<" at position "<<BosonIndex[i]<<": ";
      this->BosonicSpace[BosonIndex[i]]= new BosonOnSphereWithSpin(NbrParticles, BosonLz[BosonIndex[i]], BosonLzMax, BosonSz);
    }
  for (int i=0; i<NbrMultipletBosons; ++i)
    if (BosonicSpace[i]==NULL)
      this->BosonicState[i]=RealVector();  
  
  
  this->OneBodyBasisPol = new ParticleOnSphereFunctionBasis(PolarizedLzMax, basisType);
  this->OneBodyBasisBos = new ParticleOnSphereFunctionBasis(BosonLzMax, basisType);

  this->LastValueBosonic = new Complex[NbrCoupling];
  this->TmpValuePolarized = new Complex[NbrCoupling];
  for (int i=0; i<NbrCoupling; ++i)
    {
      this->LastValueBosonic[i]=0.0;
      this->TmpValuePolarized[i]=0.0;
    }
}

// copy constructor
//
// function = reference on the wave function to copy

SkyrmionOnSphereWaveFunction::SkyrmionOnSphereWaveFunction(const SkyrmionOnSphereWaveFunction& function)
{
  this->NbrParticles=function.NbrParticles;
  this->PolarizedLzMax=function.PolarizedLzMax;
  this->PolarizedLz=function.PolarizedLz;
  this->PolarizedL=function.PolarizedL;
  this->BosonLzMax=function.BosonLzMax;
  this->BosonLz=function.BosonLz;
  this->BosonSz=function.BosonSz;
  this->BosonL=function.BosonL;
  this->TotalL=function.TotalL;
  this->TotalLz=function.TotalLz;
  this->NbrCoupling=function.NbrCoupling;
  this->BosonIndex=function.BosonIndex;
  this->PolarizedIndex=function.PolarizedIndex;
  this->CouplingForIndex=function.CouplingForIndex;
  this->NbrMultipletPolarized=function.NbrMultipletPolarized;
  this->NbrMultipletBosons=function.NbrMultipletBosons;
  this->PolarizedState=function.PolarizedState;
  this->BosonicState=function.BosonicState;
  this->PolarizedSpace=function.PolarizedSpace;
  this->BosonicSpace=function.BosonicSpace;
  this->OneBodyBasisBos=function.OneBodyBasisBos;
  this->OneBodyBasisPol=function.OneBodyBasisPol;
  this->UseExact=function.UseExact;
  this->AnalyticPolarizedWaveFunction=function.AnalyticPolarizedWaveFunction;
  this->Flag=function.Flag;
  this->LastValueBosonic = function.LastValueBosonic;
  this->TmpValuePolarized = function.TmpValuePolarized;
}

// destructor
//
SkyrmionOnSphereWaveFunction::~SkyrmionOnSphereWaveFunction()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (UseExact)
	{
	  for (int i=0; i<NbrMultipletPolarized; ++i)
	    if (this->PolarizedSpace[i]!=NULL)
	      delete this->PolarizedSpace[i];
	}
      delete[] this->PolarizedSpace;
      delete[] this->BosonLz;
      delete[] this->PolarizedLz;
      delete[] this->BosonIndex;
      delete[] this->PolarizedIndex;
      delete[] this->CouplingForIndex;
      for (int i=0; i<NbrMultipletBosons; ++i)
	if (this->BosonicSpace[i]!=NULL)
	  delete this->BosonicSpace[i];
      delete[] this->BosonicSpace;
      delete[] this->BosonicState;
      delete[] this->PolarizedState;
      delete this->OneBodyBasisPol;
      delete this->OneBodyBasisBos;
    }
  delete [] this->LastValueBosonic;
  delete [] this->TmpValuePolarized;
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
  if ((this->NbrMultipletPolarized==1)&&(this->NbrMultipletBosons==1))
    {
      Complex ValueSpin, ValuePolarized, ValueSpin2, ValuePolarized2, ValueSpin3, ValuePolarized3;
  
      if (UseExact)
	{
	  QHEParticleWaveFunctionOperation Operation(*PolarizedSpace, &(PolarizedState[0]), &(particles->GetPositions()),
						     this->OneBodyBasisPol, /* TimeCoherence */ -1);
	  Operation.ApplyOperation(this->Architecture);      
	  ValuePolarized = Operation.GetScalar();
	}
      else
	{
	  ValuePolarized = (*(this->AnalyticPolarizedWaveFunction))(particles->GetPositions());
	}
      if (BosonicSpace[0]->GetHilbertSpaceDimension()>this->MinParallel)
	{
	  QHEParticleWaveFunctionOperation Operation(*BosonicSpace, &(BosonicState[0]), &(particles->GetPositions()),
						     this->OneBodyBasisBos, /* TimeCoherence */ -1);
	  Operation.ApplyOperation(this->Architecture);      
	  ValueSpin = Operation.GetScalar();
	}
      else ValueSpin = BosonicSpace[0]->EvaluateWaveFunction(BosonicState[0], particles->GetPositions(), *(this->OneBodyBasisBos));
  
      particles->ToggleHalfHalf();

      // recalculate:
      if (UseExact)
	{
	  QHEParticleWaveFunctionOperation Operation(*PolarizedSpace, &(PolarizedState[0]), &(particles->GetPositions()),
						     this->OneBodyBasisPol, /* TimeCoherence */ -1);
	  Operation.ApplyOperation(this->Architecture);      
	  ValuePolarized2 = Operation.GetScalar();
	}
      else
	{
	  ValuePolarized2 = (*(this->AnalyticPolarizedWaveFunction))(particles->GetPositions());
	}
      
      if (BosonicSpace[0]->GetHilbertSpaceDimension()>this->MinParallel)
	{
	  QHEParticleWaveFunctionOperation Operation2(*BosonicSpace, &(BosonicState[0]), &(particles->GetPositions()),
						      this->OneBodyBasisBos, /* TimeCoherence */ -1);
	  Operation2.ApplyOperation(this->Architecture);      
	  ValueSpin2 = Operation2.GetScalar();
	}
      else ValueSpin2 = BosonicSpace[0]->EvaluateWaveFunction(BosonicState[0], particles->GetPositions(), (*this->OneBodyBasisBos));
  
      // rotate all particles      
      particles->RotateAll(0.781723465, 2.13428571);

      // recalculate:
      if (UseExact)
	{
	  QHEParticleWaveFunctionOperation Operation(*PolarizedSpace, &(PolarizedState[0]), &(particles->GetPositions()),
						     this->OneBodyBasisPol, /* TimeCoherence */ -1);
	  Operation.ApplyOperation(this->Architecture);      
	  ValuePolarized3 = Operation.GetScalar();
	}
      else
	{
	  ValuePolarized3 = (*(this->AnalyticPolarizedWaveFunction))(particles->GetPositions());
	}
      if (BosonicSpace[0]->GetHilbertSpaceDimension()>this->MinParallel)
	{
	  QHEParticleWaveFunctionOperation Operation3(*BosonicSpace, &(BosonicState[0]), &(particles->GetPositions()),
						      this->OneBodyBasisBos, /* TimeCoherence */ -1);
	  Operation3.ApplyOperation(this->Architecture);
	  ValueSpin3 = Operation3.GetScalar();
	}
      else ValueSpin3 = BosonicSpace[0]->EvaluateWaveFunction(BosonicState[0], particles->GetPositions(), *(this->OneBodyBasisBos));
  
      cout << "Pol Before exchange: "<< ValuePolarized << endl << "After exchange:  " << ValuePolarized2 << endl;
      cout << "Pol Parity: " << ValuePolarized/ValuePolarized2 << endl;
      cout << "Pol After rotation: " << ValuePolarized3 << " ratio: "<< Norm(ValuePolarized/ValuePolarized3) << endl;
      cout << "Bos Before exchange: "<< ValueSpin  << endl << "After exchange:  " << ValueSpin2 << endl;
      cout << "Bos Parity: " << ValueSpin / ValueSpin2 << endl;
      cout << "Bos After rotation: " << ValueSpin3 << " ratio: "<< Norm(ValueSpin/ValueSpin3) << endl;
    }

}


// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex SkyrmionOnSphereWaveFunction::operator ()(RealVector& x)
{
  Complex sum=0.0;
  for (int c=0; c<NbrCoupling; ++c)
    {
      int IndexP=PolarizedIndex[c];
      int IndexB=BosonIndex[c];
      if (UseExact) 
	{	  
	  QHEParticleWaveFunctionOperation Operation(PolarizedSpace[IndexP], &(PolarizedState[IndexP]), &x,
						     this->OneBodyBasisPol, /* TimeCoherence */ -1);
	  Operation.ApplyOperation(this->Architecture);
	  TmpValuePolarized[IndexP] = Operation.GetScalar();
	}
      else // only occurs if we have a single component
	{
	  TmpValuePolarized[IndexP] = (*(this->AnalyticPolarizedWaveFunction))(x);
	}
      if (BosonicSpace[IndexB]->GetHilbertSpaceDimension()>this->MinParallel)
	{  
	  QHEParticleWaveFunctionOperation Operation(BosonicSpace[IndexB], &(BosonicState[IndexB]), &x,
						     this->OneBodyBasisBos, /* TimeCoherence */ -1);
	  Operation.ApplyOperation(this->Architecture);      
	  this->LastValueBosonic[IndexB] = Operation.GetScalar();
	}
      else this->LastValueBosonic[IndexB] = BosonicSpace[IndexB]->EvaluateWaveFunction(BosonicState[IndexB],
									       x, *(this->OneBodyBasisBos));
      sum+=(CouplingForIndex[c]*(TmpValuePolarized[IndexP]*this->LastValueBosonic[IndexB]));
    }
  
  return sum;
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
  
  return Result*this->LastValueBosonic[0];
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
  if (BosonicSpace[0]->GetHilbertSpaceDimension()>this->MinParallel)
    {  
      QHEParticleWaveFunctionOperation Operation(BosonicSpace[0], &(BosonicState[0]), &x,
						 this->OneBodyBasisBos, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      BosonicPart = Operation.GetScalar();
    }
  else BosonicPart = BosonicSpace[0]->EvaluateWaveFunction(BosonicState[0], x, *(this->OneBodyBasisBos));
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
  
  (*SkyrmionGroup) += new MultipleStringOption  ('\n', "polarized-state", "file name of polarized fermionic reference wave function (if omitted using analytic function)");
  (*SkyrmionGroup) += new SingleIntegerOption  ('L', "polarizedL", "twice the angular momentum of the polarized state", 0, true, 0);
  (*SkyrmionGroup) += new SingleIntegerOption  ('l', "polarized-lzmax", "total number of flux quanta (0 if it has to be guessed from input file name)", 0);  
  (*SkyrmionGroup) += new MultipleIntegerOption  ('z', "polarized-lz", "twice the total lz value of the polarized vector(s) (0 if it has to be guessed from input file name)", ',',',',"0");
  (*SkyrmionGroup) += new MultipleStringOption  ('\n', "bosonic-state", "file name of spinful bosonic part of wave function");
  (*SkyrmionGroup) += new MultipleIntegerOption  ('\n', "bosonic-lz", "twice the total lz value of the bosonic vector(s) (0 if it has to be guessed from input file name)", ',',',',"0");
  (*SkyrmionGroup) += new SingleIntegerOption  ('\n', "bosonic-sz", "twice the total sz value of the bosonic vector(s) (0 if it has to be guessed from input file name)", 0);
  (*SkyrmionGroup) += new SingleIntegerOption  ('\n', "bosonicL", "twice the angular momentum of the polarized state", 0);
  (*SkyrmionGroup) += new SingleIntegerOption  ('\n', "skyrmionL", "twice the angular momentum of the skyrmion trial state", 0);
  (*SkyrmionGroup) += new SingleIntegerOption  ('\n', "skyrmion-lz", "twice the angular momentum projection Lz of the skyrmion trial state", 0);
  (*SkyrmionGroup) += new SingleIntegerOption  ('\n', "min-parallel", "minimum dimension of bosonic space before parallel evaluation is applied", 5000);

  if (wfManager!=NULL)
    wfManager->AddOptionGroup(&manager);
}
