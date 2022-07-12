////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//      class for a basic Monte Carlo algorith for particles on a sphere      //
//                                                                            //
//                        last modification : 23/01/2008                      //
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


#include "MRBlockSamplingFunction.h"
#include "ParticleOnSphereCollection.h"
#include "MathTools/FactorialCoefficient.h"
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;

#ifdef __64_BITS__
#define ULENGTH 64
#define USHIFT 6
#else
#define ULENGTH 32
#define USHIFT 5
#endif


// constructor
// nbrParticlesPerBlock = number of particles in each of the two blocks
// jastrowExponent = exponent of jastrow factor multiplying state
// sqrCriticalDistance = tuning parameter for recalculation of blocks
MRBlockSamplingFunction::MRBlockSamplingFunction(int nbrParticlesPerBlock, int jastrowExponent, double sqrCriticalDistance)
{
  this->NbrParticlesPerBlock=nbrParticlesPerBlock;
  this->NbrParticles=2*NbrParticlesPerBlock;
  this->JastrowExponent=jastrowExponent;
  this->System=NULL;
  this->ElementNorm=1.0;

  this->NbrBlocks = (nbrParticlesPerBlock+2)/2;
  this->BlockPermutations = new int*[NbrBlocks];
  for (int i=0; i<NbrBlocks; ++i)
    {
      this->BlockPermutations[i] = new int[NbrParticles];
      for (int n=0; n<i; ++n)
	{
	  this->BlockPermutations[i][n]=nbrParticlesPerBlock+n;
	  this->BlockPermutations[i][nbrParticlesPerBlock+n]=n;
	}
      for (int n=i; n<nbrParticlesPerBlock; ++n)
	{
	  this->BlockPermutations[i][n]=n;
	  this->BlockPermutations[i][nbrParticlesPerBlock+n]=nbrParticlesPerBlock+n;
	}
    }
  this->BlockWeights = new double[NbrBlocks];
  FactorialCoefficient Prefactor;
  Prefactor.SetToOne();
  Prefactor.FactorialMultiply(NbrParticlesPerBlock);
  Prefactor.FactorialMultiply(NbrParticlesPerBlock);
  Prefactor.FactorialDivide(NbrParticles);
  Prefactor.FactorialMultiply(NbrParticlesPerBlock);
  Prefactor.FactorialMultiply(NbrParticlesPerBlock);
  for (int k=0; k<NbrBlocks; ++k)
    {
      FactorialCoefficient Factor;
      Factor = Prefactor;
      Factor.FactorialDivide(k);
      Factor.FactorialDivide(k);
      Factor.FactorialDivide(NbrParticlesPerBlock-k);
      Factor.FactorialDivide(NbrParticlesPerBlock-k);
      this->BlockWeights[k] = Factor.GetNumericalValue();
    }
  this->SqrCriticalDistance = sqrCriticalDistance;
  this->JastrowElements = new Complex[NbrParticles*(NbrParticles+1)/2];
  this->CriticalParticles = new unsigned long[NbrParticles>>USHIFT];
}

// virtual destructor
MRBlockSamplingFunction::~MRBlockSamplingFunction()
{
  if (NbrParticles!=0)
    {
      for (int i=0; i<NbrBlocks; ++i)
	delete [] this->BlockPermutations[i];
      delete [] this->BlockPermutations;
      delete [] this->BlockWeights;
      delete [] this->CriticalParticles;
    }
}


// register basic system of particles
// this function needs to be called before any of the other routines are functional
void MRBlockSamplingFunction::RegisterSystem(AbstractParticleCollectionOnSphere *system)
{
  this->System=system;
  if (System->GetNbrParticles() != this->NbrParticles)
    {
      cout << "Number of particles in system not compatible in sampling function";
      exit(1);
    }
  // pointers to spinor coordinates (external)
  System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
}



// method for ratio of probabilities with respect to the last configuration
// allows for more rapid calculation due to cancellation of factors
double MRBlockSamplingFunction::GetTransitionRatio()
{
  double ratio=1.0;
  int tomove = System->GetMovedNbr();
  System->GetPreviousPos(LastU,LastV);
  // block containing moved particle
  int LowerLimit=0;
  int UpperLimit=NbrParticles;
  if (tomove<NbrParticlesPerBlock)
    UpperLimit=NbrParticlesPerBlock;
  else
    LowerLimit=NbrParticlesPerBlock;
  for (int i=LowerLimit;i<tomove;i++)
    {
      ratio *= SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	SqrNorm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
    }
  for (int i=tomove+1;i<UpperLimit;i++)
    {
      ratio *= SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	SqrNorm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
    }
  ratio *= ratio;
  // global Laughlin part
  for (int i=0;i<tomove;i++)
    {
      ratio *= SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	SqrNorm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
    }
  for (int i=tomove+1;i<this->NbrParticles;i++)
    {
      ratio *= SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[tomove]-SpinorUCoordinates[tomove]*SpinorVCoordinates[i])/
	SqrNorm(SpinorUCoordinates[i]*LastV-LastU*SpinorVCoordinates[i]);
    }
  
  return sqrt(ratio);
}


  // get the estimate of the full function value calculated over all Blocks for the given system of particles
Complex MRBlockSamplingFunction::GetFunctionValue()
{
  //cout << "Attention: MRBlockSamplingFunction::GetFunctionValue() yields the sampling block, only!"<<endl;
  Complex Result = 1.0;
  for (int i=1;i<NbrParticlesPerBlock;++i)
    for (int j=0;j<i;++j)
      Result *= this->ElementNorm*(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
  for (int i=NbrParticlesPerBlock;i<NbrParticles;++i)
    for (int j=NbrParticlesPerBlock;j<i;++j)
      Result *= this->ElementNorm*(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);

  Result *= Result;

  for (int j=1; j<this->NbrParticles; ++j)
    for (int i=0;i<j;i++)
      Result *= this->ElementNorm*(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
  
  return Result;

}

// get number of flux quanta, or degree of underlying polynomial for simulation wavefunction
int MRBlockSamplingFunction::GetNbrFluxQuanta()
{
  return (2*(NbrParticlesPerBlock -1) + JastrowExponent*(NbrParticles-1));
}


// get the estimate of the sampling function value calculated over the sampling block only, for the given system of particles
Complex MRBlockSamplingFunction::GetSamplingBlockValue()
{
  Complex Result = 1.0;
  for (int i=1;i<NbrParticlesPerBlock;++i)
    for (int j=0;j<i;++j)
      Result *= this->ElementNorm*(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
  for (int i=NbrParticlesPerBlock;i<NbrParticles;++i)
    for (int j=NbrParticlesPerBlock;j<i;++j)
      Result *= this->ElementNorm*(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);

  Result *= Result;

  for (int j=1; j<this->NbrParticles; ++j)
    for (int i=0;i<j;i++)
      Result *= this->ElementNorm*(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
  
  return Result;
}


// call this method to scale the sampling function (needed to normalize the function)
// scale = total scaling factor
void MRBlockSamplingFunction::ScaleByFactor(double scale)
{
  double factors = (double)NbrParticles*(NbrParticles-1)/2.0+2.0*(double)NbrParticlesPerBlock*(NbrParticlesPerBlock-1);
  this->ElementNorm *= pow(scale,1.0/factors);
}


// get the current phase of the sampling function
Complex MRBlockSamplingFunction::GetSamplingBlockPhase()
{
  Complex Tmp, Result = 1.0;
  for (int i=1;i<NbrParticlesPerBlock;++i)
    for (int j=0;j<i;++j)
      {
	Tmp=(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	Result *=Tmp/Norm(Tmp);
      }
  for (int i=NbrParticlesPerBlock;i<NbrParticles;++i)
    for (int j=NbrParticlesPerBlock;j<i;++j)
      {
	Tmp=(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	Result *=Tmp/Norm(Tmp);
      }

  return Result;
}

// query weights of individual blocks
// weights = array reference to values of weighting factors
void MRBlockSamplingFunction::GetBlockWeights(double *weights)
{
  for (int i=0; i<NbrBlocks; ++i)
    weights[i]=BlockWeights[i];
}

// get the ratio of the block with nbrPermute particles exchanged, and the original setting
// nbrPermute = number of particles to exchange between blocks
Complex GetBlockWeight(int nbrBlock)
{
  int * Permutation=NULL;
  int NbrParticlesPerBlock=2;
  int NbrParticles=4;
  if (nbrBlock%NbrParticlesPerBlock==0)
    //return this->GetFunctionValue();
    {}
  for (int i=0; i<nbrBlock; ++i)
    {
      Permutation[i]=NbrParticlesPerBlock+i;
      Permutation[NbrParticlesPerBlock+i]=i;
    }
  for (int i=nbrBlock; i<NbrParticlesPerBlock; ++i)
    Permutation[i]=i;
  for (int i=NbrParticlesPerBlock+nbrBlock; i<NbrParticles; ++i)
    Permutation[i]=i;

  return Complex(0.0);
}


// get the Monte Carlo amplitude for the requested block with nbrPermute particles exchanged
// nbrBlock = nbrPermute = number of particles to exchange between blocks
// amplitude = reference of return argument
void MRBlockSamplingFunction::GetBlockAmplitude(int nbrBlock, Complex &amplitude)
{
  this->EvaluateTable();
  double AbsJ=Norm(JastrowFactor);

  amplitude=1.0;
  for (int i=1;i<NbrParticlesPerBlock;++i)
    for (int j=0;j<i;++j)
      amplitude *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
  for (int i=NbrParticlesPerBlock;i<NbrParticles;++i)
    for (int j=NbrParticlesPerBlock;j<i;++j)
      amplitude *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
  amplitude*=amplitude;
  if (nbrBlock == 0)
    {
      amplitude.Re = BlockWeights[0] * Norm(amplitude) * AbsJ;
      amplitude.Im = 0.0;
      return;
    }

  //complex phase of Block 0
  Complex ArgPhi0 = amplitude;
  ArgPhi0 *= 1.0/Norm(ArgPhi0);

  amplitude=1.0;
  // Block I
  int ToSwap = nbrBlock-1;
  for (int i=ToSwap+2;i<NbrParticlesPerBlock;++i)
    for (int j=ToSwap+1;j<i;++j)
      amplitude *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
  for (int i=NbrParticlesPerBlock+1;i<NbrParticlesPerBlock+ToSwap+1;++i)
    for (int j=NbrParticlesPerBlock;j<i;++j)
      amplitude *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
  for (int i=NbrParticlesPerBlock;i<NbrParticlesPerBlock+ToSwap+1;++i)
    for (int j=ToSwap+1;j<NbrParticlesPerBlock;++j)
      amplitude *= -this->ElementNorm*JastrowElements[SymIndex(j, i)];
  // Block II
  for (int i=NbrParticlesPerBlock+ToSwap+2;i<NbrParticles;++i)
    for (int j=NbrParticlesPerBlock+ToSwap+1;j<i;++j)
      amplitude *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
  for (int i=1;i<ToSwap+1;++i)
    for (int j=0;j<i;++j)
      amplitude *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
  for (int i=0;i<ToSwap+1;++i)
    for (int j=NbrParticlesPerBlock+ToSwap+1;j<NbrParticles;++j)
      amplitude *= this->ElementNorm*JastrowElements[SymIndex(i, j)];

  amplitude *= amplitude;

  amplitude = BlockWeights[nbrBlock] * Conj(amplitude)*ArgPhi0*AbsJ;


#ifdef TESTING_MRBLOCKSAMPLINGFUNCTION
  // testing using the explicit permutations
  Complex NewAmplitude=1.0;
  int Sign=0,Index;
  for (int i=1;i<NbrParticlesPerBlock;++i)
    for (int j=0;j<i;++j)
      {
	Index = AsymIndex(BlockPermutations[nbrBlock][j], BlockPermutations[nbrBlock][i],Sign);
	NewAmplitude *= Sign*this->ElementNorm*JastrowElements[Index];
      }
  for (int i=NbrParticlesPerBlock;i<NbrParticles;++i)
    for (int j=NbrParticlesPerBlock;j<i;++j)
      {
	Index = AsymIndex(BlockPermutations[nbrBlock][j], BlockPermutations[nbrBlock][i],Sign);
	NewAmplitude *= Sign*this->ElementNorm*JastrowElements[Index];
      }
  NewAmplitude*=NewAmplitude;
  NewAmplitude = BlockWeights[nbrBlock] * Conj(NewAmplitude)*ArgPhi0*AbsJ;

  if (Norm(NewAmplitude-amplitude)>1e-10)
    cout << "serious discrepancy B"<<nbrBlock<<": "<<amplitude<<" vs "<<NewAmplitude<<endl;
#endif
  
}

// get the Monte Carlo amplitude for the requested block with nbrPermute particles exchanged
// amplitudes = pointer to array of return arguments
void MRBlockSamplingFunction::GetAllBlockAmplitudes(Complex *amplitudes)
{
  this->EvaluateTable();
  double AbsJ=Norm(JastrowFactor);

  // first step: evaluate block 0 (without permutations)
  amplitudes[0]=1.0;
  for (int i=1;i<NbrParticlesPerBlock;++i)
    for (int j=0;j<i;++j)
      amplitudes[0] *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
  for (int i=NbrParticlesPerBlock;i<NbrParticles;++i)
    for (int j=NbrParticlesPerBlock;j<i;++j)
      amplitudes[0] *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
  //complex phase of Block 0
  Complex ArgPhi0 = amplitudes[0]*amplitudes[0];
  ArgPhi0 *= 1.0/Norm(ArgPhi0);
  //cout <<"AbsJ="<<AbsJ<<", amplitudes[0]="<< amplitudes[0] << endl;
  
  bool ProblemWithPrevious=false;
  // to proceed: evaluate blocks one by one from previous - k: position of last particle to be swapped
  for (int k=0; k<NbrBlocks-1; ++k)
    {
      // test if there is any issue with the block - then recalculate
      if (CriticalFlag && ( ProblemWithPrevious || (this->CriticalParticles[k>>USHIFT]&(0x1ul<<(k%ULENGTH)))
			    ||(this->CriticalParticles[(k+NbrParticlesPerBlock)>>USHIFT]&(0x1ul<<((k+NbrParticlesPerBlock)%ULENGTH)))))
	{
	  ProblemWithPrevious=true;
	  amplitudes[k+1]=1.0;
	  // Block I
	  for (int i=k+2;i<NbrParticlesPerBlock;++i)
	    for (int j=k+1;j<i;++j)
	      amplitudes[k+1] *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
	  for (int i=NbrParticlesPerBlock+1;i<NbrParticlesPerBlock+k+1;++i)
	    for (int j=NbrParticlesPerBlock;j<i;++j)
	      amplitudes[k+1] *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
	  for (int i=NbrParticlesPerBlock;i<NbrParticlesPerBlock+k+1;++i)
	    for (int j=k+1;j<NbrParticlesPerBlock;++j)
	      amplitudes[k+1] *= -this->ElementNorm*JastrowElements[SymIndex(j, i)];
	  // Block II
	  for (int i=NbrParticlesPerBlock+k+2;i<NbrParticles;++i)
	    for (int j=NbrParticlesPerBlock+k+1;j<i;++j)
	      amplitudes[k+1] *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
	  for (int i=1;i<k+1;++i)
	    for (int j=0;j<i;++j)
	      amplitudes[k+1] *= this->ElementNorm*JastrowElements[SymIndex(j, i)];
	  for (int i=0;i<k+1;++i)
	    for (int j=NbrParticlesPerBlock+k+1;j<NbrParticles;++j)
	      amplitudes[k+1] *= this->ElementNorm*JastrowElements[SymIndex(i, j)];
	}
      else
	{
	  Complex Num(1.0), Den(1.0);
	  // terms to remove from Block I
	  // terms of z_k, with z_(NbrParticlesPerBlock+i), i<k
	  for(int i=0; i<k; ++i)
	    Den *= -JastrowElements[SymIndex(k,NbrParticlesPerBlock+i)];
	  // terms of z_k, with z_i, i>k
	  for(int i=k+1; i<NbrParticlesPerBlock; ++i)
	    Den *= JastrowElements[SymIndex(k,i)];

	  // terms to remove from Block II
	  // terms of z_i, with z_(NbrParticlesPerBlock+k), i<k
	  for(int i=0; i<k; ++i)
	    Den *= JastrowElements[SymIndex(i,NbrParticlesPerBlock+k)];
	  // terms of z_k, with z_i, i>k
	  for(int i=NbrParticlesPerBlock+k+1; i<NbrParticles; ++i)
	    Den *= JastrowElements[SymIndex(k,i)];

	  // terms to add to Block I
	  // terms of z_(NbrParticlesPerBlock+k), with z_i, i>k
	  for(int i=k+1; i<NbrParticlesPerBlock; ++i)
	    Num *= -JastrowElements[SymIndex(i,NbrParticlesPerBlock+k)];
	  // terms of z_(NbrParticlesPerBlock+k), with z_(NbrParticlesPerBlock+i), i<k
	  for(int i=0; i<k; ++i)
	    Num *= JastrowElements[SymIndex(NbrParticlesPerBlock+i,NbrParticlesPerBlock+k)];
	  // terms to add to Block II
	  // terms of z_(k), with z_(i), i<k
	  for(int i=0; i<k; ++i)
	    Num *= JastrowElements[SymIndex(i,k)];
	  // terms of z_k, with z_(NbrParticlesPerBlock+i), i>k
	  for(int i=k+1; i<NbrParticlesPerBlock; ++i)
	    Num *= -JastrowElements[SymIndex(k,NbrParticlesPerBlock+i)];
	  
	  amplitudes[k+1]=amplitudes[k]*Num/Den;

	  ProblemWithPrevious=false;
	}
    }

  for (int i=0; i<NbrBlocks; ++i)
    amplitudes[i] = Conj(amplitudes[i]*amplitudes[i]) * ArgPhi0 * (AbsJ * BlockWeights[i]);

#ifdef TESTING_MRBLOCKSAMPLINGFUNCTION
  Complex TmpC;
  for (int i=0; i<NbrBlocks; ++i)
    {
      this->GetBlockAmplitude(i, TmpC);
      //if (Norm(TmpC-amplitudes[i])>1e-10)
	//cout << "discrepancy B"<<i<<": "<<TmpC<<" vs "<<amplitudes[i]<<endl;
      amplitudes[i] = TmpC;
    }

#endif  
}


// query permutation of particles applied to given block
// nbrBlock = block index
// permutations = pointer to array to be filled with return argument
void MRBlockSamplingFunction::GetBlockPermutations(int nbrBlock, int* permutations)
{
  for (int i=0; i<NbrParticles; ++i)
    permutations[i]=BlockPermutations[nbrBlock%NbrBlocks][i];
}

// query permutation of particles applied to given block
// nbrBlock = block index
// permutations = pointer to array to be filled with return argument
void MRBlockSamplingFunction::GetAllBlockPermutations(int** permutations)
{
  for (int b=0; b<NbrBlocks; ++b)
    for (int i=0; i<NbrParticles; ++i)
      permutations[b][i]=BlockPermutations[b][i];
}


// precalculate Jastrow factor elements
void MRBlockSamplingFunction::EvaluateTable()
{
  for (int i=0; i<NbrParticles>>USHIFT; ++i)
    this->CriticalParticles[i] = 0x0ul;
  this->CriticalFlag=false;
  JastrowFactor=1.0;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      for (int j = 0; j < i; ++j)
	{
	  int Index = SymIndex(j, i);
	  
	  JastrowElements[Index] = ((this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]));
	  //check if critical value was encountered
	  if (SqrNorm(JastrowElements[Index])<this->SqrCriticalDistance) 
	    {
	      this->CriticalFlag=true;
	      this->CriticalParticles[i>>USHIFT]|=(0x1ul<<(i%ULENGTH));
	      this->CriticalParticles[j>>USHIFT]|=(0x1ul<<(j%ULENGTH));
	    }
	  JastrowFactor *= JastrowElements[Index]*ElementNorm;
	}
    }
}
