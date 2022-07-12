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


#include "BlockMonteCarloOnSphereAlgorithm.h"
#include "TrivialSamplingFunction.h"
#include "ParticleOnSphereCollection.h"
#include "ParticleOnSphereCollectionSouthPole.h"
#include "ParticleOnSphereCollectionGauged.h"
#include "Options/Options.h"


#include <iostream>
using std::cout;
using std::endl;

// default constructor
BlockMonteCarloOnSphereAlgorithm::BlockMonteCarloOnSphereAlgorithm()
{
  this->NbrParticles=0;
}

// set up for basic monte-carlo scheme
// nbrParticles = number of particles in system
// targetFunction = block sampling function describing the target wavefunction
// manager = pointer to option manager
// maxNbrObservables = maximum number of observables to be assigned
BlockMonteCarloOnSphereAlgorithm::BlockMonteCarloOnSphereAlgorithm(int nbrParticles,
		    AbstractMCBlockSamplingFunctionOnSphere *targetFunction,
		    OptionManager *manager, int maxNbrObservables)
{
  if (targetFunction==0)
    {
      cout << "Error: no target function in BlockMonteCarloOnSphereAlgorithm" << endl;
      exit(1);
    }
  this->TargetFunction = targetFunction;
  this->NbrParticles=nbrParticles;
  this->NbrBlocks = this->TargetFunction->GetNbrBlocks();

  this->Options=manager;  
  long Seed = Options->GetInteger("randomSeed");
  if (manager->GetBoolean("no-southpole"))
    this->System=new ParticleOnSphereCollection(this->NbrParticles, Seed);
  else
    {
      if (manager->GetBoolean("gauged"))
	this->System=new ParticleOnSphereCollectionGauged(this->NbrParticles, Seed);
      else
	this->System=new ParticleOnSphereCollectionSouthPole(this->NbrParticles, Seed);
    }

  this->TargetFunction->RegisterSystem(this->System);
  this->TargetFunction->AdaptAverageMCNorm(Options->GetInteger("thermalize")); // this also relaxes the particle positions in System
  // renormalize wavefunction, as this has led to problems
  //this->NormalizePsi();
  this->MaxNbrObservables=maxNbrObservables;
  this->Observables= new AbstractObservable*[MaxNbrObservables];
  this->Frequencies= new int[MaxNbrObservables];
  this->NbrObservables=0;
  this->MaxNbrBlockObservables=maxNbrObservables;
  this->BlockObservables= new AbstractBlockObservable*[MaxNbrObservables];
  this->BlockFrequencies= new int[MaxNbrObservables];
  this->NbrBlockObservables=0;
}
  
// destructor
BlockMonteCarloOnSphereAlgorithm::~BlockMonteCarloOnSphereAlgorithm()
{
  if (this->NbrParticles!=0)
    {
      delete this->System;      
      delete [] this->Frequencies;
      for (int i=0; i<NbrObservables; ++i)
	delete this->Observables[i];
      delete [] this->Observables;
      delete [] this->BlockFrequencies;
      for (int i=0; i<NbrBlockObservables; ++i)
	delete this->BlockObservables[i];
      delete [] this->BlockObservables;
    }
}

// add a regular observable
// O = pointer to the observable to be added
// frequency = integer value indicating on which microsteps observations are being made
void BlockMonteCarloOnSphereAlgorithm::AddObservable(AbstractObservable *O, int frequency)
{
  if (this->NbrObservables<this->MaxNbrObservables)
    {
      this->Frequencies[this->NbrObservables]=frequency;
      this->Observables[this->NbrObservables++]=O;
    }
  else
    {
      AbstractObservable **newObservables = new AbstractObservable*[this->MaxNbrObservables+10];
      int *newFrequencies = new int[this->MaxNbrObservables+10];
      for (int i=0; i<NbrObservables; ++i)
	{
	  newObservables[i]=this->Observables[i];
	  newFrequencies[i]=this->Frequencies[i];
	}
      newFrequencies[this->NbrObservables]=frequency;
      newObservables[NbrObservables++]=O;
      delete [] this->Observables;
      this->Observables = newObservables;
      this->Frequencies=newFrequencies;
    }
  O->SetParticleCollection(this->System);
}

// add a block observable for detailed views of the contribution of individual terms in the symmetric expansion
// O = pointer to the observable to be added
// frequency = integer value indicating on which microsteps observations are being made
void BlockMonteCarloOnSphereAlgorithm::AddBlockObservable(AbstractBlockObservable *O, int frequency)
{
  if (this->NbrBlockObservables<this->MaxNbrObservables)
    {
      this->BlockFrequencies[this->NbrBlockObservables]=frequency;
      this->BlockObservables[this->NbrBlockObservables++]=O;
    }
  else
    {
      AbstractBlockObservable **newObservables = new AbstractBlockObservable*[this->MaxNbrBlockObservables+10];
      int *newFrequencies = new int[this->MaxNbrBlockObservables+10];
      for (int i=0; i<NbrBlockObservables; ++i)
	{
	  newObservables[i]=this->BlockObservables[i];
	  newFrequencies[i]=this->BlockFrequencies[i];
	}
      newFrequencies[this->NbrBlockObservables]=frequency;
      newObservables[NbrBlockObservables++]=O;
      delete [] this->BlockObservables;
      this->BlockObservables = newObservables;
      this->BlockFrequencies=newFrequencies;
    }
  O->SetParticleCollection(this->System);
}

// thermalize system with a number of microsteps
// time = number of microsteps
// startFromRandom = flag indicating if we want to restart from a random configuration
void BlockMonteCarloOnSphereAlgorithm::Thermalize(int time, bool startFromRandom)
{
  if (startFromRandom)
    System->Randomize();
  this->PerformMicroSteps(time);
}

// renormalize wavefunction
// time = number of points to average
void BlockMonteCarloOnSphereAlgorithm::NormalizePsi(int time)
{
  double SumSqrPsiValues=0.0;
  for (int t=0; t<time; ++t)
    {
      PerformMicroSteps(1);
      SumSqrPsiValues+=SqrNorm(this->TargetFunction->GetFunctionValue());
    }
  this->TargetFunction->ScaleByFactor(sqrt((double)time/SumSqrPsiValues));  
}

// run simulation
// all options are included via the AddOptionGroup method
void BlockMonteCarloOnSphereAlgorithm::Simulate(ostream &Output)
{
  this->NbrAcceptedMoves = 0l;
  this->NbrAttemptedMoves = 0l;
  int NbrSteps = this->Options->GetInteger("nbr-iter");
  int DensityOfSamples = this->Options->GetInteger("sample-density");
  if (DensityOfSamples < 0)
    DensityOfSamples = this->NbrParticles;
  int NbrDisplay = this->Options->GetInteger("nbr-display");
  if (NbrDisplay==0) NbrDisplay=1;
  int DisplaySteps = NbrSteps / NbrDisplay;
  Complex *BlockAmplitudes = new Complex[NbrBlocks+1];
  Complex *TotalAmplitude = &BlockAmplitudes[NbrBlocks];
  Output <<"Step";
  for (int i=0; i<NbrObservables; ++i)
    {
      Output << "\t";
      Observables[i]->PrintLegend(Output);
    }
  Output << endl;
  int s=0;
  //cout << "NbrBlocks="<<NbrBlocks<<endl;
  for (int d=0; d<NbrDisplay; ++d)
    {
      for (/* s */; s<(d+1)*DisplaySteps; ++s)
	{
	  this->PerformMicroSteps(DensityOfSamples);
	  this->TargetFunction->GetAllBlockAmplitudes(BlockAmplitudes);
	  //if (s%500) cout << "BlockAmplitudes[0]="<<BlockAmplitudes[0]<<endl;
	  *TotalAmplitude=0.0;
	  for (int i=0; i<NbrBlocks; ++i)
	    *TotalAmplitude += BlockAmplitudes[i];
	  for (int i=0; i<NbrObservables; ++i)
	    if (s%Frequencies[i]==0) Observables[i]->RecordValue(TotalAmplitude->Re);
	  for (int i=0; i<NbrBlockObservables; ++i)
	    if (s%BlockFrequencies[i]==0) BlockObservables[i]->RecordValue(BlockAmplitudes);
	}      
      Output <<(d+1)*DisplaySteps;
      for (int i=0; i<NbrObservables; ++i)
	{
	  Output << "\t";
	  Observables[i]->PrintStatus(Output);	  
	}
      for (int i=0; i<NbrBlockObservables; ++i)
	{
	  Output << "\t";
	  BlockObservables[i]->PrintStatus(Output);	  
	}
      Output << endl;
    }
  cout << "====== accepted moves: "<< (double)this->NbrAcceptedMoves / (double)this->NbrAttemptedMoves *100.0
       <<"% ======"<<endl;
  delete [] BlockAmplitudes;
}


// perform a number of Monte-Carlo microsteps
// nbrSteps = number of steps
//
void BlockMonteCarloOnSphereAlgorithm::PerformMicroSteps(int nbrSteps)
{
  double acceptanceProbability;
  for (int i = 0; i < nbrSteps; ++i)
    {
      System->Move();      
      ++this->NbrAttemptedMoves;
      acceptanceProbability = TargetFunction->GetTransitionRatio();
      if ((acceptanceProbability < 1.0) &&  (System->GetRandomNumber() > acceptanceProbability))
	{
	  System->RestoreMove();
	}
      else
	{
	  this->TargetFunction->AcceptedMove();
	  ++this->NbrAcceptedMoves;
	}
    }
}


// add an option group containing all options related to the wave functions
//
// manager = pointer to the option manager
void BlockMonteCarloOnSphereAlgorithm::AddOptionGroup(OptionManager* manager)
{
 OptionGroup* MCGroup  = new OptionGroup ("Monte-Carlo options");
  (*(manager)) += MCGroup;  
  (*MCGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MCGroup) += new SingleIntegerOption  ('d', "sample-density", "number of Monte Carlo microsteps between iterations (-1=NbrParticles)", -1);
  (*MCGroup) += new SingleIntegerOption  ('\n',
					  "thermalize", "number of microsteps used for thermalization", 2500);  
  (*MCGroup) += new SingleIntegerOption  ('\n', "nbr-display",
					  "number of intermediate results displayed", 10);
  (*MCGroup) += new SingleIntegerOption  ('\n', "randomSeed", "random seed for internal random number generator", -1);
    
}

// write measurements of all observables to a given stream
// str = output stream
//
void BlockMonteCarloOnSphereAlgorithm::WriteObservations(ostream &str)
{
  for (int i=0; i<NbrObservables; ++i)
    {
      if (Observables[i]->IncludeInPrint())
	Observables[i]->WriteDataFile(str);
    }
  for (int i=0; i<NbrBlockObservables; ++i)
    {
      if (BlockObservables[i]->IncludeInPrint())
	BlockObservables[i]->WriteDataFile(str);
    }

}
