////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of main task for qhe overlap on sphere                //
//                                                                            //
//                        last modification : 21/04/2005                      //
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
#include "MainTask/QHEOnSphereOverlapMainTask.h"

#include "MathTools/RandomNumber/AbstractRandomNumberGenerator.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"

#include "Vector/RealVector.h"

#include "Options/OptionManager.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// constructor
//  
// options = pointer to the options managers containing all running options
// referenceWaveFunction = pointer to the reference wave function
// testWaveFunction = pointer to the test wave function
// nbrParticles = number of particles
// randomNumberGenerator = pointer to the random number generator

QHEOnSphereOverlapMainTask::QHEOnSphereOverlapMainTask(OptionManager* options, Abstract1DComplexFunction* referenceWaveFunction, 
						       Abstract1DComplexFunction* testWaveFunction, int nbrParticles, 
						       AbstractRandomNumberGenerator* randomNumberGenerator)
{
  this->ReferenceWaveFunction = referenceWaveFunction;
  this->TestWaveFunction = testWaveFunction;
  this->NbrParticles = nbrParticles;
  this->RandomNumberGenerator = randomNumberGenerator;

  this->ShowDetailFlag = ((BooleanOption*) (*options)["show-details"])->GetBoolean();
  this->DisplayStep = ((SingleIntegerOption*) (*options)["display-step"])->GetInteger();
  this->NbrIterations = ((SingleIntegerOption*) (*options)["nbr-iter"])->GetInteger();
}  
 
// destructor
//  

QHEOnSphereOverlapMainTask::~QHEOnSphereOverlapMainTask()
{
}
  
// execute the main task
// 
// return value = 0 if no error occurs, else return error code

int QHEOnSphereOverlapMainTask::ExecuteMainTask()
{
  RealVector Location(2 * this->NbrParticles, true);
  double Factor = 1.0;
  for (int j = 0; j < this->NbrParticles; ++j)
    {
      Factor *= 4.0 * M_PI;
    }
  Complex Overlap;
  Complex ErrorOverlap;
  double Normalization = 0.0;
  double ErrorNormalization = 0.0;
  double NormalizationExact = 0.0;
  double ErrorNormalizationExact = 0.0;
  Complex Tmp;
  Complex Tmp3;
  double Tmp2;
  double Tmp2bis;
  int NextCoordinates = 0;
  for (int j = 0; j < this->NbrParticles; ++j)
    {
      Location[j << 1] = acos (1.0- (2.0 * this->RandomNumberGenerator->GetRealRandomNumber()));
      Location[1 + (j << 1)] = 2.0 * M_PI * this->RandomNumberGenerator->GetRealRandomNumber();
    }
  Tmp = (*TestWaveFunction)(Location);
  double PreviousProbabilities = Norm(Tmp);
  double CurrentProbabilities = PreviousProbabilities;
  double TotalProbability = PreviousProbabilities;
  for (int i = 0; i < this->NbrIterations; ++i)
    {
      double PreviousCoordinates1 = Location[NextCoordinates << 1];
      double PreviousCoordinates2 = Location[1 + (NextCoordinates << 1)];
      Location[NextCoordinates << 1] = acos (1.0- (2.0 * this->RandomNumberGenerator->GetRealRandomNumber()));	  
      Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * this->RandomNumberGenerator->GetRealRandomNumber();
      Complex TmpMetropolis = (*(this->TestWaveFunction))(Location);
      CurrentProbabilities = Norm(TmpMetropolis);
      if ((CurrentProbabilities > PreviousProbabilities) || ((this->RandomNumberGenerator->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	{
	  PreviousProbabilities = CurrentProbabilities;
	  Tmp = TmpMetropolis;
	}
      else
	{
	  Location[NextCoordinates << 1] = PreviousCoordinates1;
	  Location[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
	  CurrentProbabilities = PreviousProbabilities;
	}
      TotalProbability += CurrentProbabilities;
      NextCoordinates = (int) (((double) this->NbrParticles) * this->RandomNumberGenerator->GetRealRandomNumber());
      if (NextCoordinates == this->NbrParticles)
	--NextCoordinates;

      Complex ValueExact = (*(this->ReferenceWaveFunction))(Location);
      Tmp2 = (Tmp.Re * Tmp.Re) + (Tmp.Im * Tmp.Im);
      Tmp2bis = (ValueExact.Re * ValueExact.Re) + (ValueExact.Im * ValueExact.Im);
      Tmp3 = (Conj(Tmp) * ValueExact);
      Tmp2 /= CurrentProbabilities;
      Tmp3 /= CurrentProbabilities;      
      Tmp2bis /= CurrentProbabilities;  
      Overlap += Tmp3;
      ErrorOverlap.Re += Tmp3.Re * Tmp3.Re;
      ErrorOverlap.Im += Tmp3.Im * Tmp3.Im;
      Normalization += Tmp2;
      ErrorNormalization += Tmp2 * Tmp2;
      NormalizationExact += Tmp2bis;
      ErrorNormalizationExact += Tmp2bis * Tmp2bis;
      if ((i > 0) && ((i % this->DisplayStep) == 0))
	{
	  cout << " i = " << i << endl;
	  Complex Tmp4 = Overlap / ((double) i);
	  Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
			sqrt( ((ErrorOverlap.Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
	  double Tmp6 = Normalization  / ((double) i);
	  double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	  double Tmp8 = NormalizationExact  / ((double) i);
	  double Tmp9 = sqrt( ((ErrorNormalizationExact / ((double) i))  -  (Tmp8 * Tmp8)) / ((double) i) );	  

	  if (this->ShowDetailFlag == true)
	    {
	      cout << Tmp4;
	      cout << " +/- " << Tmp5 << endl;
	      cout << Tmp6;
	      cout << " +/- " << Tmp7 << endl;	  
	      cout << Tmp8;
	      cout << " +/- " << Tmp9 << endl;	  
	    }
	  Tmp5.Re /= Tmp4.Re;
	  Tmp5.Im /= Tmp4.Im;
	  Tmp5.Re = fabs(Tmp5.Re);
	  Tmp5.Im = fabs(Tmp5.Im);
	  Tmp5.Re += (Tmp7 / Tmp6);
	  Tmp5.Im += (Tmp7 / Tmp6);
	  Tmp5.Re += (Tmp9 / Tmp8);
	  Tmp5.Im += (Tmp9 / Tmp8);
	  Tmp4 /= sqrt(Tmp6 * Tmp8);	  
	  Tmp5.Re *= Tmp4.Re;
	  Tmp5.Im *= Tmp4.Im;
	  cout << Tmp4 << " +/- " << Tmp5 << endl;
	  cout << "-----------------------------------------------" << endl;
	}
    } 
  cout << " final results :" << endl;
  Complex Tmp4 = Overlap / ((double) this->NbrIterations);
  cout << Tmp4;
  Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) this->NbrIterations)) - (Tmp4.Re * Tmp4.Re)) / ((double) this->NbrIterations) ),
		sqrt( ((ErrorOverlap.Im / ((double) this->NbrIterations)) - (Tmp4.Im * Tmp4.Im)) / ((double) this->NbrIterations) ));
  cout << " +/- " << Tmp5 << endl;
  double Tmp6 = Normalization  / ((double) this->NbrIterations);
  cout << Tmp6;
  double Tmp7 = sqrt( ((ErrorNormalization / ((double) this->NbrIterations))  -  (Tmp6 * Tmp6)) / ((double) this->NbrIterations) );	  
  cout << " +/- " << Tmp7 << endl;	  
  double Tmp8 = NormalizationExact  / ((double) this->NbrIterations);
  cout << Tmp8;
  double Tmp9 = sqrt( ((ErrorNormalizationExact / ((double) this->NbrIterations))  -  (Tmp8 * Tmp8)) / ((double) this->NbrIterations) );	  
  cout << " +/- " << Tmp9 << endl;	  
  
  Tmp5.Re /= Tmp4.Re;
  Tmp5.Im /= Tmp4.Im;
  Tmp5.Re = fabs(Tmp5.Re);
  Tmp5.Im = fabs(Tmp5.Im);
  Tmp5.Re += (Tmp7 / Tmp6);
  Tmp5.Im += (Tmp7 / Tmp6);
  Tmp5.Re += (Tmp9 / Tmp8);
  Tmp5.Im += (Tmp9 / Tmp8);
  Tmp4 /= sqrt(Tmp6 * Tmp8);	  
  Tmp5.Re *= Tmp4.Re;
  Tmp5.Im *= Tmp4.Im;
  cout << Tmp4 << " +/- " << Tmp5 << endl;
  cout << "-----------------------------------------------" << endl;
  return 0;
}

