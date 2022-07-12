////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                      Copyright (C) 2007 Gunnar Möller                      //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 18/05/2007                      //
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


#ifndef WAVEFUNCTIONOVERLAPOPTIMIZER_H
#define WAVEFUNCTIONOVERLAPOPTIMIZER_H

#include "config.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "MCObservables/WeightedRealObservable.h"
#include "MCObservables/WeightedRealVectorObservable.h"
#include "MCObservables/WeightedComplexVectorObservable.h"
#include "MCObservables/MCHistoryRecord.h"
#include <fstream>
using std::ofstream;

class WaveFunctionOverlapOptimizer
{
 protected:  
  int NbrParticles;
  int NbrParameters;
  int EffectiveNbrParameters;
  int LimitSamples;
  Abstract1DComplexTrialFunction *TrialState;
  RealVector Positions;
  RealVector Gradient;
  ComplexVector ManyValues;
  double *InitialParameters;
  double InitialSqrOverlap;
  double **NewParameters;
  double *NormObservation;
  double *Differentials;
  double MinDifferential;
  double StepLength;
  double *LastOverlaps;
  Complex *OverlapObservation;
  int MaxPoints; 
  int LinearPoints;
  int CloudyPoints;
  int MaxParameters;
  double Precision;
  double NormExactWF;
  double ErrorNormExactWF;
  double OutlierLimit;
  WeightedRealVectorObservable *NormTrialObs;
  WeightedComplexVectorObservable *OverlapObs;
  MCHistoryRecord *History;
  bool LastParameterExcluded;
  double typicalSA;
  double typicalTV;
  double typicalWF;
  ofstream LogFile;

  NumRecRandomGenerator *Generator;
  
 public:

  // set up an optimizing algorithm based on a Newton type steepest descent procedure
  // trialState = trial wavefunction
  // historyFileName = precalculated record of values of the wavefunction for given particle positions
  // nbrParticles = number of particles
  // excludeLastParameter = exclude one parameters from optimization
  // cloudyPoints = number of random parameters values calculated simultaeusly along with steepest descent
  // linearPoints = number of points calculated at once along direction of steepest descent
  // limitSamples = upper limit to number of samples used from history-record
  // logFileName = name of an (optional) logfile
  WaveFunctionOverlapOptimizer( Abstract1DComplexFunction *trialState, char *historyFileName, int nbrParticles, bool excludeLastParameter = true, int linearPoints = 20, int cloudyPoints = 30, int limitSamples = 10000000, char* logFileName = NULL);
  ~WaveFunctionOverlapOptimizer();

  // launch optimization procedure
  // optimalParameters = initial, and upon return final variational parameters
  // Overlap = corresponding overlap, upon return
  double GetMaximumSqrOverlap(RealVector &optimalParameters, Complex &Overlap,
			      double toleranceFinal=1e-6, double toleranceIteration =0.01);

 private:

  void EvaluateTrialOverlaps();
  void DetermineGradientAndDifferentials(double *parameters);
  void CalculateLinearSeries(RealVector &startParameters, RealVector &stepDirection, RealVector &overlaps, RealMatrix &gradients);
  double GetRandomOffset(int parameter);
  
};

#endif // WAVEFUNCTIONOVERLAPOPTIMIZER_H
