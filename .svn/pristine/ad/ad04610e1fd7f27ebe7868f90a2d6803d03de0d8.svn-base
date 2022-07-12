/***************************************************************************
 *          WeightedComplexVectorObservable.h
 *
 *  Sat Dec  17 12:30:16 2005
 *  Copyright  2005  Gunnar Möller
 *  moller@lptms.u-psud.fr
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef WEIGHTEDCOMPLEXVECTOROBSERVABLE_H
#define WEIGHTEDCOMPLEXVECTOROBSERVABLE_H

#include "config.h"

#include <vector>
#include <string>
#include <ostream>
#include "MathTools/Complex.h"

using std::string;
using std::ostream;

typedef Complex WCVOIndividualType;
typedef Complex* WCVOType;
typedef double WCVOIndividualWeightType;
typedef double* WCVOWeightType;

class WeightedComplexVectorObservable
{
 public:
  WeightedComplexVectorObservable (int length, int maxBins=256 , WCVOIndividualWeightType minObsPerBin=1.0);
  ~WeightedComplexVectorObservable();
  
  void operator<<(const WCVOType x); // makes observation with default weight 1.0 and different observables
  void Observe(const WCVOType x, WCVOWeightType w); // makes observation with different weight and different observables
  void operator<<(const WCVOIndividualType x); // makes observation with default weight 1.0 on the same observable
  void Observe(const WCVOIndividualType x, WCVOWeightType w=NULL); // observation with different weight on the same observable
  void Observe(const WCVOType allX, WCVOIndividualWeightType w=1.0);
  
  int MaxBinNumber() const { return MaxBinNum;}
  WCVOWeightType AverageBinSize();
  WCVOWeightType AverageBinWeight();
  WCVOIndividualWeightType AverageBinWeight(int field);
  int* BinNumber() const;
  int SingleBinNumber(int i) const;
  int* FilledBinNumber() const;
  unsigned Measurements(); 
  WCVOWeightType TotalWeight();
  WCVOIndividualWeightType TotalWeight(int field);
  WCVOWeightType AverageWeight();
  WCVOIndividualWeightType AverageWeight(int field);
  void Rescale(double factor);
  void Rescale(int field, double factor);
  WCVOType Average();  // returns the average of the measurements
  WCVOIndividualType Average(int field);  // returns the average of the measurements
  WCVOIndividualType PresentBinAverage(int);  // returns the average of the measurements in the present bin
  WCVOIndividualType PresentBinAverage2(int);  // returns the respective average of the measurements of squares
  double Variance(int field, WCVOIndividualWeightType typicalWgt=1.0);
  double* Variance(WCVOWeightType typical_wgt=NULL); // returns the variance of the measurements
  double* VarianceOfBins();
  double VarianceOfBins(int field);
  double* ErrorEstimate();
  double ErrorEstimate(int field);
  void SetBinNumber(int BinNum); // set maximal bin number

  void SetName(string &newName);
  void SetName(char *newName);
  
  void SetFieldName(int field, string &newName);
  void SetFieldName(int field, char *newName);

  int GetLength(){return this->NumFields;}

  void Rebin(int field,  unsigned int sequence);
  
  
  const WCVOIndividualType& BinValue(int entry, int i) const  { return (Values[entry])[i]; }
  const double& BinValue2(int entry,int i) const { return (Values2[entry])[i];}
  const WCVOIndividualType& HistoryValue(int entry, int i) const  { return (HistoryAverages[entry])[i]; }
  const double& HistoryValue2(int entry, int i) const { return (HistoryBinVariances[entry])[i];}
  
  friend ostream& operator << (ostream& str, WeightedComplexVectorObservable &O);
  
 protected:
  string Name;
  string *Names;
  int NumFields; // number of observables in string
  int StartBins;        //arguments of constructor stored for reset
  WCVOIndividualWeightType StartObsPerBin; //dito
  double *BinSize;       // number of measurements per bin
  unsigned int MaxBinNum;    // maximum number of bins
  double *BinWeight;
  unsigned *BinEntries;    // number of measurements in last bin
  double *PreviousWeight; // weight of all observations to this point
  unsigned *PreviousEntries;
  std::vector<WCVOIndividualType>* Values;  // bin values
  std::vector<double>* Values2; // bin values of weighted squares
  std::vector<int>* Observations;  // number of observations in bin
  std::vector<WCVOIndividualWeightType>* Weights;  // bin weights
  std::vector<WCVOIndividualType>* HistoryAverages;
  std::vector<double>* HistoryBinVariances;


};


#endif // WEIGHTEDREALVECOBS_H


