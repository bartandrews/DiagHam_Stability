/***************************************************************************
 *          ComplexVectorObservable.h
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

#ifndef COMPLEXVECTOROBSERVABLE_H
#define COMPLEXVECTOROBSERVABLE_H

#include "config.h"

#include <vector>
#include <string>
#include <ostream>
#include "MathTools/Complex.h"

using std::string;
using std::ostream;

typedef Complex WCVOIndividualType;
typedef Complex* WCVOType;
#define DEF_OBS_PER_BIN_CVO 2

class ComplexVectorObservable
{
 public:
  ComplexVectorObservable (int length, int maxBins=256, int minObsPerBin=DEF_OBS_PER_BIN_CVO);
  ~ComplexVectorObservable();
  
  void operator<<(const WCVOType x); // makes observation for all observables (vector valued)
  void operator<<(const WCVOIndividualType x); // makes observationon of the same scalar observable value
  void Observe(const WCVOIndividualType x); // observation with different weight on the same observable
  void Observe(const WCVOType allX);
  
  int MaxBinNumber() const { return MaxBinNum;}
  int BinNumber() const;
  int FilledBinNumber() const;
  unsigned long Measurements(); 
  void Rescale(double factor);
  void Rescale(int field, double factor);
  WCVOType Average();  // returns the average of the measurements
  WCVOIndividualType Average(int field);  // returns the average of the measurements
  WCVOIndividualType PresentBinAverage(int);  // returns the average of the measurements in the present bin
  WCVOIndividualType PresentBinAverage2(int);  // returns the respective average of the measurements of squares
  double Variance(int field);
  double* Variance(); // returns the variance of the measurements
  double* VarianceOfBins();
  double VarianceOfBins(int field);
  double* ErrorEstimate();
  double ErrorEstimate(int field);
  // variance and error for real part
  double VarianceRe(int field);
  double* VarianceRe(); // returns the variance of the measurements
  double* VarianceOfBinsRe();
  double VarianceOfBinsRe(int field);
  double* ErrorEstimateRe();
  double ErrorEstimateRe(int field);
  // variance and error for imaginary part
  double VarianceIm(int field);
  double* VarianceIm(); // returns the variance of the measurements
  double* VarianceOfBinsIm();
  double VarianceOfBinsIm(int field);
  double* ErrorEstimateIm();
  double ErrorEstimateIm(int field);
  
  void SetBinNumber(int BinNum); // set maximal bin number

  void SetName(string &newName);
  void SetName(char *newName);
  
  void SetFieldName(int field, string &newName);
  void SetFieldName(int field, char *newName);

  int GetLength(){return this->NumFields;}

  void Rebin(unsigned int sequence);
  
  
  const WCVOIndividualType& BinValue(int entry, int i) const  { return (Values[entry])[i]; }
  double BinValue2(int entry,int i) const { return (Values2Re[entry][i]+Values2Im[entry][i]) ;}
  const WCVOIndividualType& HistoryValue(int entry, int i) const  { return (HistoryAverages[entry])[i]; }
  const double& HistoryValue2(int entry, int i) const { return (HistoryBinVariances[entry])[i];}
  
  friend ostream& operator << (ostream& str, ComplexVectorObservable &O);
  
 protected:
  string Name;
  string *Names;
  int NumFields; // number of observables in string
  int StartBins;        //arguments of constructor stored for reset
  int StartObsPerBin; //dito
  unsigned long BinSize;       // number of measurements per bin
  unsigned int MaxBinNum;    // maximum number of bins
  unsigned long BinEntries;    // number of measurements in last bin
  unsigned long PreviousEntries;
  std::vector<WCVOIndividualType>* Values;  // bin values
  std::vector<double>* Values2Re; // bin values of weighted squares
  std::vector<double>* Values2Im; // bin values of weighted squares
  std::vector<WCVOIndividualType>* HistoryAverages;
  std::vector<double>* HistoryBinVariances;
};

inline void ComplexVectorObservable::operator<<(const WCVOType x)
{
  this->Observe(x);
}

// makes observation with default weight 1.0 on the same observable
inline void ComplexVectorObservable::operator<<(const WCVOIndividualType x)
{
  this->Observe(x);
}


#endif // WEIGHTEDREALVECOBS_H


