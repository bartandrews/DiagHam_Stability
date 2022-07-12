/***************************************************************************
 *            WeightedRealObservable.h
 *
 *  Sat Dec  17 12:30:16 2005
 *  Author:  Gunnar Möller
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

#ifndef WEIGHTEDREALOBSERVABLE_H
#define WEIGHTEDREALOBSERVABLE_H

#include <vector>
#include <string>
#include <ostream>
using std::string;
using std::ostream;

#define WRO_DEF_MAX_BINS  256
#define WRO_DEF_OBS_PER_BIN  1.0

class WeightedRealObservable
{
 public:
  WeightedRealObservable (int MaxBins=WRO_DEF_MAX_BINS , double MinObsPerBin=WRO_DEF_OBS_PER_BIN);
  
  void operator<<(const double& x); // makes observation with default weight 1.0
  void Observe(const double& x, double w=1.0); // makes observation with different weight
  
  int MaxBinNumber() const { return MaxBinNum;}
  double AverageBinSize();
  double AverageBinWeight();
  int BinNumber() const;
  int FilledBinNumber() const;
  unsigned Measurements(); 
  double TotalWeight();
  double AverageWeight();
  
  void Rescale(double factor);
  
  double Average();  // returns the average of the measurements
  double PresentBinAverage();  // returns the average of the measurements in the present bin
  double PresentBinAverage2();  // returns the respective average of the measurements of squares 
  double Variance(double TypicalWgt=0.0); // returns the variance of the measurements
  double VarianceOfBins();
  double ErrorEstimate();
  void SetBinNumber(int binnum); // set maximal bin number

  void SetName(string &newName);
  void SetName(char *newName);
  void Rebin(unsigned int sequence);
  
  const double& BinValue(int i) const  { return Values[i]; }
  const double& BinValue2(int i) const { return Values2[i];}
  const double& HistoryValue(int i) const  { return HistoryAverages[i]; }
  const double& HistoryValue2(int i) const { return HistoryBinVariances[i];}

  void PrintInternal(); // some testing output...
  friend ostream& operator << (ostream& str, WeightedRealObservable &O);
  
 protected:
  string Name;
  int StartBins;        //arguments of constructor stored for reset
  double StartObsPerBin; //dito
  double BinSize;       // number of measurements per bin
  unsigned int MaxBinNum;    // maximum number of bins
  double BinWeight;
  unsigned BinEntries;    // number of measurements in last bin
  double PreviousWeight; // weight of all observations to this point
  unsigned PreviousEntries;
  std::vector<double> Values;  // bin values
  std::vector<double> Values2; // bin values of weighted squares
  std::vector<int> Observations;  // number of observations in bin
  std::vector<double> Weights;  // bin weights

  std::vector<double> HistoryAverages;
  std::vector<double> HistoryBinVariances;


};


#endif // WEIGHTEDREALOBS_H


