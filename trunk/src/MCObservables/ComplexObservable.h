////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//           class implementing a paired CF wave function on the sphere          //
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


#ifndef COMPLEXOBSERVABLE_H
#define COMPLEXOBSERVABLE_H

#include "config.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include "MathTools/Complex.h"


#define LARGE_N 256  // defines the bin size from which a partly filled bin is used in Variance
#define DEF_OBS_PER_BIN 2
#define DEF_MAX_BINS 256

using std::ostream;
using std::string;
using std::cout;
using std::endl;
using std::abs;


class ComplexObservable{
  public:

  // constructor
  ComplexObservable(int Maxbins=DEF_MAX_BINS , int minObsPerBin=DEF_OBS_PER_BIN);
  
  void operator<<(const Complex& x);

  friend ostream& operator << (ostream& str, ComplexObservable &B);
  
  void SetName(string &newName);
  void SetName(char *newName);
  
  int MaxBinNumber() const { return MaxBinNum;}
  int GetBinSize() const { return BinSize;}
  int BinNumber() const;
  int FilledBinNumber() const;
  int Measurements(); 

  void Rescale(double factor);
  
  Complex Average();  // returns the Average of the Measurements
  Complex PresentBinAverage();  // returns the Average of the Measurements in the present bin
  double PresentBinAverage2();  // returns the respective Average of the Measurements of squares 
  double Variance(); // returns the Variance of the Measurements
  double VarianceOfBins();
  double ErrorEstimate();
  void SetBinNumber(unsigned binnum); // set maximal bin number
  void SetMinMaxbin(int newMaxbins, int newMinObsPerBin); // Reset initial parameters
  
  void Rebin(unsigned int sequence);
  void Reset();
  
  const Complex& BinValue(int i) const  { return Values[i]; }
  const double& BinValue2(int i) const { return Values2[i];}
  const Complex& HistoryValue(int i) const  { return HistoryAverages[i]; }
  const double& HistoryValue2(int i) const { return HistoryBinVariances[i];}
  
  
  private:
  string name;
  int StartBins;        //arguments of constructor stored for Reset
  int StartObsPerBin; //dito
  int BinSize;       // number of Measurements per bin
  unsigned int MaxBinNum;    // maximum number of bins 
  int BinEntries;    // number of Measurements in last bin
  
  std::vector<Complex> Values;  // bin values
  std::vector<double> Values2; // bin values of squares
  std::vector<Complex> HistoryAverages;
  std::vector<double> HistoryBinVariances;
  
};
    

#endif // COMPLEXOBSERVABLE_H
