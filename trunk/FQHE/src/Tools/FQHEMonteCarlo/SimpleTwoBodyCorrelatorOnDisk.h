////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//        class for a binned correlation function on the sphere geometry      //
//                                                                            //
//                        last modification : 19/10/2009                      //
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


#ifndef SIMPLETWOBODYCORRELATORONDISK_H
#define SIMPLETWOBODYCORRELATORONDISK_H

#include "config.h"
#include "AbstractObservable.h"
#include "ParticleOnDiskCollection.h"
#include <iostream>
#include <cmath>

class SimpleTwoBodyCorrelatorOnDisk : public AbstractObservable
{
 protected:
    // maximum radius up to which samples need to be recorded
  double MaxRadius;

  // number of bins
  int Bins;

  // overall resolution (number of bins on regular axis)
  int Resolution;

  // increased resolution for the interval close to the origin (in place of Range bins)
  int Highres;

  // number of datapoints to be included in that interval
  int Range;

  // ratio of highres area
  double HighResRatio;

  // total weight of observations
  double Measures;

  // array of measurements
  double *Correlations;

  // Number of particles
  int NbrParticles;

   // system that the observable operates on
  ParticleOnDiskCollection *System;

  // pointers to coordinates (external)
  Complex *CoordinatesZ;

  // Flag indicating whether variable is printed by default
  bool PrintFlag;
  
 public:

  // standard constructor
  SimpleTwoBodyCorrelatorOnDisk();
  
  // constructor
  // rMax = maximum radius to represent
  // resolution = total number of bins
  // highres = number of points in high resolution interval at small r
  // range =  ranger over which high resolution is implemented
  SimpleTwoBodyCorrelatorOnDisk(double rMax, int resolution, int highres, int range=0);
  
  // destructor
  virtual ~SimpleTwoBodyCorrelatorOnDisk();

  // call to make an observation
  // weight = relative weight of this sample
  virtual void RecordValue(double weight);

  // print legend to the given stream
  // all = flag indicating whether to print all, or shortened information
  virtual void PrintLegend(std::ostream &output, bool all = false);

  // print status to the given stream
  // all = flag indicating whether to print all, or shortened information
  virtual void PrintStatus(std::ostream &output, bool all = false);

  // request whether observable should be printed
  //
  virtual bool IncludeInPrint();

  // set print status
  //
  virtual void IncludeInPrint(bool newStatus);

  // print formatted data suitable for plotting
  // ouput = the target stream
  virtual void WriteDataFile(std::ostream &output);

  // write binary data 
  // ouput = the target stream
  virtual void WriteBinaryData(std::ostream &output);

  // set particle collection that the observable operates on
  // system = particle collection
  virtual void SetParticleCollection(AbstractParticleCollection *system);

  // accessor function to return the legend and numerical values for legend
  virtual void GetVectorLegend(std::string &legendParameters, std::string &legendValue, RealVector &parameterValues);

  // accessor function for average and error for variables with real measurements
  virtual void GetRealVectorMeasurement(RealVector &values, RealVector &errors);
  
 private:
  // get bin index for a given radius
  int GetIndex(double radius)
  {
    if (radius > this->MaxRadius) return this->Resolution;
    else
      return (int)(std::pow(radius/this->MaxRadius,2)*this->Resolution);
  }

  // get bin index for a radius on the high resolution grid
  int GetHighResIndex(double radius)
  {
    if (radius > this->MaxRadius) return this->Bins;
    else
      {
	// std::cout << "radius="<<radius<<", maxR="<<MaxRadius<<", HighResRatio="<<HighResRatio<<", Res="<<Resolution <<std::endl;
	return (int)(std::pow(radius/this->MaxRadius,2)*this->HighResRatio*this->Resolution);
      }
  }
    
  // get bin radius for a (continuous) bin index on regular scale
  double GetBinRadius(double binNbr)
  {
    return std::sqrt( binNbr/this->Resolution) * this->MaxRadius;
  }

  // get bin radius for a (continuous) bin index on high resolution scale
  double GetHighResBinRadius(double binNbr)
  {
    return std::sqrt( binNbr/this->Resolution) * this->MaxRadius * this->HighResRatio;
  }

  
};

#endif
