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


#ifndef ABSTRACTOBSERVABLE_H
#define ABSTRACTOBSERVABLE_H

#include "config.h"
#include "GeneralTools/Warnings.h"
#include "Vector/RealVector.h"
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include <iostream>

class AbstractParticleCollection;

class AbstractObservable
{
 protected:
  unsigned Type;
  
 public:
  
  enum Properties{
    RealObservableT = 0x1u,
    ComplexObservableT = 0x2u,
    VectorValued = 0x4u
  };
  
  // destructor
  virtual ~AbstractObservable();

  // call to make an observation
  // weight = relative weight of this sample
  virtual void RecordValue(double weight) = 0;

  // print legend to the given stream
  // all = flag indicating whether to print all, or shortened information
  virtual void PrintLegend(std::ostream &output, bool all = false) = 0;

  // print status to the given stream
  // all = flag indicating whether to print all, or shortened information
  virtual void PrintStatus(std::ostream &output, bool all = false) = 0;

  // request whether observable should be printed
  //
  virtual bool IncludeInPrint();

  // set print status
  //
  virtual void IncludeInPrint(bool newStatus);

  // print formatted data suitable for plotting
  // ouput = the target stream
  virtual void WriteDataFile(std::ostream &output) = 0;

  // set particle collection that the observable operates on
  // system = particle collection
  virtual void SetParticleCollection(AbstractParticleCollection *system) = 0;

  // check for real data
  bool IsReal(){return (this->Type & RealObservableT)!=0;}

  // check for complex data
  bool IsComplex(){return (this->Type & ComplexObservableT)!=0;}

  // check for real data
  bool IsVectorValued(){return (this->Type & VectorValued)!=0;}

  // accessor function to query a unique identifier for the observable
  virtual std::string GetIdentifier()  { NoOverload(); return "";}

  // accessor function to return the legend corresponding to the value of the observable
  virtual std::string GetLegend()  { NoOverload(); return "";}

  // accessor function for average and error for variables with real measurements
  virtual void GetRealMeasurement(double &value, double &error)  { NoOverload(); }

  // accessor function for average and error for variables with Complex measurements
  virtual void GetComplexMeasurement(Complex &value, double &error)  { NoOverload(); }

  // accessor function to return the legend and numerical values for legend
  virtual void GetVectorLegend(std::string &legendParameters, std::string &legendValue, RealVector &parameterValues)  { NoOverload(); }

  // accessor function for average and error for variables with real measurements
  // errors is returned as a vector of length zero if the observable does not provide error estimates.
  virtual void GetRealVectorMeasurement(RealVector &values, RealVector &errors)  { NoOverload(); }

  // accessor function for average and error for variables with Complex measurements
  // errors is returned as a vector of length zero if the observable does not provide error estimates.
  virtual void GetComplexVectorMeasurement(ComplexVector &values, RealVector &errors)  { NoOverload(); }
  
};

#endif
