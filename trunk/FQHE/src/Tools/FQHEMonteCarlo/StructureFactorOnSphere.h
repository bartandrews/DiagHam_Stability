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


#ifndef STRUCTUREFACTORONSPHERE_H
#define STRUCTUREFACTORONSPHERE_H

#include "config.h"
#include "AbstractObservable.h"
#include "ParticleOnSphereCollection.h"
#include "Polynomial/LegendrePolynomials.h"
#include "MCObservables/WeightedRealVectorObservable.h"
#include <iostream>


class StructureFactorOnSphere : public AbstractObservable
{
 protected:
    // Radius of the sphere in units of magnetic length
  double Radius;

  // Number of Flux quanta piercing the sphere
  int NbrFlux;

  // number of structure factors
  int NbrStructureFactor;

  // total weight of observations
  double Measures;

  // Number of particles
  int NbrParticles;

  // object for generating legendre polynomials
  LegendrePolynomials *LegendreBasis;

  // current set of legendre values
  double *CurrentLegendrePolynomials;

  // number of structure factors to calculate
  int NbrStructureFactors;

  // observable to perform statistics
  WeightedRealVectorObservable *StructureFactors;

   // system that the observable operates on
  ParticleOnSphereCollection *System;

  // pointers to spinor coordinates (external)
  Complex *SpinorUCoordinates;
  Complex *SpinorVCoordinates;

  // Flag indicating whether variable is printed by default
  bool PrintFlag;
  
 public:

  // standard constructor
  StructureFactorOnSphere();
  
  // constructor
  // nbrFlux = number of flux piercing the sphere
  StructureFactorOnSphere(int nbrFlux, int nbrStructureFactors);
  
  // destructor
  virtual ~StructureFactorOnSphere();

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
  
};

#endif // STRUCTUREFACTORONSPHERE_H
