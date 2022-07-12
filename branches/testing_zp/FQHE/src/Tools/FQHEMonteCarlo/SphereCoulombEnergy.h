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


#ifndef SPHERECOULOMBENERGY_H
#define SPHERECOULOMBENERGY_H

#include "config.h"

#include "AbstractObservable.h"
#include "ParticleOnSphereCollection.h"
#include "MCObservables/WeightedRealObservable.h"

class SphereCoulombEnergy : public AbstractObservable
{
 protected:
  // Radius of the sphere in units of magnetic length
  double Radius;

  // Number of Flux quanta piercing the sphere
  int NbrFlux;
  
  // square width;
  double WidthSqr;

  // number of observations made
  int NbrObservations;

  // system that the observable operates on
  ParticleOnSphereCollection *System;

  // number of particles in System:
  int NbrParticles;

  // pointers to spinor coordinates (external)
  Complex *SpinorUCoordinates;
  Complex *SpinorVCoordinates;

  // core observable
  WeightedRealObservable *Values;
  
 public:

  // default constructor
  SphereCoulombEnergy();

  // constructor
  // nbrFlux = Number of Flux piercing sphere
  // width = simple model of finite layer width for 1/sqrt(w^2+r^2)
  SphereCoulombEnergy(int nbrFlux, double width=0.0);


  // destructor
  virtual ~SphereCoulombEnergy();

  // call to make an observation
  virtual void RecordValue(double weight);

  // print legend to the given stream
  // all = flag indicating whether to print all, or shortened information
  virtual void PrintLegend(std::ostream &output, bool all = false);

  // print status to the given stream
  // all = flag indicating whether to print all, or shortened information
  virtual void PrintStatus(std::ostream &output, bool all = false);

  // print formatted data suitable for plotting
  // ouput = the target stream
  virtual void WriteDataFile(std::ostream &output);

  // set particle collection that the observable operates on
  // system = particle collection
  virtual void SetParticleCollection(AbstractParticleCollection *system);
  
};

#endif
