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


#ifndef SPHEREBLOCKCOULOMBENERGY_H
#define SPHEREBLOCKCOULOMBENERGY_H

#include "config.h"

#include "AbstractBlockObservable.h"
#include "ParticleOnSphereCollection.h"
#include "MCObservables/ComplexVectorObservable.h"
#include "MCObservables/RealObservable.h"


class SphereBlockCoulombEnergy : public AbstractBlockObservable
{
 protected:
  // Radius of the sphere in units of magnetic length
  double Radius;

  // Number of Flux quanta piercing the sphere
  int NbrFlux;

  // Number of Blocks
  int NbrBlocks;
  
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

  // flag indicating whether to disregard imaginary part of weighhts
  bool UseComplex;

  // core observable
  ComplexVectorObservable *ComplexValues;
  // alternative format avoiding complex numbers
  RealObservable *RealValues;

  // temporary storage
  Complex *TempComplex;
  
 public:

  // default constructor
  SphereBlockCoulombEnergy();

  // constructor
  // nbrBlockEntries = total number of entries for block vector observable  
  // nbrFlux = Number of Flux piercing sphere
  // width = simple model of finite layer width for 1/sqrt(w^2+r^2)
  SphereBlockCoulombEnergy(int nbrBlockEntries, int nbrFlux, double width=0.0, bool useComplex=true);


  // destructor
  virtual ~SphereBlockCoulombEnergy();

  // call to make an observation
  virtual void RecordValue(double weight);

    // call to make an observation
  // weight = relative weight of this sample
  virtual void RecordRealValue(Complex *weights);

  // call to make an observation
  // weight = relative weight of this sample
  virtual void RecordValue(Complex *weights);

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

  // additional routines for energy observables:
  // returns the total background energy
  double GetTotalBackgroundEnergy();

 private:

  // calculate energy
  inline double GetEnergy();

  
};
// calculate energy
inline double SphereBlockCoulombEnergy::GetEnergy()
{
  int N = this->NbrParticles;
  ++NbrObservations;
  double E=0.0;
  if (WidthSqr>0.0)
    {
      for (int i=1;i<N;i++)
	for(int j=0;j<i;j++)
	  E+=1.0/sqrt(WidthSqr+SqrNorm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]));
    }
  else
    {
      for (int i=1;i<N;i++)
	for(int j=0;j<i;j++)
	  E+=1.0/Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
    }
  return 0.5*E/Radius;
}

#endif
