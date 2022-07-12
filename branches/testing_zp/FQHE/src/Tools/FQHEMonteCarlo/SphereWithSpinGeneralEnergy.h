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


#ifndef SPHEREGENERALENERGY_H
#define SPHEREGENERALENERGY_H

#include "config.h"

#include "AbstractObservable.h"
#include "ParticleOnSphereCollection.h"
#include "MCObservables/WeightedRealObservable.h"

class SphereWithSpinGeneralEnergy : public AbstractObservable
{
 public:
  enum InteractionTypes
    {
      Unknown = 0x0000,
      Polynomial = 0x0001,
      AsymptoticExp = 0x0002
    };
  
 protected:
  // Type of interaction
  int InteractionType;
  
  // Radius of the sphere in units of magnetic length
  double Radius;

  // Number of Flux quanta piercing the sphere
  int NbrFlux;
  
  // Number of parameters
  int NbrParametersInter;
  int NbrParametersIntra;
  int NbrParametersIntra2;

  // Flag indicating if second interaction is being used
  bool HaveIntra2;

  // Coefficients of polynomial interaction
  double *CoefficientsInter;
  double *CoefficientsIntra;
  double *CoefficientsIntra2;

  // Parameters for asymptotic form of interaction
  // (number of values)
  int NbrAsymptoticsInter;
  int NbrAsymptoticsIntra;
  int NbrAsymptoticsIntra2;
  // (values for prefactors of asymptotic terms)
  double *AsymptoticsInter;
  double *AsymptoticsIntra;
  double *AsymptoticsIntra2;
  // (values for regularization of asymptotic terms)
  double *AsymptoticsInterReg;
  double *AsymptoticsIntraReg;
  double *AsymptoticsIntra2Reg;

  // number of observations made
  int NbrObservations;

  // system that the observable operates on
  ParticleOnSphereCollection *System;

  // number of particles in System:
  int NbrParticles;
  // number of particles with up spin:
  int NbrUp;

  // pointers to spinor coordinates (external)
  Complex *SpinorUCoordinates;
  Complex *SpinorVCoordinates;

  // core observable
  WeightedRealObservable *Values;

  // tables for storage of distances and their powers
  double **RijSq;
  double ***RijSqPowers;
  int NumSqPowers;
  double **GaussianIJ;
  
  
 public:

  // default constructor
  SphereWithSpinGeneralEnergy();

  // constructor (old style, multiple parameter files)
  // nbrUp = Number of particles with up spin
  // nbrFlux = Number of Flux piercing sphere
  // parametersInter = file describing parameters of the interaction (inter-spin)
  // parametersIntra = file describing parameters of the interaction (intra-spin)
  // parametersIntra2 = file describing parameters of the interaction (intra-spin, other spin-channel)
  SphereWithSpinGeneralEnergy(int nbrUp, int nbrFlux, const char* parametersInter, const char* parametersIntra, const char* parametersIntra2=NULL);

  // constructor (new style, single parameter file)
  // nbrUp = Number of particles with up spin
  // nbrFlux = Number of Flux piercing sphere
  // parameters = file describing parameters of the interaction (all spin channels UpUp, UpDown, DownDown)
  SphereWithSpinGeneralEnergy(int nbrUp, int nbrFlux, const char* parameters);

  // destructor
  virtual ~SphereWithSpinGeneralEnergy();

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

  // additional routines for energy observables:
  
  // returns the total background energy
  double GetTotalBackgroundEnergy();

 private:

  // evaluate exponentials and powers of r^2
  void EvaluateGaussianTables();
  
};

#endif
