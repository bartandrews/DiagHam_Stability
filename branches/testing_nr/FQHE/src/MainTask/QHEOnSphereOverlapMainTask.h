////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of main task for qhe overlap on sphere                //
//                                                                            //
//                        last modification : 21/04/2005                      //
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


#ifndef QHEONSPHEREOVERLAPMAINTASK_H
#define QHEONSPHEREOVERLAPMAINTASK_H


#include "config.h"

#include "MainTask/AbstractMainTask.h"


class OptionManager;
class AbstractRandomNumberGenerator;
class Abstract1DComplexFunction;


class QHEOnSphereOverlapMainTask: public AbstractMainTask
{

 protected:

  // number of particles
  int NbrParticles;

  // pointer to the random number generator
  AbstractRandomNumberGenerator* RandomNumberGenerator;

  // pointer to the reference wave function
  Abstract1DComplexFunction* ReferenceWaveFunction;
  // pointer to the test wave function
  Abstract1DComplexFunction* TestWaveFunction;

  // number of iteration between two consecutive result displays
  int DisplayStep;
  // flag to indicate if intermediate values used for overlap calculation have to be shown
  bool ShowDetailFlag;
  // number of Monte Carlo iterations
  int NbrIterations;

 public:

  // constructor
  //  
  // options = pointer to the options managers containing all running options
  // referenceWaveFunction = pointer to the reference wave function
  // testWaveFunction = pointer to the test wave function
  // nbrParticles = number of particles
  // randomNumberGenerator = pointer to the random number generator
  QHEOnSphereOverlapMainTask(OptionManager* options, Abstract1DComplexFunction* referenceWaveFunction,
			     Abstract1DComplexFunction* testWaveFunction, int nbrParticles, 
			     AbstractRandomNumberGenerator* randomNumberGenerator);
  
  // destructor
  //  
  ~QHEOnSphereOverlapMainTask();
  
  // execute the main task
  // 
  // return value = 0 if no error occurs, else return error code
  int ExecuteMainTask();

};

#endif
