////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of FQHE wave function manager                    //
//                                                                            //
//                        last modification : 18/01/2005                      //
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


#ifndef QHESAMPLINGFUNCTIONMANAGER_H
#define QHESAMPLINGFUNCTIONMANAGER_H


#include "config.h"

#include <iostream>


using std::ostream;


class OptionManager;
class Abstract1DComplexFunction;
class AbstractMCSamplingFunction;


class QHESamplingFunctionManager
{

 protected:

  // pointer to the option manager
  OptionManager* Options;

  // id of the geometry to use
  int GeometryID;

  int SamplingfunctionID;

 public:

  // list of avalaible geometries
  enum Geometries
    {
      SphereGeometry = 0x01,
      DiskGeometry = 0x02,
      SphereWithSpinGeometry = 0x04
    };

  // list of available wavefunctions:
  enum WaveFunctions
    {
      InvalidWaveFunction = 0x00000,
      Laughlin = 0x00103,
    };
  
  // constructor
  //
  // geometry = id of the geometry to use
  QHESamplingFunctionManager(int geometry = QHESamplingFunctionManager::SphereGeometry);

  // destructor
  //
  ~QHESamplingFunctionManager();

  // add an option group containing all options related to the wave functions
  //
  // manager = pointer to the option manager
  void AddOptionGroup(OptionManager* manager);

  // get list of all available sampling functions
  // 
  // str = reference on the output stream
  ostream& ShowAvalaibleSamplingFunctions (ostream& str);

  // test if list of sampling functions should be displayed
  void TestForShow(ostream& str=std::cout);
  
  // get the wave function corresponding to the option constraints
  //
  // return value = pointer to the wave function (null if an error occurs)
  AbstractMCSamplingFunction* GetSamplingFunction();

  // get a description of the returned wave function
  //
  char* GetDescription();

  // get type of WaveFunction

  int GetWaveFunctionType();
  
};

#endif
