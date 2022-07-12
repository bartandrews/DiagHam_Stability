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


#ifndef QHEWAVEFUNCTIONMANAGER_H
#define QHEWAVEFUNCTIONMANAGER_H


#include "config.h"

#include <iostream>


using std::ostream;


class OptionManager;
class Abstract1DComplexFunction;


class QHEWaveFunctionManager
{

 protected:

  // pointer to the option manager
  OptionManager* Options;

  // id of the geometry to use
  int GeometryID;

  int WavefunctionID;

 public:

  // list of avalaible geometries
  enum Geometries
    {
      SphereGeometry = 0x01,
      DiskGeometry = 0x02,
      SphereWithSpinGeometry = 0x04,
      SphereWithSU3SpinGeometry = 0x08
    };

  // list of available wavefunctions:
  enum WaveFunctions
    {
      InvalidWaveFunction = 0x00000,
      Laughlin = 0x00103,
      Pfaffian = 0x00203,
      Pfaffian2QH = 0x00401,
      ReadRezayi = 0x00803,
      FilledCF = 0x01101,
      GenericCF =0x01001,
      UnprojectedCF = 0x02001,
      PairedCF =  0x14005,
      PairedCFCB = 0x18004,
      OneOneOne = 0x04004,
      OneS = 0x100000,
      TwoThirdsS = 0x200000,
      TwoThirdsUnpolarized = 0x200001,
      HaldaneRezayi = 0x400000,
      ExtendedHalperin = 0x800000,
      HundRuleSinglet = 0xf00000,
      TrialWaveFunction = 0x10000,
      Halperin = 0x1000000,
      SLBS = 0x2000000,
      SLBSV = 0x2010000
    };
  
  // constructor
  //
  // geometry = id of the geometry to use
  QHEWaveFunctionManager(int geometry = QHEWaveFunctionManager::SphereGeometry);

  // destructor
  //
  ~QHEWaveFunctionManager();

  // add an option group containing all options related to the wave functions
  //
  // manager = pointer to the option manager
  void AddOptionGroup(OptionManager* manager);

  // get list of all available wave functions
  // 
  // str = reference on the output stream
  ostream& ShowAvalaibleWaveFunctions (ostream& str);
  
  // get the wave function corresponding to the option constraints
  //
  // return value = pointer to the wave function (null if an error occurs)
  Abstract1DComplexFunction* GetWaveFunction();

  // get a description of the returned wave function
  //
  char* GetDescription();

  // get type of WaveFunction

  int GetWaveFunctionType();
  
};

#endif
