////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of abstract DMRG algorithm                     //
//                                                                            //
//                        last modification : 16/04/2001                      //
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


#ifndef ABSTRACTDMRGALGORITHM_H
#define ABSTRACTDMRGALGORITHM_H


class AbstractQuantumNumber;


class AbstractDMRGAlgorithm
{

 public:

  // virtual destructor
  //
  virtual ~AbstractDMRGAlgorithm ();

  // force constraint on global quantum number
  //
  // quantumNumber = pointer to the global quantum number to use
  virtual void Constraint(AbstractQuantumNumber* quantumNumber) = 0;

  // run DMRG algorithm
  //
  virtual void RunDMRG(int currentBlockIndex = 0) = 0;
 
};

#endif
