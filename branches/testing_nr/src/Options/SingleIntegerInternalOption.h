////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of single integer option                      //
//                                                                            //
//                        last modification : 19/08/2001                      //
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


#ifndef SINGLEINTEGERINTERNALOPTION_H
#define SINGLEINTEGERINTERNALOPTION_H


#include "config.h"
#include "Options/SingleIntegerOption.h"


class SingleIntegerInternalOption : public SingleIntegerOption
{

 protected:

  bool extFlag;

 public:


  // constructor from default datas
  //
  // optionCode = character associated to the option
  // optionName = string corresponding to option name
  // optionDescription = string describing option (used for -h option)
  // defaultValue = integer default value
  // external = indicates if visible to the exterior
  // minValueFlag = flag to indicates an integer minimum value
  // minValue = integer minimum value 
  // maxValueFlag = flag to indicates an integer maximum value
  // maxValue = integer maximum value (no maximum value if lower or equal to minValue) 
  SingleIntegerInternalOption(char optionCode, char* optionName, char* optionDescription, int defaultValue = 0,
			      bool external=false, 
		      bool minValueFlag = false, int minValue = 0, 
		      bool maxValueFlag = false, int maxValue = 0);

  // destructor
  //
  ~SingleIntegerInternalOption();


  // print the current option and its values
  //  
  // output = reference on output stream;
  // shortVersion = true if return only option code and the option value, false if return option description in addition
  // return value = reference on current output stream
  virtual ostream& DisplayOption (ostream& output, bool shortVersion = false);

  // print help concerning current option
  //
  // output = reference on output stream;
  // return value = reference on current output stream
  ostream& DisplayHelp (ostream& output);

  // Get read integer
  //
  // return value = integer value
  int GetInteger();

  // Set read integer
  //
  void SetInteger(int newValue);

  // Increase value
  void Increase();

};

#endif
