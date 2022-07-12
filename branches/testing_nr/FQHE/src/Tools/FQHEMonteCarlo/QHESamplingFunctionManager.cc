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


#include "config.h"
#include "Options/Options.h"

#include "QHESamplingFunctionManager.h"
#include "AbstractMCSamplingFunction.h"
#include "LaughlinWithSpinSamplingFunction.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include <iostream>


using std::endl;
using std::cout;

// constructor
//
// geometry = id of the geometry to use

QHESamplingFunctionManager::QHESamplingFunctionManager(int geometry)
{
  this->GeometryID = geometry;
  this->Options = 0;
}

// destructor
//

QHESamplingFunctionManager::~QHESamplingFunctionManager()
{
}

// add an option group containing all options related to the wave functions
//
// manager = pointer to the option manager

void QHESamplingFunctionManager::AddOptionGroup(OptionManager* manager)
{
  this->Options = manager;
  OptionGroup* SamplingFunctionGroup  = new OptionGroup ("sampling function options");
  (*(this->Options)) += SamplingFunctionGroup;
  (*SamplingFunctionGroup) += new BooleanOption ('\n', "list-sampling-functions", "list all available sampling fuctions");
  (*SamplingFunctionGroup) += new SingleStringOption  ('\n', "sampler", "name of the test wave fuction",0);
  if (this->GeometryID & QHESamplingFunctionManager::SphereGeometry)
    {
      //(*SamplingFunctionGroup) += new SingleStringOption  ('\n', "xxx", "description",0);
    }
  else if (this->GeometryID & QHESamplingFunctionManager::SphereWithSpinGeometry)
    {
      // PairedCF(CB)Options:
      (*SamplingFunctionGroup) += new SingleIntegerOption  ('\n', "laughlin-exponent", "power to which the jastrow factors in sampling function are raised",2);
    }
}

// get list of all available sampling functions
// 
// str = reference on the output stream

ostream& QHESamplingFunctionManager::ShowAvalaibleSamplingFunctions (ostream& str)
{
  str << "list of avalaible sampling functions:" << endl;
  if (this->GeometryID == QHESamplingFunctionManager::SphereGeometry)
    {
      str << "  no sampling functions, yet" << endl;	
    }
  else
    if (this->GeometryID == QHESamplingFunctionManager::DiskGeometry)
      {
	str << "  no sampling functions, yet" << endl;	
      }
    else
      if (this->GeometryID == QHESamplingFunctionManager::SphereWithSpinGeometry)
	{	  
	  str << "  * laughlin : uncorrelated layer each with a laughlin-type wavefunction" << endl;
	}
  return str;
}


void QHESamplingFunctionManager::TestForShow(ostream& str)
{
  if (Options->GetBoolean("list-samplingfunctions") == true)
    {
      this->ShowAvalaibleSamplingFunctions(str);
      exit(0);
    }
}
  
// get the wave function corresponding to the option constraints
//
// return value = pointer to the wave function (null if an error occurs)
AbstractMCSamplingFunction* QHESamplingFunctionManager::GetSamplingFunction()
{
  if ((*(this->Options))["sampler"] == 0)
    {
      return 0;
    }
  if (this->Options->GetString("sampler") == 0)
    {
      return 0;
    }
  if (this->GeometryID == QHESamplingFunctionManager::SphereGeometry)
    {
      // none implemented for the moment
      return 0;
    }
  else
    if (this->GeometryID == QHESamplingFunctionManager::DiskGeometry)
      {
	// none implemented for the moment
	return 0;
      }
    else
      if (this->GeometryID == QHESamplingFunctionManager::SphereWithSpinGeometry)
	{
	  if ((strcmp (this->Options->GetString("sampler"), "laughlin") == 0))
	    {
	      int N = this->Options->GetInteger("nbr-particles");	      
	      int m = this->Options->GetInteger("laughlin-exponent");
	      AbstractMCSamplingFunction* rst  = new LaughlinWithSpinSamplingFunction(N,m);
	      return rst;
	    }
	  return 0;
	}
  return 0;
}

char* QHESamplingFunctionManager::GetDescription()
{
  if ((*(this->Options))["sampler"] == 0)
    {
      return 0;
    }
  char * buffer = new char[1000];
  sprintf(buffer,"%s N=%d",this->Options->GetString("sampler"), this->Options->GetInteger("nbr-particles"));
  char *rst = new char[strlen(buffer)+1];
  strcpy(rst,buffer);
  delete [] buffer;
  return rst;
}


int QHESamplingFunctionManager::GetWaveFunctionType()
{
  if ((*(this->Options))["sampler"] == 0)
    return QHESamplingFunctionManager::InvalidWaveFunction;
  if (strcmp (this->Options->GetString("sampler"), "laughlin") == 0)
    return QHESamplingFunctionManager::Laughlin;  
  return QHESamplingFunctionManager::InvalidWaveFunction;
}
