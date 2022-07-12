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
#include "LaughlinSamplingFunction.h"
#include "HalperinSamplingFunction.h"
#include "MRBlockSamplingFunction.h"
#include "LaughlinSamplingFunctionOnDisk.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include <iostream>
#include <cstring>
#include <stdio.h>

using std::endl;
using std::cout;

// constructor
//
// geometry = id of the geometry to use

QHESamplingFunctionManager::QHESamplingFunctionManager(int geometry, bool blockAlgorithm)
{
  this->GeometryID = geometry;
  this->Options = 0;
  this->BlockAlgorithm = blockAlgorithm;
}

// destructor
//

QHESamplingFunctionManager::~QHESamplingFunctionManager()
{
}

// add an option group containing all options related to the wave functions
//
// manager = pointer to the option manager
// blockAlgorithm = flag indicating whether block symmetric Monte Carlo shall be used
void QHESamplingFunctionManager::AddOptionGroup(OptionManager* manager)
{
  this->Options = manager;
  OptionGroup* SamplingFunctionGroup  = new OptionGroup ("sampling function options");
  (*(this->Options)) += SamplingFunctionGroup;
  (*SamplingFunctionGroup) += new BooleanOption ('\n', "list-sampling-functions", "list all available sampling fuctions");
  (*SamplingFunctionGroup) += new SingleStringOption  ('\n', "sampler", "name of the test wave fuction",0);
  if (this->BlockAlgorithm)
    {
      if (this->GeometryID & QHESamplingFunctionManager::SphereGeometry)
	{
	  (*SamplingFunctionGroup) += new SingleIntegerOption  ('\n', "jastrow-exponent", "power to which the jastrow factor in sampling function is raised",1);
	  //(*SamplingFunctionGroup) += new SingleStringOption  ('\n', "xxx", "description",0);
	}
    }
  else
    {
      if (this->GeometryID & QHESamplingFunctionManager::SphereGeometry)
	{
	  (*SamplingFunctionGroup) += new SingleIntegerOption  ('\n', "laughlin-exponent", "power to which the jastrow factors in sampling function are raised",2);
	  //(*SamplingFunctionGroup) += new SingleStringOption  ('\n', "xxx", "description",0);
	}
      else if (this->GeometryID & QHESamplingFunctionManager::SphereWithSpinGeometry)
	{
	  (*SamplingFunctionGroup) += new SingleIntegerOption  ('\n', "laughlin-exponent", "power to which the jastrow factors in sampling function are raised",2);
	  (*SamplingFunctionGroup) += new MultipleIntegerOption  ('\n', "SHC", "coefficients (k,l,m) of sampling Halperin wavefunction",',' ,',', "1,1,1");
	}
      else if (this->GeometryID & QHESamplingFunctionManager::DiskGeometry)
	{
	  (*SamplingFunctionGroup) += new SingleIntegerOption  ('\n', "laughlin-exponent", "power to which the jastrow factors in sampling function are raised",2);
	  (*SamplingFunctionGroup) += new SingleDoubleOption  ('\n', "defect-angle", "defect angle for disk geometry (units of 2pi)", 0.0, true, 0.0, true, 1.0);
	  (*SamplingFunctionGroup) += new SingleDoubleOption  ('\n', "geometric-spin", "spin j for Laughlin state", 0.0);
	}

    }
}

// get list of all available sampling functions
// 
// str = reference on the output stream

ostream& QHESamplingFunctionManager::ShowAvalaibleSamplingFunctions (ostream& str)
{
  str << "list of avalaible sampling functions (indicate one with option --sampler):" << endl;
  if (this->BlockAlgorithm)
    {
      if (this->GeometryID == QHESamplingFunctionManager::SphereGeometry)
	{
	  str << "  * MR: Moore-Read wavefunction" << endl;
	}
    }
  else
    {
      if (this->GeometryID == QHESamplingFunctionManager::SphereGeometry)
	{
	  str << "  * laughlin : laughlin-type wavefunction" << endl;
	}
      else
	if (this->GeometryID == QHESamplingFunctionManager::DiskGeometry)
	  {
	    str << "  * laughlin : laughlin-type wavefunction" << endl;
	  }
	else
	  if (this->GeometryID == QHESamplingFunctionManager::SphereWithSpinGeometry)
	    {	  
	      str << "  * laughlin : uncorrelated layer each with a laughlin-type wavefunction" << endl;
	      str << "  * halperin : general halperin wavefunctions (z-z)^k (w-w)^l (z-w)^m (set by SHC)" << endl;
	    }
    }
  return str;
}


void QHESamplingFunctionManager::TestForShow(ostream& str)
{
  if (Options->GetBoolean("list-sampling-functions") == true)
    {
      this->ShowAvalaibleSamplingFunctions(str);
      exit(0);
    }
}


// get the wave function corresponding to the option constraints
//
// return value = pointer to the wave function (null if an error occurs)
AbstractMCBlockSamplingFunction* QHESamplingFunctionManager::GetBlockSamplingFunction()
{
  if (this->BlockAlgorithm)
    return (AbstractMCBlockSamplingFunction*)(this->GetSamplingFunction());
  else
    return 0;
}

// get the wave function corresponding to the option constraints
//
// return value = pointer to the wave function (null if an error occurs)
AbstractMCBlockSamplingFunctionOnSphere* QHESamplingFunctionManager::GetBlockSamplingFunctionOnSphere()
{
  if ((this->BlockAlgorithm)&&(this->GeometryID & QHESamplingFunctionManager::SphereGeometry))
    return (AbstractMCBlockSamplingFunctionOnSphere*)(this->GetSamplingFunction());
  else
    return 0;
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
  if (this->BlockAlgorithm)
    {
      if (this->GeometryID == QHESamplingFunctionManager::SphereGeometry)
	{
	  if ((strcmp (this->Options->GetString("sampler"), "MR") == 0))
	    {
	      int N = this->Options->GetInteger("nbr-particles");
	      int J = this->Options->GetInteger("jastrow-exponent");
	      if (N&1)
		{
		  cout << "Require even number of particles in MR wavefunction"<<endl;
		  exit(1);
		}
	      AbstractMCSamplingFunction* rst  = new MRBlockSamplingFunction(N/2, J /*, sqrCriticalDistance*/);
	      return rst;
	    }
	  return 0;
	}
    }
  else // no block algorithm
    { 
      if (this->GeometryID == QHESamplingFunctionManager::SphereGeometry)
	{
	  if ((strcmp (this->Options->GetString("sampler"), "laughlin") == 0))
	    {
	      int N = this->Options->GetInteger("nbr-particles");	      
	      int m = this->Options->GetInteger("laughlin-exponent");
	      AbstractMCSamplingFunction* rst  = new LaughlinSamplingFunction(N,m);
	      return rst;
	    }
	  return 0;
	}
      else
	if (this->GeometryID == QHESamplingFunctionManager::DiskGeometry)
	  {
	  if ((strcmp (this->Options->GetString("sampler"), "laughlin") == 0))
	    {
	      int N = this->Options->GetInteger("nbr-particles");	      
	      int m = this->Options->GetInteger("laughlin-exponent");
	      double alpha = this->Options->GetDouble("defect-angle");
	      double spin = this->Options->GetDouble("geometric-spin");
	      AbstractMCSamplingFunction* rst  = new LaughlinSamplingFunctionOnDisk(N,m,alpha,spin);
	      return rst;
	    }
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
	      if ((strcmp (this->Options->GetString("sampler"), "halperin") == 0))
		{
		  int N= this->Options->GetInteger("nbr-particles");
		  int length;
		  int *Params = Options->GetIntegers("SHC",length);
		  if (length != 3)
		    {
		      cout << "the Halperin sampling function requires precisely three parameters -SHC k,l,m !"<<endl;
		      exit(-1);
		    }
		  int K,L,M;
		  K = Params[0];
		  L = Params[1];
		  M = Params[2];
		  delete [] Params;
		  int NbrUp=-1;
		  if ((*this->Options)["total-sz"]!=NULL)
		    NbrUp = (N + this->Options->GetInteger("total-sz"))/2;
		  HalperinSamplingFunction* rst = new HalperinSamplingFunction(N, K, L, M, NbrUp);
		  return rst;
		}
	      return 0;
	    }
    }
  return 0;
}

char* QHESamplingFunctionManager::GetDescription()
{
  if ((*(this->Options))["sampler"] == 0)
    {
      return 0;
    }
  if (this->Options->GetString("sampler") == 0)
    return 0;
  char * buffer = new char[1000];
  sprintf(buffer,"%s N=%ld",this->Options->GetString("sampler"), this->Options->GetInteger("nbr-particles"));
  char *rst = new char[strlen(buffer)+1];
  strcpy(rst,buffer);
  if ((this->GeometryID == QHESamplingFunctionManager::SphereWithSpinGeometry)
      &&(strcmp (this->Options->GetString("sampler"), "laughlin") == 0))
    {
      int length;
      int *Params = Options->GetIntegers("SHC",length);
      if (length != 3)
	{
	  cout << "the Halperin sampling function requires precisely three parameters -SHC k,l,m !"<<endl;
	  exit(-1);
	}
      int K,L,M;
      K = Params[0];
      L = Params[1];
      M = Params[2];
      delete [] Params;      
      sprintf(buffer,"%s [k,l,m]=[%d,%d,%d]",rst,K,L,M);
      delete [] rst;
      rst = new char[strlen(buffer)+1];
      strcpy(rst,buffer);
    }
  delete [] buffer;
  return rst;
}


int QHESamplingFunctionManager::GetWaveFunctionType()
{
  if ((*(this->Options))["sampler"] == 0)
    return QHESamplingFunctionManager::InvalidWaveFunction;
  if (BlockAlgorithm)
    {
      if (strcmp (this->Options->GetString("sampler"), "MR") == 0)
	return QHESamplingFunctionManager::MooreReadBlock;
    }
  else
    {
      if (strcmp (this->Options->GetString("sampler"), "laughlin") == 0)
	return QHESamplingFunctionManager::Laughlin;
      if (strcmp (this->Options->GetString("sampler"), "halperin") == 0)
	return QHESamplingFunctionManager::Halperin;
    }
  return QHESamplingFunctionManager::InvalidWaveFunction;
}
