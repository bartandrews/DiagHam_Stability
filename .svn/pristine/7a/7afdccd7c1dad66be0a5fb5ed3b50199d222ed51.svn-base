////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of Architecture Manager                      //
//                                                                            //
//                        last modification : 28/05/2004                      //
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

#include "Architecture/AbstractArchitecture.h"

#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithGroundState.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedBlockLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithGroundStateDiskStorage.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"


// default constructor
//
// complexFlag = use complex Lanczos algorithms if true

LanczosManager::LanczosManager(bool complexFlag)
{
  this->LanczosAlgorithm = 0;
  this->Options = 0;
  this->ComplexFlag = complexFlag;
}

// destructor
//

LanczosManager::~LanczosManager()
{
  if (this->LanczosAlgorithm != 0)
    {
      delete this->LanczosAlgorithm;
    }
}

// add an option group containing all options related to the architecture 
//
// manager = pointer to the option manager

void LanczosManager::AddOptionGroup(OptionManager* manager)
{
  this->Options = manager;
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  (*(this->Options)) += LanczosGroup;
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "full-diag", 
						"maximum Hilbert space dimension for which full diagonalization is applied", 
						500, true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup) += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "limit-time", "use limit in time instead of a number of lanczos iteration (0 if none, time in seconds)", 0);
  (*LanczosGroup) += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new SingleStringOption  ('\n', "set-reorthogonalize", 
					      "force reorthogonalization with a set of vectors describe in the name of file passed through this option");
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-blockvectors", "use file that describe a set of initial vectors for the block Lanczos algorithm (syntax : InitialVectors=vec0.vec vec1.vec ...)", 0);
  (*LanczosGroup) += new BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);
  (*LanczosGroup) += new SingleDoubleOption ('\n', "lanczos-precision", "define Lanczos precision for eigenvalues (0 if automatically defined by the program)", 0);
  (*LanczosGroup) += new BooleanOption ('\n', "fast-disk", "use disk storage to increase speed of ground state calculation and decrease memory footprint when using Lanczos algorithm");
  (*LanczosGroup) += new BooleanOption ('\n', "resume-fastdisk", "resume the fast-disk mode Lanczos algorithm from a stopped one (for example due to computer crash)");
  
}

// get the best avalaible Lanczos algorithm in agrement with the option constraints
//
// architecture = pointer to the architecture to use within the Lanczos algorithm
// forceEigenstateComputation = if true, force computation of eigenstates
// return value = pointer to the Lanczos algorithm


AbstractLanczosAlgorithm* LanczosManager::GetLanczosAlgorithm(AbstractArchitecture* architecture, bool forceEigenstateComputation)
{
  if ((this->Options != 0) && (this->LanczosAlgorithm == 0))
    {      
      bool ResumeFlag = ((BooleanOption*) (*(this->Options))["resume"])->GetBoolean();
      bool DiskFlag = ((BooleanOption*) (*(this->Options))["disk"])->GetBoolean();
      int MaxNbrIterLanczos = ((SingleIntegerOption*) (*(this->Options))["iter-max"])->GetInteger();
      int NbrIterLanczos = ((SingleIntegerOption*) (*(this->Options))["nbr-iter"])->GetInteger();
      int NbrEigenvalue = ((SingleIntegerOption*) (*(this->Options))["nbr-eigen"])->GetInteger();
      bool BlockLanczosFlag = ((BooleanOption*) (*(this->Options))["block-lanczos"])->GetBoolean();
      bool SizeBlockLanczos = ((SingleIntegerOption*) (*(this->Options))["block-size"])->GetInteger();
      bool VectorMemory = ((SingleIntegerOption*) (*(this->Options))["nbr-vector"])->GetInteger();
      bool EvaluateEigenvectors = ((BooleanOption*) (*(this->Options))["eigenstate"])->GetBoolean();
      char* InitialVectorFileName = ((SingleStringOption*) (*(this->Options))["initial-vector"])->GetString();
      bool PartialLanczos = ((BooleanOption*) (*(this->Options))["partial-lanczos"])->GetBoolean();
      double LanczosPrecision = ((SingleDoubleOption*) (*(this->Options))["lanczos-precision"])->GetDouble();
      bool FastDiskFlag = ((BooleanOption*) (*(this->Options))["fast-disk"])->GetBoolean();
      bool ResumeFastDiskFlag = ((BooleanOption*) (*(this->Options))["resume-fastdisk"])->GetBoolean();
      bool FullReorthogonalizationFlag = ((BooleanOption*) (*(this->Options))["force-reorthogonalize"])->GetBoolean();
      if ((NbrEigenvalue > 1) || (EvaluateEigenvectors == false))
	{
	  FastDiskFlag = false;
	  ResumeFastDiskFlag = false;
	}
      if (this->ComplexFlag == false)
	{
	  if ((NbrEigenvalue == 1) && (FullReorthogonalizationFlag == false))
	    {
	      if (DiskFlag == false)
		if ((EvaluateEigenvectors == true) || (forceEigenstateComputation == true))
		  this->LanczosAlgorithm = new BasicLanczosAlgorithmWithGroundState(architecture, MaxNbrIterLanczos, FastDiskFlag, ResumeFastDiskFlag);
		else
		  this->LanczosAlgorithm = new BasicLanczosAlgorithm(architecture, NbrEigenvalue, MaxNbrIterLanczos);
	      else
		if ((EvaluateEigenvectors == true) || (forceEigenstateComputation == true))
		  this->LanczosAlgorithm = new BasicLanczosAlgorithmWithGroundStateDiskStorage(architecture, NbrIterLanczos, MaxNbrIterLanczos);
		else
		  this->LanczosAlgorithm = new BasicLanczosAlgorithmWithDiskStorage(architecture, NbrEigenvalue, MaxNbrIterLanczos);
	    }
	  else
	    {
	      if (DiskFlag == false)
		{
		  if (BlockLanczosFlag == true)
		    this->LanczosAlgorithm = new FullReorthogonalizedBlockLanczosAlgorithm (architecture, NbrEigenvalue, SizeBlockLanczos, MaxNbrIterLanczos);
		  else
		    this->LanczosAlgorithm = new FullReorthogonalizedLanczosAlgorithm (architecture, NbrEigenvalue, MaxNbrIterLanczos);
		}
	      else
		this->LanczosAlgorithm = new FullReorthogonalizedLanczosAlgorithmWithDiskStorage (architecture, NbrEigenvalue, VectorMemory, MaxNbrIterLanczos);
	    }
	}
    }
  return this->LanczosAlgorithm;
}
