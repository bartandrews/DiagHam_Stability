#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#ifdef HAVE_GSL  
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

#include "Hamiltonian/SpinChainWithDisorderHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("GenericOpenSpinChainWithDisorder" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 1);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  SingleDoubleOption ('j', "j-value", "isotropic coupling constant value", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('z', "djz-value", "delta compare to the coupling constant value along z", 0.0);
  (*SystemGroup) += new  BooleanOption ('\n', "random-field", "assume random field sum_i h_i S_i^z");
  (*SystemGroup) += new  SingleDoubleOption ('w', "randomfield-uniform", "random field h_i is drawn from uniform distribution [-W,W]", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('g', "randomfield-gaussian", "random field h_i is drawn from Gaussian distribution with a given variance", 0.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "disorder-file", "file describing the disorder potential");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericOpenSpinChain -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int SpinValue = Manager.GetInteger("spin");
  int NbrSpins = Manager.GetInteger("nbr-spin");

  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  if ((SpinValue & 1) == 0)
    {
      sprintf (OutputFileName, "spin_%d_openchain_n_%d", (SpinValue / 2), NbrSpins);
      sprintf (CommentLine, " open spin %d chain with %d sites \n# 2Sz ", (SpinValue / 2), NbrSpins);
    }
  else
    {
      sprintf (OutputFileName, "spin_%d_2_openchain_n_%d", SpinValue, NbrSpins);
      sprintf (CommentLine, " open spin %d/2 chain with %d sites \n# 2Sz", SpinValue, NbrSpins);
    }
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);
  double* JValues = new double [NbrSpins - 1];
  JValues[0] = Manager.GetDouble("j-value");
  for (int i = 1; i < (NbrSpins - 1); ++i)
    JValues[i] = JValues[0];
  double* JzValues = new double [NbrSpins - 1];
  double TmpDeltaJz = Manager.GetDouble("djz-value");
  for (int i = 0; i < (NbrSpins - 1); ++i)
    JzValues[i] = JValues[i] + TmpDeltaJz;

  double* HValues = new double [NbrSpins];
  int TmpNbrDisorderTerms;

#ifndef HAVE_GSL  
  cout << "GSL required. " << endl;
  exit(1);
#endif

  if (Manager.GetString("disorder-file") == 0)
    {
#ifdef HAVE_GSL  
       const gsl_rng_type * T;
       gsl_rng * rng;
  
       // select random number generator 
       rng = gsl_rng_alloc (gsl_rng_mt19937);

       long seed = time (NULL) * getpid();
       gsl_rng_set (rng, seed);

      for (int i = 0; i < NbrSpins; i++)
       {
         if (Manager.GetBoolean("random-field"))
           {
             if (fabs(Manager.GetDouble("randomfield-uniform")) > 0.0)
               {
                  HValues[i] = gsl_ran_flat (rng, -fabs(Manager.GetDouble("randomfield-uniform")), fabs(Manager.GetDouble("randomfield-uniform")));
               }
             else if (fabs(Manager.GetDouble("randomfield-gaussian")) > 0.0)  
              {
                 HValues[i] = gsl_ran_gaussian (rng, fabs(Manager.GetDouble("randomfield-gaussian")));
              } 
          }
         else
           HValues[i] = 0.0;
         cout<<"h["<<i<<"]="<<HValues[i]<<endl;
        } 

      gsl_rng_free (rng);
#endif
    }
  else
    {
      ConfigurationParser DisorderDefinition;
      if (DisorderDefinition.Parse(Manager.GetString("disorder-file")) == false)
	{
	  DisorderDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      if (DisorderDefinition.GetAsDoubleArray("Disorder", ' ', HValues, TmpNbrDisorderTerms) == true)
	{
	  if (TmpNbrDisorderTerms != NbrSpins)
	    {
	      cout << "Invalid number of disorder terms" << endl;
	      return -1;	  
	    }
	}

      for (int i = 0; i < NbrSpins; i++)
        cout<<"h["<<i<<"]="<<HValues[i]<<endl;
    }

  int MaxSzValue = NbrSpins * SpinValue;
  int InitalSzValue = MaxSzValue & 1;
  if (Manager.GetInteger("initial-sz") > 1)
    {
      InitalSzValue += (Manager.GetInteger("initial-sz") & ~1);
    }
  if (Manager.GetInteger("nbr-sz") > 0)
    {
      MaxSzValue = InitalSzValue + ((Manager.GetInteger("nbr-sz") - 1) * 2);
    }
  bool FirstRun = true;
  for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
    {
      AbstractSpinChain* Chain = 0;
      switch (SpinValue)
	{
	case 1 :
	  Chain = new Spin1_2Chain (NbrSpins, InitalSzValue, 1000000);
	  break;
	case 2 :
	  Chain = new Spin1Chain (NbrSpins, InitalSzValue, 1000000);
	  break;
	default :
	  {
	    if ((SpinValue & 1) == 0)
	      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
	    else 
	      cout << "spin " << SpinValue << "/2 are not available" << endl;
	    return -1;
	  }
	}

      SpinChainWithDisorderHamiltonian Hamiltonian (Chain, NbrSpins, JValues, JzValues, HValues);

      char* TmpSzString = new char[64];
      sprintf (TmpSzString, "%d", InitalSzValue);
      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
      sprintf (TmpEigenstateString, "%s_sz_%d", OutputFileName, InitalSzValue);

      GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
			       FirstRun, TmpEigenstateString);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      FirstRun = false;
      delete Chain;
      delete[] TmpSzString;
      delete[] TmpEigenstateString;
    }
  delete[] OutputFileName;
  delete[] CommentLine;
  delete[] JValues;
  delete[] FullOutputFileName;
  delete[] JzValues;
  delete[] HValues; 
  return 0;
}
