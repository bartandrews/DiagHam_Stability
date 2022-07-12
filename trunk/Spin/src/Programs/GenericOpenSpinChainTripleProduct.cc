#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/SpinChainTripleProductHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainParitySymmetry.h"
#include "HilbertSpace/Spin1Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("GenericOpenSpinChainTripleProduct" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);

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
  (*SystemGroup) += new  SingleDoubleOption ('c', "chi-value", "coupling constant  in front of the triple product S_i . (S_{i+1} x S_{i+2})", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hz-value", "amplitude of the Zeeman term along the z axis", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-hzvalue", "amplitude of the random Zeeman term on each site", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-gaussianhzvalue", "amplitude of the random Zeeman term on each site, using a gaussian disrtibution with zero mean value and a given standard deviation", 0.0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "run-id", "add an additional run id to the file name when using the --random-hzvalue option", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "fullhz-values", "name of the file that contains the Zeeman term amplitudes for each site");
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  BooleanOption ('\n', "use-parity", "use the parity symmetry");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericOpenSpinChainTripleProduct -h" << endl;
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
  char* BoundaryName = new char [16];
  if (Manager.GetBoolean("use-periodic") == false)
    sprintf (BoundaryName, "open");
  else
    sprintf (BoundaryName, "closed");
  if ((SpinValue & 1) == 0)
    {
      sprintf (OutputFileName, "spin_%d_%schain_n_%d", (SpinValue / 2), BoundaryName, NbrSpins);
      if (Manager.GetBoolean("use-parity") == true)
	{
	  sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz P_parity ", BoundaryName, (SpinValue / 2), NbrSpins);
	}
      else
	{
	  sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz ", BoundaryName, (SpinValue / 2), NbrSpins);
	}
    }
  else
    {
      sprintf (OutputFileName, "spin_%d_2_%schain_n_%d", SpinValue, BoundaryName, NbrSpins);
      if (Manager.GetBoolean("use-parity") == true)
	{
	  sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz P_parity", BoundaryName, SpinValue, NbrSpins);
	}
      else
	{
	  sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz P_parity", BoundaryName, SpinValue, NbrSpins);
	}
    }
  char* OutputParameterFileName = new char [256];
  if (Manager.GetDouble("djz-value") == 0.0)
    {      
      if ((Manager.GetDouble("hz-value") == 0.0) && (Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	{
	  sprintf (OutputParameterFileName, "j_%.6f_chi_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("chi-value"));
	}
      else
	{
	  if ((Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	    {
	      sprintf (OutputParameterFileName, "j_%.6f_chi_%.6f_hz_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("chi-value"), Manager.GetDouble("hz-value"));
	    }
	  else
	    {
	      if (Manager.GetDouble("random-gaussianhzvalue") == 0.0)
		{
		  sprintf (OutputParameterFileName, "j_%.6f_chi_%.6f_hz_%.6f_randomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), Manager.GetDouble("chi-value"), 
			   Manager.GetDouble("hz-value"), Manager.GetDouble("random-hzvalue"), Manager.GetInteger("run-id"));
		}
	      else
		{
		  sprintf (OutputParameterFileName, "j_%.6f_chi_%.6f_hz_%.6f_gaussianrandomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), Manager.GetDouble("chi-value"), 
			   Manager.GetDouble("hz-value"), Manager.GetDouble("random-gaussianhzvalue"), Manager.GetInteger("run-id"));
		}
	    }
	}
    }
  else
    {
      if ((Manager.GetDouble("hz-value") == 0.0) && (Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	{
	  sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_chi_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("djz-value"), Manager.GetDouble("chi-value"));
	}
      else
	{
	  if ((Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	    {
	      sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_chi_%.6f_hz_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("djz-value"), 
		       Manager.GetDouble("chi-value"), Manager.GetDouble("hz-value"));
	    }
	  else
	    {
	      if (Manager.GetDouble("random-gaussianhzvalue") == 0.0)
		{
		  sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_chi_%.6f_hz_%.6f_randomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), 
			   Manager.GetDouble("djz-value"), Manager.GetDouble("chi-value"), 
			   Manager.GetDouble("hz-value"), Manager.GetDouble("random-hzvalue"), Manager.GetInteger("run-id"));
		}
	      else
		{
		  sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_chi_%.6f_hz_%.6f_gaussianrandomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), 
			   Manager.GetDouble("djz-value"), Manager.GetDouble("chi-value"), 
			   Manager.GetDouble("hz-value"), Manager.GetDouble("random-gaussianhzvalue"), Manager.GetInteger("run-id"));
		}
	    }
	}
    }
    
  char* FullOutputFileName = new char [strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
  sprintf (FullOutputFileName, "%s_%s.dat", OutputFileName, OutputParameterFileName);
  double* JValues = new double [NbrSpins];
  JValues[0] = Manager.GetDouble("j-value");
  for (int i = 1; i < NbrSpins; ++i)
    JValues[i] = JValues[0];
  double* JzValues = new double [NbrSpins];
  double TmpDeltaJz = Manager.GetDouble("djz-value");
  for (int i = 0; i < NbrSpins; ++i)
    JzValues[i] = JValues[i] + TmpDeltaJz;
  double* ChiValues = new double [NbrSpins];
  ChiValues[0] = Manager.GetDouble("chi-value");
  for (int i = 1; i < NbrSpins; ++i)
    ChiValues[i] = ChiValues[0];
  double* HzValues = 0;
  if (Manager.GetString("fullhz-values") != 0)
    {
      MultiColumnASCIIFile HFieldFile;
      if (HFieldFile.Parse(Manager.GetString("fullhz-values")) == false)
	{
	  HFieldFile.DumpErrors(cout);
	  return -1;
	}
      if (HFieldFile.GetNbrLines() == NbrSpins)
	{
	  HzValues = HFieldFile.GetAsDoubleArray(0);
	}
      else
	{
	  if (HFieldFile.GetNbrLines() > NbrSpins)
	    {
	      cout << "warning, " << Manager.GetString("fullhz-values") << " has more hz values than the number of sites" << endl;
	      HzValues = HFieldFile.GetAsDoubleArray(0);
	    }
	  else
	    {
	      cout << "error, " << Manager.GetString("fullhz-values") << " has less hz values than the number of sites" << endl;
	      return 0;
	    }	  
	}
    }
  else
    {
      if ((Manager.GetDouble("hz-value") != 0.0) || (Manager.GetDouble("random-hzvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  HzValues = new double [NbrSpins];
	  HzValues[0] = Manager.GetDouble("hz-value");
	  for (int i = 1; i < NbrSpins; ++i)
	    HzValues[i] = HzValues[0];
	  if ((Manager.GetDouble("random-hzvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	    {
	      AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	      RandomNumber->UseTimeSeed();
	      char* HzOutputFileName = new char [strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
	      sprintf (HzOutputFileName, "%s_%s.hzvalues", OutputFileName, OutputParameterFileName);
	      ofstream File;
	      File.open(HzOutputFileName, ios::binary | ios::out); 
	      File.precision(14); 
	      for (int i = 0; i < NbrSpins; ++i)
		{
		  double Tmp;
		  if (Manager.GetDouble("random-hzvalue") != 0.0)
		    {
		      Tmp = Manager.GetDouble("random-hzvalue") * (2.0 * RandomNumber->GetRealRandomNumber() - 1.0);
		    }
		  else
		    {
		      Tmp = RandomNumber->GetGaussianRandomNumber(0.0, Manager.GetDouble("random-gaussianhzvalue"));
		    }
		  HzValues[i] += Tmp;
		  File << Tmp << endl;
		}
	      File.close();
	      
	    }
      
	}
    }
  int MaxSzValue = NbrSpins * SpinValue;
  int InitalSzValue = MaxSzValue & 1;
  if  ((Manager.GetDouble("hz-value") != 0) || (Manager.GetDouble("random-hzvalue") != 0.0)|| (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
    InitalSzValue = -MaxSzValue;
  if ((Manager.GetInteger("initial-sz") != 0) || (Manager.GetInteger("nbr-sz") > 0))
    {
      InitalSzValue = Manager.GetInteger("initial-sz");
    }
  if (Manager.GetInteger("nbr-sz") > 0)
    {
      MaxSzValue = InitalSzValue + ((Manager.GetInteger("nbr-sz") - 1) * 2);
    }
  bool FirstRun = true;
  for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
    {
      int MaxParity = 0;
      if ((InitalSzValue == 0) && (Manager.GetBoolean("use-parity") == true))
	{
	  MaxParity = 1;
	}
      for (int Parity = 0; Parity <= MaxParity; ++Parity)
	{ 
	  AbstractSpinChain* Chain = 0;
	  switch (SpinValue)
	    {
	    case 1 :
	      {
		if (MaxParity == 1)
		  {
		    Chain = new Spin1_2ChainParitySymmetry (NbrSpins, Parity, 1000000);
		  }
		else
		  {
		    Chain = new Spin1_2ChainNew (NbrSpins, InitalSzValue, 1000000);
		  }
	      }
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
	  
	  SpinChainTripleProductHamiltonian* Hamiltonian = 0;
	  if (HzValues == 0)
	    Hamiltonian = new SpinChainTripleProductHamiltonian(Chain, NbrSpins, JValues, JzValues, ChiValues, Manager.GetBoolean("use-periodic"));
	  else
	    Hamiltonian = new SpinChainTripleProductHamiltonian(Chain, NbrSpins, JValues, JzValues, ChiValues, HzValues, Manager.GetBoolean("use-periodic"));
	  char* TmpSzString = new char[64];
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
	  if (Manager.GetBoolean("use-parity") == true)
	    {
	      sprintf (TmpSzString, "%d %d", InitalSzValue, Parity);
	    }
	  else
	    {
	      sprintf (TmpSzString, "%d", InitalSzValue);
	    }
	  if (MaxParity == 1)
	    {
	      sprintf (TmpEigenstateString, "%s_%s_sz_%d_p_%d", OutputFileName, OutputParameterFileName, InitalSzValue, Parity);
	    }
	  else
	    {
	      sprintf (TmpEigenstateString, "%s_%s_sz_%d", OutputFileName, OutputParameterFileName, InitalSzValue);
	    }
	  GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
				      FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;
	  delete Hamiltonian;
	  delete Chain;
	  delete[] TmpSzString;
	  delete[] TmpEigenstateString;
	}
    }
  delete[] OutputFileName;
  delete[] CommentLine;
  delete[] JValues;
  delete[] FullOutputFileName;
  delete[] JzValues;
  return 0;
}
