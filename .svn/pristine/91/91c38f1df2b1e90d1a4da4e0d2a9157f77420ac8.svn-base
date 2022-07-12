#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/SpinChainHamiltonian.h"
#include "Hamiltonian/SpinChainFullHamiltonian.h"
#include "Hamiltonian/SpinChainRealFullHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainMirrorSymmetry.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"
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
  OptionManager Manager ("FullGenericOpenSpinChain" , "0.01");
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
  (*SystemGroup) += new  SingleDoubleOption ('x', "jx-value", "coupling constant in the x direction between neighboring sites", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('y', "jy-value", "coupling constant in the x direction between neighboring sites", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('z', "jz-value", "coupling constant in the x direction between neighboring sites", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hx-value", "amplitude of the Zeeman term along the x axis", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hy-value", "amplitude of the Zeeman term along the y axis", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hz-value", "amplitude of the Zeeman term along the z axis", 0.0);
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  BooleanOption ('\n', "use-mirror", "use the mirror symmetry");
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-hxvalue", "amplitude of the random Zeeman term along the x direction on each site", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-gaussianhxvalue", "amplitude of the random Zeeman term along the x direction on each site, using a gaussian disrtibution with zero mean value and a given standard deviation", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-hyvalue", "amplitude of the random Zeeman term along the y direction on each site", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-gaussianhyvalue", "amplitude of the random Zeeman term along the y direction on each site, using a gaussian disrtibution with zero mean value and a given standard deviation", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-hzvalue", "amplitude of the random Zeeman term along the z direction on each site", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-gaussianhzvalue", "amplitude of the random Zeeman term along the z direction on each site, using a gaussian disrtibution with zero mean value and a given standard deviation", 0.0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "run-id", "add an additional run id to the file name when using the --random-hzvalue option", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "fullh-values", "name of the file that contains the Zeeman term amplitudes for each site");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FullGenericOpenSpinChain -h" << endl;
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
      sprintf (CommentLine, " %s spin %d chain with %d sites \n#", BoundaryName, (SpinValue / 2), NbrSpins);
    }
  else
    {
      sprintf (OutputFileName, "spin_%d_2_%schain_n_%d", SpinValue, BoundaryName, NbrSpins);
      sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n#", BoundaryName, SpinValue, NbrSpins);
    }
  char* OutputParameterFileName = new char [256];
  sprintf (OutputParameterFileName, "jx_%.6f_jy_%.6f_jz_%.6f", Manager.GetDouble("jx-value"), Manager.GetDouble("jy-value"), Manager.GetDouble("jz-value"));    
  if ((Manager.GetDouble("hx-value") != 0.0) || (Manager.GetDouble("hy-value") != 0.0) || (Manager.GetDouble("hz-value") != 0.0))
    {
      sprintf (OutputParameterFileName + strlen(OutputParameterFileName), "_hx_%.6f_hy_%.6f_hz_%.6f", Manager.GetDouble("hx-value"), 
	       Manager.GetDouble("hy-value"), Manager.GetDouble("hz-value"));
    }
  if ((Manager.GetDouble("random-hxvalue") != 0.0) || (Manager.GetDouble("random-hyvalue") != 0.0) || (Manager.GetDouble("random-hzvalue") != 0.0))
    {
      sprintf (OutputParameterFileName + strlen(OutputParameterFileName), "_randomhx_%.6f_randomhy_%.6f_randomhz_%.6f_runid_%ld", Manager.GetDouble("random-hxvalue"), 
	       Manager.GetDouble("random-hyvalue"), Manager.GetDouble("random-hzvalue"), Manager.GetInteger("run-id"));
    }
  else
    { 
      if ((Manager.GetDouble("random-gaussianhxvalue") != 0.0) || (Manager.GetDouble("random-gaussianhyvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  sprintf (OutputParameterFileName + strlen(OutputParameterFileName), "_grandomhx_%.6f_grandomhy_%.6f_grandomhz_%.6f_runid_%ld", Manager.GetDouble("random-gaussianhxvalue"), 
		   Manager.GetDouble("random-gaussianhyvalue"), Manager.GetDouble("random-gaussianhzvalue"), Manager.GetInteger("run-id"));
	}
    }
    
  char* FullOutputFileName = new char [strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
  sprintf (FullOutputFileName, "%s_%s.dat", OutputFileName, OutputParameterFileName);
  double* JxValues = new double [NbrSpins];
  JxValues[0] = Manager.GetDouble("jx-value");
  for (int i = 1; i < NbrSpins; ++i)
    JxValues[i] = JxValues[0];
  double* JyValues = new double [NbrSpins];
  JyValues[0] = Manager.GetDouble("jy-value");
  for (int i = 1; i < NbrSpins; ++i)
    JyValues[i] = JyValues[0];
  double* JzValues = new double [NbrSpins];
  JzValues[0] = Manager.GetDouble("jz-value");
  for (int i = 1; i < NbrSpins; ++i)
    JzValues[i] = JzValues[0];
  double* HxValues = 0;
  double* HyValues = 0;
  double* HzValues = 0;
  if (Manager.GetString("fullh-values") != 0)
    {
      MultiColumnASCIIFile HFieldFile;
      if (HFieldFile.Parse(Manager.GetString("fullh-values")) == false)
	{
	  HFieldFile.DumpErrors(cout);
	  return -1;
	}
      if (HFieldFile.GetNbrLines() == NbrSpins)
	{
	  HxValues = HFieldFile.GetAsDoubleArray(0);
	  HyValues = HFieldFile.GetAsDoubleArray(1);
	  HzValues = HFieldFile.GetAsDoubleArray(2);
	}
      else
	{
	  if (HFieldFile.GetNbrLines() > NbrSpins)
	    {
	      cout << "warning, " << Manager.GetString("fullh-values") << " has more hz values than the number of sites" << endl;
	      HzValues = HFieldFile.GetAsDoubleArray(0);
	    }
	  else
	    {
	      cout << "error, " << Manager.GetString("fullh-values") << " has less hz values than the number of sites" << endl;
	      return 0;
	    }	  
	}
    }
  else
    {
      HxValues = new double [NbrSpins];
      HxValues[0] = Manager.GetDouble("hx-value");
      for (int i = 1; i < NbrSpins; ++i)
	HxValues[i] = HxValues[0];
      if ((Manager.GetDouble("random-hxvalue") != 0.0) || (Manager.GetDouble("random-gaussianhxvalue") != 0.0))
	{
	  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	  RandomNumber->UseTimeSeed();
	  for (int i = 0; i < NbrSpins; ++i)
	    {
	      double Tmp;
	      if (Manager.GetDouble("random-hxvalue") != 0.0)
		{
		  Tmp = Manager.GetDouble("random-hxvalue") * (2.0 * RandomNumber->GetRealRandomNumber() - 1.0);
		}
	      else
		{
		  Tmp = RandomNumber->GetGaussianRandomNumber(0.0, Manager.GetDouble("random-gaussianhxvalue"));
		}
	      HxValues[i] += Tmp;
	    }
	}
      HyValues = new double [NbrSpins];
      HyValues[0] = Manager.GetDouble("hy-value");
      for (int i = 1; i < NbrSpins; ++i)
	HyValues[i] = HyValues[0];
      if ((Manager.GetDouble("random-hyvalue") != 0.0) || (Manager.GetDouble("random-gaussianhyvalue") != 0.0))
	{
	  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	  RandomNumber->UseTimeSeed();
	  for (int i = 0; i < NbrSpins; ++i)
	    {
	      double Tmp;
	      if (Manager.GetDouble("random-hyvalue") != 0.0)
		{
		  Tmp = Manager.GetDouble("random-hyvalue") * (2.0 * RandomNumber->GetRealRandomNumber() - 1.0);
		}
	      else
		{
		  Tmp = RandomNumber->GetGaussianRandomNumber(0.0, Manager.GetDouble("random-gaussianhyvalue"));
		}
	      HyValues[i] += Tmp;
	    }
	}
      HzValues = new double [NbrSpins];
      HzValues[0] = Manager.GetDouble("hz-value");
      for (int i = 1; i < NbrSpins; ++i)
	HzValues[i] = HzValues[0];
      if ((Manager.GetDouble("random-hzvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	  RandomNumber->UseTimeSeed();
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
	    }
	}
      if ((Manager.GetDouble("random-hxvalue") != 0.0) || (Manager.GetDouble("random-hyvalue") != 0.0) || (Manager.GetDouble("random-hzvalue") != 0.0) ||
	  (Manager.GetDouble("random-gaussianhxvalue") != 0.0) || (Manager.GetDouble("random-gaussianhyvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  if (Architecture.GetArchitecture()->CanWriteOnDisk())
	    {
	      char* HOutputFileName = new char [strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
	      sprintf (HOutputFileName, "%s_%s.hvalues", OutputFileName, OutputParameterFileName);
	      ofstream File;
	      File.open(HOutputFileName, ios::binary | ios::out); 
	      File.precision(14); 
	      for (int i = 0; i < NbrSpins; ++i)
		{
		  File << HxValues[i] << " " << HyValues[i] << " " << HzValues[i] << endl;
		}
	      File.close();	      
	    }
	}
      if ((Architecture.GetArchitecture()->GetArchitectureID() & AbstractArchitecture::SimpleMPI) != 0)
	{
	  SimpleMPIArchitecture* TmpArchitecture = (SimpleMPIArchitecture*) Architecture.GetArchitecture();
	  TmpArchitecture->BroadcastToSlaves(HxValues, NbrSpins);
	  TmpArchitecture->BroadcastToSlaves(HyValues, NbrSpins);
	  TmpArchitecture->BroadcastToSlaves(HzValues, NbrSpins);
	}
    }

  bool FirstRun = true;
  AbstractSpinChain* Chain = 0;
  switch (SpinValue)
    {
    case 1 :
      {
	Chain = new Spin1_2ChainFull (NbrSpins);
      }
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
  if (Chain->GetHilbertSpaceDimension() > 0)
    {
      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
      char* TmpEigenstateString = new char[strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
      sprintf (TmpEigenstateString, "%s_%s", OutputFileName, OutputParameterFileName);
      char* TmpString = new char[1];
      TmpString[0] = '\0';
      double TestHyValues = 0.0;
      for (int i = 0; i < NbrSpins; ++i)
	TestHyValues += fabs(HyValues[i]);
      if (TestHyValues == 0.0)
	{
	  Lanczos.SetRealAlgorithms();
	  SpinChainRealFullHamiltonian* Hamiltonian = new SpinChainRealFullHamiltonian(Chain, NbrSpins, JxValues, JyValues, JzValues, 
										       HxValues, HzValues, Manager.GetBoolean("use-periodic"));
	  GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  delete Hamiltonian;
	}
      else
	{
	  SpinChainFullHamiltonian* Hamiltonian = 0;
	  Hamiltonian = new SpinChainFullHamiltonian(Chain, NbrSpins, JxValues, JyValues, JzValues, 
						     HxValues, HyValues, HzValues, 
						     Manager.GetBoolean("use-periodic"));
	  GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
				      FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  delete Hamiltonian;
	  FirstRun = false;
	}
      delete[] TmpString;
      delete[] TmpEigenstateString;
    }
  delete Chain;

  delete[] OutputFileName;
  delete[] CommentLine;
  delete[] JxValues;
  delete[] JyValues;
  delete[] JzValues;
  delete[] FullOutputFileName;
  return 0;
}
