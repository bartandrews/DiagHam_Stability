#include "Hamiltonian/TwoDimensionalTransverseFieldIsingHamiltonian.h"
#include "Hamiltonian/TwoDimensionalTransverseFieldAnd2DTranslationIsingHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainMirrorSymmetry.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainFullAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainFullInversionAnd2DTranslation.h"

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
  OptionManager Manager ("TwoDimensionalTransverseFieldIsingModel" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "jz-value", "coupling constant in the z direction between neighboring sites", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hx-value", "amplitude of the Zeeman term along the x axis", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hz-value", "amplitude of the Zeeman term along the z axis", 0.0);
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new SingleStringOption ('\n', "selected-points", "provide a two column ascii file that indicates which momentum sectors have to be computed");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-momentum", "disable momentum quantum numbers even if the system is translation invariant");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-inversion", "disable the inversion symmetry quantum number");
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-hxvalue", "amplitude of the random Zeeman term along the x direction on each site", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-gaussianhxvalue", "amplitude of the random Zeeman term along the x direction on each site, using a gaussian disrtibution with zero mean value and a given standard deviation", 0.0);
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
      cout << "see man page for option syntax or type TwoDimensionalTransverseFieldIsingModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSitesX = Manager.GetInteger("nbr-sitex");
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSpins = NbrSitesX * NbrSitesY;
  
  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  char* BoundaryName = new char [16];
  bool RandomFieldFlag = false;
  if (Manager.GetString("fullh-values") != 0)
    RandomFieldFlag = true;
  if (Manager.GetBoolean("use-periodic") == false)
    sprintf (BoundaryName, "open");
  else
    sprintf (BoundaryName, "closed");
  sprintf (OutputFileName, "spin_1_2_ising_transversefield_%s_n_%d_x_%d_y_%d", BoundaryName, (NbrSitesX * NbrSitesY), NbrSitesX, NbrSitesY);
  sprintf (CommentLine, " ising with %s boundary conditions and %d sites in the x direction, %d sites in the y direction \n#", BoundaryName, NbrSitesX, NbrSitesY);

  char* OutputParameterFileName = new char [256];
  sprintf (OutputParameterFileName, "jz_%.6f", Manager.GetDouble("jz-value"));    
  if ((Manager.GetDouble("hx-value") != 0.0) || (Manager.GetDouble("hz-value") != 0.0))
    {
      sprintf (OutputParameterFileName + strlen(OutputParameterFileName), "_hx_%.6f_hz_%.6f", Manager.GetDouble("hx-value"), Manager.GetDouble("hz-value"));
    }
  if ((Manager.GetDouble("random-hxvalue") != 0.0) || (Manager.GetDouble("random-hzvalue") != 0.0))
    {
      RandomFieldFlag = true;
      sprintf (OutputParameterFileName + strlen(OutputParameterFileName), "_randomhx_%.6f_randomhz_%.6f_runid_%ld", Manager.GetDouble("random-hxvalue"), 
	       Manager.GetDouble("random-hzvalue"), Manager.GetInteger("run-id"));
    }
  else
    { 
      if ((Manager.GetDouble("random-gaussianhxvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  RandomFieldFlag = true;
	  sprintf (OutputParameterFileName + strlen(OutputParameterFileName), "_grandomhx_%.6f_grandomhz_%.6f_runid_%ld", Manager.GetDouble("random-gaussianhxvalue"), 
		   Manager.GetDouble("random-gaussianhzvalue"), Manager.GetInteger("run-id"));
	}
    }
  if ((Manager.GetBoolean("use-periodic") == true) && (RandomFieldFlag == false) && (Manager.GetBoolean("disable-momentum") == false))
    Lanczos.SetComplexAlgorithms();
  
  char* FullOutputFileName = new char [strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
  sprintf (FullOutputFileName, "%s_%s.dat", OutputFileName, OutputParameterFileName);
  double JzValue = Manager.GetDouble("jz-value");
  double** HxValues = new double* [NbrSitesX];
  double** HzValues = new double* [NbrSitesX];
  for (int i = 0; i < NbrSitesX; ++i)
    {
      HxValues[i] = new double [NbrSitesY];
      HzValues[i] = new double [NbrSitesY];
      for (int j = 0; j < NbrSitesY; ++j)
	{
	  HxValues[i][j] = Manager.GetDouble("hx-value");
	  HzValues[i][j] = Manager.GetDouble("hz-value");
	}
    }
  if (Manager.GetString("fullh-values") != 0)
    {
      MultiColumnASCIIFile HFieldFile;
      if (HFieldFile.Parse(Manager.GetString("fullh-values")) == false)
	{
	  HFieldFile.DumpErrors(cout);
	  return -1;
	}
      if (HFieldFile.GetNbrLines() >= NbrSpins)
	{
	  double* TmpHxValues = HFieldFile.GetAsDoubleArray(2);
	  double* TmpHzValues = HFieldFile.GetAsDoubleArray(3);
	  int Index = 0;
	  for (int i = 0; i < NbrSitesX; ++i)
	    {
	      for (int j = 0; j < NbrSitesY; ++j)
		{
		  HxValues[i][j] = TmpHxValues[Index];
		  HzValues[i][j] = TmpHzValues[Index];
		  ++Index;
		}
	    }
	  delete[] TmpHxValues;
	  delete[] TmpHzValues;
	  if (HFieldFile.GetNbrLines() > NbrSpins)
	    {
	      cout << "warning, " << Manager.GetString("fullh-values") << " has more hz values than the number of sites" << endl;
	    }
	}
      else
	{
	  cout << "error, " << Manager.GetString("fullh-values") << " has less hz values than the number of sites" << endl;
	  return 0;
	}
    }
  else
    {
      if ((Manager.GetDouble("random-hxvalue") != 0.0) || (Manager.GetDouble("random-gaussianhxvalue") != 0.0))
	{
	  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	  RandomNumber->UseTimeSeed();
	  for (int i = 0; i < NbrSitesX; ++i)
	    {
	      for (int j = 0; j < NbrSitesY; ++j)
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
		  HxValues[i][j] += Tmp;
		}
	    }
	}
      if ((Manager.GetDouble("random-hzvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	  RandomNumber->UseTimeSeed();
	  for (int i = 0; i < NbrSitesX; ++i)
	    {
	      for (int j = 0; j < NbrSitesY; ++j)
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
		  HxValues[i][j] += Tmp;
		}
	    }
	}
      if ((Manager.GetDouble("random-hxvalue") != 0.0) || (Manager.GetDouble("random-hzvalue") != 0.0) ||
	  (Manager.GetDouble("random-gaussianhxvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0))
	{
	  if (Architecture.GetArchitecture()->CanWriteOnDisk())
	    {
	      char* HOutputFileName = new char [strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
	      sprintf (HOutputFileName, "%s_%s.hvalues", OutputFileName, OutputParameterFileName);
	      ofstream File;
	      File.open(HOutputFileName, ios::binary | ios::out); 
	      File.precision(14); 
	      for (int i = 0; i < NbrSitesX; ++i)
		{
		  for (int j = 0; j < NbrSitesY; ++j)
		    {
		      File << i << " " << j << " " << HxValues[i][j] << " " << HzValues[i][j] << endl;
		    }
		}
	      File.close();	      
	    }
	}
      if ((Architecture.GetArchitecture()->GetArchitectureID() & AbstractArchitecture::SimpleMPI) != 0)
	{
	  SimpleMPIArchitecture* TmpArchitecture = (SimpleMPIArchitecture*) Architecture.GetArchitecture();
	  for (int i = 0; i < NbrSitesX; ++i)
	    {
	      TmpArchitecture->BroadcastToSlaves(HxValues[i], NbrSitesY);
	      TmpArchitecture->BroadcastToSlaves(HzValues[i], NbrSitesY);
	    }
	}
    }

  if (((Manager.GetBoolean("use-periodic") == true) && (RandomFieldFlag == false)) && (Manager.GetBoolean("disable-momentum") == false))
    {
      int NbrMomenta = 0;
      int* XMomenta = 0;
      int* YMomenta = 0;
      if (Manager.GetString("selected-points") == 0)
	{
	  int MinXMomentum = Manager.GetInteger("only-kx");
	  int MaxXMomentum = MinXMomentum;
	  if (MinXMomentum < 0)
	    {
	      MinXMomentum = 0;
	      MaxXMomentum = NbrSitesX - 1;
	    }
	  int MinYMomentum = Manager.GetInteger("only-ky");
	  int MaxYMomentum = MinYMomentum;
	  if (MinYMomentum < 0)
	    {
	      MinYMomentum = 0;
	      MaxYMomentum = NbrSitesY - 1;
	    }
	  NbrMomenta = (MaxXMomentum - MinXMomentum + 1) * (MaxYMomentum - MinYMomentum + 1);
	  XMomenta = new int[NbrMomenta];
	  YMomenta = new int[NbrMomenta];
	  int TmpIndex = 0;
	  for (int i = MinXMomentum; i <= MaxXMomentum; ++i)
	    {
	      for (int j = MinYMomentum; j <= MaxYMomentum; ++j)
		{
		  XMomenta[TmpIndex] = i;
		  YMomenta[TmpIndex] = j;
		  ++TmpIndex;
		}
	    }
	}
      else
	{
	  MultiColumnASCIIFile MomentumFile;
	  if (MomentumFile.Parse(Manager.GetString("selected-points")) == false)
	    {
	      MomentumFile.DumpErrors(cout);
	      return -1;
	    }
	  NbrMomenta = MomentumFile.GetNbrLines();
	  XMomenta = MomentumFile.GetAsIntegerArray(0);
	  YMomenta = MomentumFile.GetAsIntegerArray(1);
	}
      bool FirstRun = true;
      for (int MomentumSector = 0; MomentumSector < NbrMomenta; ++MomentumSector)
	{
	  int TmpInvXMomenta = (NbrSitesX - XMomenta[MomentumSector]) % NbrSitesX;
	  int TmpInvYMomenta = (NbrSitesY - YMomenta[MomentumSector]) % NbrSitesY;
	  int MaxInversionSector = -1;
	  if ((TmpInvXMomenta == XMomenta[MomentumSector]) && (TmpInvYMomenta == YMomenta[MomentumSector]) && (Manager.GetBoolean("disable-inversion") == false))
	    {
	      MaxInversionSector = 1;
	    }
	  for (int InversionSector = -1; InversionSector <= MaxInversionSector; InversionSector += 2)
	    {
	      cout << "-------------------------------------------" << endl;
	      cout << "kx=" << XMomenta[MomentumSector] << "  ky=" << YMomenta[MomentumSector] << endl; 
	      AbstractSpinChain* Chain = 0;
	      if (MaxInversionSector == -1)
		Chain = new Spin1_2ChainFullAnd2DTranslation (XMomenta[MomentumSector], NbrSitesX, YMomenta[MomentumSector], NbrSitesY);
	      else
		Chain = new Spin1_2ChainFullInversionAnd2DTranslation (InversionSector, XMomenta[MomentumSector], NbrSitesX, YMomenta[MomentumSector], NbrSitesY);
	      if (Chain->GetHilbertSpaceDimension() > 0)
		{
		  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
		  // 	      for (int i = 0; i < Chain->GetHilbertSpaceDimension(); ++i)
		  // 		Chain->PrintState(cout, i) << endl;
		  TwoDimensionalTransverseFieldAnd2DTranslationIsingHamiltonian* Hamiltonian = 0;
		  Hamiltonian = new TwoDimensionalTransverseFieldAnd2DTranslationIsingHamiltonian(Chain, XMomenta[MomentumSector], NbrSitesX, 
												  YMomenta[MomentumSector], NbrSitesY, 
												  JzValue, HxValues[0][0], HzValues[0][0]);
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
		  if (MaxInversionSector == -1)
		    sprintf (TmpEigenstateString, "%s_%s_kx_%d_ky_%d", OutputFileName, OutputParameterFileName, XMomenta[MomentumSector], YMomenta[MomentumSector]);
		  else
		    sprintf (TmpEigenstateString, "%s_%s_kx_%d_ky_%d_i_%d", OutputFileName, OutputParameterFileName, XMomenta[MomentumSector], YMomenta[MomentumSector], InversionSector);
		  
		  char* TmpString = new char[16];
		  if (Manager.GetBoolean("disable-inversion") == true)
		    {
		      sprintf (TmpString, "%d %d", XMomenta[MomentumSector], YMomenta[MomentumSector]);
		    }
		  else
		    {
		      if (MaxInversionSector == -1)
			{
			  sprintf (TmpString, "%d %d 0", XMomenta[MomentumSector], YMomenta[MomentumSector]);
			}
		      else
			{
			  sprintf (TmpString, "%d %d %d", XMomenta[MomentumSector], YMomenta[MomentumSector], InversionSector);
			}
		    }
		  GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
					      FirstRun, TmpEigenstateString);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  FirstRun = false;
		  delete Hamiltonian;
		  delete[] TmpString;
		  delete[] TmpEigenstateString;
		}
	      delete Chain;
	    }
	}
     
    }
  else
    {
      bool FirstRun = true;
      AbstractSpinChain* Chain = new Spin1_2ChainFull (NbrSpins);
      if (Chain->GetHilbertSpaceDimension() > 0)
	{
	  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
	  TwoDimensionalTransverseFieldIsingHamiltonian* Hamiltonian = 0;
	  Hamiltonian = new TwoDimensionalTransverseFieldIsingHamiltonian(Chain, NbrSitesX, NbrSitesY, JzValue, HxValues, HzValues, Manager.GetBoolean("use-periodic"));
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
	  sprintf (TmpEigenstateString, "%s_%s", OutputFileName, OutputParameterFileName);
	  char* TmpString = new char[1];
	  TmpString[0] = '\0';
	  GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;
	  delete Hamiltonian;
	  delete[] TmpString;
	  delete[] TmpEigenstateString;
	}
      delete Chain;
    }

  delete[] OutputFileName;
  delete[] CommentLine;
  delete[] HxValues;
  delete[] HzValues;
  delete[] FullOutputFileName;
  return 0;
}
