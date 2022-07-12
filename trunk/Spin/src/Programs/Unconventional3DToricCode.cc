#include "Hamiltonian/Unconventional3DToricCodeHamiltonian.h"

#include "HilbertSpace/Spin1_2ChainFixedParity.h"
#include "HilbertSpace/Spin1_2ChainFull.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("Unconventional3DToricCode" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 2);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 2);
  (*SystemGroup) += new SingleIntegerOption  ('z', "nbr-sitez", "number of sites along the z direction", 2);
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodicx", "use periodic boundary conditionsalong the x direction");
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodicy", "use periodic boundary conditionsalong the y direction");
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodicz", "use periodic boundary conditionsalong the z direction");
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
      cout << "see man page for option syntax or type Unconventional3DToricCode -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSitesX = Manager.GetInteger("nbr-sitex");
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSitesZ = Manager.GetInteger("nbr-sitez");
  int NbrSpins = NbrSitesX * NbrSitesY * NbrSitesZ;
  
  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  char* BoundaryName = new char [64];

  if (Manager.GetBoolean("use-periodicx") == false)
    {
      if (Manager.GetBoolean("use-periodicy") == false)
	{
	  if (Manager.GetBoolean("use-periodicz") == false)
	    {
	      sprintf (BoundaryName, "openx_openy_openz");
	    }
	  else
	    {
	      sprintf (BoundaryName, "openx_openy_periodicz");
	    }
	}
      else
	{
	  if (Manager.GetBoolean("use-periodicz") == false)
	    {
	      sprintf (BoundaryName, "openx_periodicy_openz");
	    }
	  else
	    {
	      sprintf (BoundaryName, "openx_periodicy_periodicz");
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("use-periodicy") == false)
	{
	  if (Manager.GetBoolean("use-periodicz") == false)
	    {
	      sprintf (BoundaryName, "periodicx_openy_openz");
	    }
	  else
	    {
	      sprintf (BoundaryName, "periodicx_openy_periodicz");
	    }
	}
      else
	{
	  if (Manager.GetBoolean("use-periodicz") == false)
	    {
	      sprintf (BoundaryName, "periodicx_periodicy_openz");
	    }
	  else
	    {
	      sprintf (BoundaryName, "periodicx_periodicy_periodicz");
	    }
	}
    }
  sprintf (OutputFileName, "unconventional3dtoriccode_%s_n_%d_x_%d_y_%d_z_%d", BoundaryName, NbrSpins, NbrSitesX, NbrSitesY, NbrSitesZ);
  sprintf (CommentLine, " Unconventional 3D torix coce with %s boundary conditions and %d sites in the x direction, %d sites in the y direction, %d sites in the z direction \n#", BoundaryName, NbrSitesX, NbrSitesY, NbrSitesZ);

  
  char* FullOutputFileName = new char [strlen(OutputFileName) + 64];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);
  bool FirstRun = true;
  AbstractSpinChain* Chain = new Spin1_2ChainFixedParity (NbrSpins, 0);
  if (Chain->GetHilbertSpaceDimension() > 0)
    {
      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
      Unconventional3DToricCodeHamiltonian* Hamiltonian = 0;
      Hamiltonian = new Unconventional3DToricCodeHamiltonian(Chain, NbrSitesX, NbrSitesY, NbrSitesZ, 
							     Manager.GetBoolean("use-periodicx"), Manager.GetBoolean("use-periodicy"), Manager.GetBoolean("use-periodicz"));
      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
      sprintf (TmpEigenstateString, "%s", OutputFileName);
      char* TmpString = new char[16];
      sprintf (TmpString, "-1");
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
  Chain = new Spin1_2ChainFixedParity (NbrSpins, 1);
  if (Chain->GetHilbertSpaceDimension() > 0)
    {
      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
      Unconventional3DToricCodeHamiltonian* Hamiltonian = 0;
      Hamiltonian = new Unconventional3DToricCodeHamiltonian(Chain, NbrSitesX, NbrSitesY, NbrSitesZ, 
							     Manager.GetBoolean("use-periodicx"), Manager.GetBoolean("use-periodicy"), Manager.GetBoolean("use-periodicz"));
      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
      sprintf (TmpEigenstateString, "%s", OutputFileName);
      char* TmpString = new char[16];
      sprintf (TmpString, "1");
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

  delete[] OutputFileName;
  delete[] CommentLine;
  delete[] FullOutputFileName;
  return 0;
}
