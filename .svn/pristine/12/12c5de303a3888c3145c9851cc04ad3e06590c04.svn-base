#include "Hamiltonian/SpinChainLongRangeXYZHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"

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
  OptionManager Manager ("SpinChainXYZ" , "0.01");
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

  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 8);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "momentum", "if non negative, only consider a given momentum sector", -1);
  (*SystemGroup) += new  SingleDoubleOption('x', "jx-value", "coupling constant along the x axis", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('y', "jy-value", "coupling constant along the y axis", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('z', "jz-value", "coupling constant along the z axis", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('\n', "powerlaw-xx", "power-law decay exponent for XX interactions", 1.0);
   (*SystemGroup) += new  SingleDoubleOption('\n', "powerlaw-yy", "power-law decay exponent for YY interactions", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('\n', "powerlaw-zz", "power-law decay exponent for ZZ interactions", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('f', "h-value", "Zeeman strength along the z axis", 0.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*OutputGroup) += new  SingleStringOption ('\n', "output-suffix", "apprend and extra suffix to the string describing the system in the output file name");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainXYZ -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSpins = Manager.GetInteger("nbr-spin");
  int InitialMomentum = 0;
  int MaxMomentum = NbrSpins;
  if ((Manager.GetInteger("momentum") >= 0) && (Manager.GetInteger("momentum") < MaxMomentum))
    {
      InitialMomentum =  Manager.GetInteger("momentum");
      MaxMomentum = InitialMomentum + 1;
    }

  double JxValue = Manager.GetDouble("jx-value");
  double JyValue = Manager.GetDouble("jy-value");
  double JzValue = Manager.GetDouble("jz-value");
  double PowerLawXX = Manager.GetDouble("powerlaw-xx");
  double PowerLawYY = Manager.GetDouble("powerlaw-yy");
  double PowerLawZZ = Manager.GetDouble("powerlaw-zz");
  double HValue = Manager.GetDouble("h-value");

  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  char* FileNamePrefix = new char [256];
  sprintf (FileNamePrefix, "spin_1_2_periodicchain");
  if (Manager.GetString("output-suffix") != 0)
    {
      char* TmpString = new char [strlen(FileNamePrefix) + strlen(Manager.GetString("output-suffix")) + 2];
      sprintf (TmpString, "%s_%s", FileNamePrefix, Manager.GetString("output-suffix"));
      delete[] FileNamePrefix;
      FileNamePrefix = TmpString;
    }
  if (HValue == 0.0)
    {
      sprintf (OutputFileName, "%s_x_%.6f_y_%.6f_z_%.6f_alphax_%.6f_alphay_%.6f_alphaz_%.6f_n_%d_k_%d", FileNamePrefix, JxValue, JyValue, JzValue, PowerLawXX, PowerLawYY, PowerLawZZ, NbrSpins, Manager.GetInteger("momentum"));
      sprintf (CommentLine, "XYZ chain with %d sites, Jx=%.6f, Jy= %.6f, Jz=%.6f, AlphaXX= %.6f, AlphaYY= %.6f, AlphaZZ= %6.f, K = %d\n", NbrSpins, JxValue, JyValue, JzValue, PowerLawXX, PowerLawYY, PowerLawZZ, Manager.GetInteger("momentum"));
    }
  else
    {
       sprintf (OutputFileName, "%s_x_%.6f_y_%.6f_z_%.6f_alphax_%.6f_alphay_%.6f_alphaz_%.6f_h_%6.f_n_%d_k_%d", FileNamePrefix, JxValue, JyValue, JzValue, PowerLawXX, PowerLawYY, PowerLawZZ, HValue, NbrSpins, Manager.GetInteger("momentum"));
      sprintf (CommentLine, "XYZ chain with %d sites, Jx=%.6f, Jy= %.6f, Jz=%.6f, AlphaXX= %.6f, AlphaYY= %.6f, AlphaZZ= %6.f, H= %6.f, K = %d\n", NbrSpins, JxValue, JyValue, JzValue, PowerLawXX, PowerLawYY, PowerLawZZ, HValue, Manager.GetInteger("momentum"));
    }
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  bool FirstRun = true;
  for (int Momentum = InitialMomentum; Momentum < MaxMomentum; ++Momentum)
	{
		Spin1_2ChainWithTranslations* Chain = new Spin1_2ChainWithTranslations (NbrSpins, Momentum, 1, 1000000, 1000000);
      
     Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());  


   		if (Chain->GetHilbertSpaceDimension() > 0)
			{
	 			 SpinChainLongRangeXYZHamiltonian Hamiltonian (Chain, NbrSpins, JxValue, JyValue, JzValue, PowerLawXX, PowerLawYY, PowerLawZZ, HValue);
	  
				char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	  			sprintf (TmpEigenstateString, "%s", OutputFileName);

          char* TmpString = new char[64];
          sprintf (TmpString, "%d ", Momentum);


	  			GenericComplexMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, TmpEigenstateString);

	  			MainTaskOperation TaskOperation (&Task);
	  			TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  			FirstRun = false;
	  			delete[] TmpEigenstateString;
			}
		 delete Chain;
	}	 
  return 0;
}
