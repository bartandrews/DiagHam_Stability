#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/IntegerMatrix.h"
#include "Matrix/LongIntegerMatrix.h"

#include "Hamiltonian/SpinChainHamiltonian.h"
#include "Hamiltonian/SpinChainJz2Hamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainMirrorSymmetry.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1ChainWithInversionSymmetry.h"

#include "HilbertSpace/Spin3_2Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"


#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "Options/Options.h"

#include "Operator/SpinS2Operator.h"
#include "Operator/SpinNonLocalHeisenbergOperator.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


// compute the characteristic polynomial for the real hamiltonians
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture
// sResolvedFlag = true if the characteristic polynomial should be computed per S sector
// minSValue = twice the minimum S value to consider (when sResolvedFlag is true)
// minSValue = twice the maximum S value to consider (when sResolvedFlag is true)
// sU2Degeneracy = array providing the Hilbert space dimension per S value (when sResolvedFlag is true)
// discardFourFactor = discard a global four factor used to ensure integer numbers
void SpinChainComputeCharacteristicPolynomial(SpinChainHamiltonian* hamiltonian, AbstractSpinChain* chain, char* outputFileName, AbstractArchitecture* architecture, bool sResolvedFlag, int minSValue, int maxSValue, int* sU2Degeneracy, bool discardFourFactor = false);

int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("GenericOpenSpinChain" , "0.01");
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
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hz-value", "amplitude of the Zeeman term along the z axis", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "jzjz-value", "amplitude of an additional local jz^2 term", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "jxy4-value", "amplitude of x-y couplingh term between spins distant by 4 sites", 0.0);
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  BooleanOption ('\n', "use-mirror", "use the mirror symmetry");
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-hzvalue", "amplitude of the random Zeeman term on each site", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "random-gaussianhzvalue", "amplitude of the random Zeeman term on each site, using a gaussian disrtibution with zero mean value and a given standard deviation", 0.0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "run-id", "add an additional run id to the file name when using the --random-hzvalue option", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "fullhz-values", "name of the file that contains the Zeeman term amplitudes for each site");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "export-charpolynomial", "export the hamiltonian characteristic polynomial");  
  (*ToolsGroup) += new BooleanOption  ('\n', "export-sresolvedcharpolynomial", "export the hamiltonian characteristic polynomial resolved in S quantum number");  
  (*ToolsGroup) += new BooleanOption  ('\n', "charpolynomial-nofour", "when exporting the hamiltonian characteristic polynomial, do not include a global 4 factor");  
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
  char* BoundaryName = new char [16];
  if (Manager.GetBoolean("use-periodic") == false)
    sprintf (BoundaryName, "open");
  else
    sprintf (BoundaryName, "closed");
  if ((SpinValue & 1) == 0)
    {
      sprintf (OutputFileName, "spin_%d_%schain_n_%d", (SpinValue / 2), BoundaryName, NbrSpins);
      if (Manager.GetBoolean("use-mirror") == true)
	{
	  sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz P_mirror ", BoundaryName, (SpinValue / 2), NbrSpins);
	}
      else
	{
	  sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz ", BoundaryName, (SpinValue / 2), NbrSpins);
	}
    }
  else
    {
      sprintf (OutputFileName, "spin_%d_2_%schain_n_%d", SpinValue, BoundaryName, NbrSpins);
      if (Manager.GetBoolean("use-mirror") == true)
	{
	  sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz P_mirror", BoundaryName, SpinValue, NbrSpins);
	}
      else
	{
	  sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz", BoundaryName, SpinValue, NbrSpins);
	}
    }
  char* OutputParameterFileName = new char [256];
  if (Manager.GetDouble("djz-value") == 0)
    {      
      if (Manager.GetDouble("jzjz-value") == 0.0)
	{
	  if ((Manager.GetDouble("hz-value") == 0.0) && (Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	    {
	      sprintf (OutputParameterFileName, "j_%.6f", Manager.GetDouble("j-value"));
	    }
	  else
	    {
	      if ((Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
		{
		  sprintf (OutputParameterFileName, "j_%.6f_hz_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("hz-value"));
		}
	      else
		{
		  if (Manager.GetDouble("random-gaussianhzvalue") == 0.0)
		    {
		      sprintf (OutputParameterFileName, "j_%.6f_hz_%.6f_randomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), 
			       Manager.GetDouble("hz-value"), Manager.GetDouble("random-hzvalue"), Manager.GetInteger("run-id"));
		    }
		  else
		    {
		      sprintf (OutputParameterFileName, "j_%.6f_hz_%.6f_gaussianrandomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), 
			       Manager.GetDouble("hz-value"), Manager.GetDouble("random-gaussianhzvalue"), Manager.GetInteger("run-id"));
		    }
		}
	    }
	}
      else
	{
	  if ((Manager.GetDouble("hz-value") == 0.0) && (Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	    {
	      sprintf (OutputParameterFileName, "j_%.6f_jzjz_%.6f_jxy4_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("jzjz-value"), Manager.GetDouble("jxy4-value"));
	    }
	  else
	    {
	      if ((Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
		{
		  sprintf (OutputParameterFileName, "j_%.6f_jzjz_%.6f_jxy4_%.6f_hz_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("jzjz-value"), Manager.GetDouble("jxy4-value"), Manager.GetDouble("hz-value"));
		}
	      else
		{
		  if (Manager.GetDouble("random-gaussianhzvalue") == 0.0)
		    {
		      sprintf (OutputParameterFileName, "j_%.6f_jzjz_%.6f_jxy4_%.6f_hz_%.6f_randomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), Manager.GetDouble("jzjz-value"), Manager.GetDouble("jxy4-value"), 
			       Manager.GetDouble("hz-value"), Manager.GetDouble("random-hzvalue"), Manager.GetInteger("run-id"));
		    }
		  else
		    {
		      sprintf (OutputParameterFileName, "j_%.6f_jzjz_%.6f_jxy4_%.6f_hz_%.6f_gaussianrandomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), Manager.GetDouble("jzjz-value"), Manager.GetDouble("jxy4-value"), 
			       Manager.GetDouble("hz-value"), Manager.GetDouble("random-gaussianhzvalue"), Manager.GetInteger("run-id"));
		    }
		}
	    }
	}
    }
  else
    {
      if (Manager.GetDouble("jzjz-value") == 0.0)
	{
	  if ((Manager.GetDouble("hz-value") == 0.0) && (Manager.GetDouble("random-hzvalue") == 0.0)
	      && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	    {
	      sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("djz-value"));
	    }
	  else
	    {
	      if ((Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
		{
		  sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_hz_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("djz-value"), 
			   Manager.GetDouble("hz-value"));
		}
	      else
		{
		  if (Manager.GetDouble("random-gaussianhzvalue") == 0.0)
		    {
		      sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_hz_%.6f_randomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), 
			       Manager.GetDouble("djz-value"), 
			       Manager.GetDouble("hz-value"), Manager.GetDouble("random-hzvalue"), Manager.GetInteger("run-id"));
		    }
		  else
		    {
		      sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_hz_%.6f_gaussianrandomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), 
			       Manager.GetDouble("djz-value"), 
			       Manager.GetDouble("hz-value"), Manager.GetDouble("random-gaussianhzvalue"), Manager.GetInteger("run-id"));
		    }
		}
	    }
	}
      else
	{
	  if ((Manager.GetDouble("hz-value") == 0.0) && (Manager.GetDouble("random-hzvalue") == 0.0)
	      && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
	    {
	      sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_jzjz_%.6f_jxy4_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("djz-value"), Manager.GetDouble("jzjz-value"), Manager.GetDouble("jxy4-value"));
	    }
	  else
	    {
	      if ((Manager.GetDouble("random-hzvalue") == 0.0) && (Manager.GetDouble("random-gaussianhzvalue") == 0.0))
		{
		  sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_jzjz_%.6f_jxy4_%.6f_hz_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("djz-value"), Manager.GetDouble("jzjz-value"), Manager.GetDouble("jxy4-value"), 
			   Manager.GetDouble("hz-value"));
		}
	      else
		{
		  if (Manager.GetDouble("random-gaussianhzvalue") == 0.0)
		    {
		      sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_jzjz_%.6f_jxy4_%.6f_hz_%.6f_randomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), 
			       Manager.GetDouble("djz-value"), Manager.GetDouble("jzjz-value"), Manager.GetDouble("jxy4-value"), 
			       Manager.GetDouble("hz-value"), Manager.GetDouble("random-hzvalue"), Manager.GetInteger("run-id"));
		    }
		  else
		    {
		      sprintf (OutputParameterFileName, "j_%.6f_djz_%.6f_jzjz_%.6f_jxy4_%.6f_hz_%.6f_gaussianrandomhz_%.6f_runid_%ld", Manager.GetDouble("j-value"), 
			       Manager.GetDouble("djz-value"), Manager.GetDouble("jzjz-value"), Manager.GetDouble("jxy4-value"), 
			       Manager.GetDouble("hz-value"), Manager.GetDouble("random-gaussianhzvalue"), Manager.GetInteger("run-id"));
		    }
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
  double* Jz2Values = 0;
  double* Jxy4Values = 0;
  if (Manager.GetDouble("jzjz-value") != 0.0)
    {
      Jz2Values = new double [NbrSpins];
      Jxy4Values = new double [NbrSpins];
      double TmpJz2 = Manager.GetDouble("jzjz-value");
      double TmpJxy4 = Manager.GetDouble("jxy4-value");
      for (int i = 0; i < NbrSpins; ++i)
	{
	  Jz2Values[i] = TmpJz2;
	  Jxy4Values[i] = TmpJxy4;
	}
    }
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
  if  ((Manager.GetDouble("hz-value") != 0) || (Manager.GetDouble("random-hzvalue") != 0.0) || (Manager.GetDouble("random-gaussianhzvalue") != 0.0) || (Manager.GetString("fullhz-values") != 0))
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
  if (Manager.GetBoolean("use-mirror") == true)
    {
      for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
	{
	  for (int Mirror = 0; Mirror < 2; ++Mirror)
	    { 
	      AbstractSpinChain* Chain = 0;
	      switch (SpinValue)
		{
		case 1:
		  {
		    Chain = new Spin1_2ChainMirrorSymmetry (NbrSpins, InitalSzValue, Mirror, 1000000);
		  }
		  break;
		case 2:
		  {
		    Chain = new Spin1ChainWithInversionSymmetry (NbrSpins, (2 * Mirror) - 1, InitalSzValue, 1000000);
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
		  SpinChainHamiltonian* Hamiltonian = 0;
		  if (Jz2Values == 0)
		    {
		      if (HzValues == 0)
			{
			  Hamiltonian = new SpinChainHamiltonian(Chain, NbrSpins, JValues, JzValues, Manager.GetBoolean("use-periodic"));
			}
		      else
			{
			  Hamiltonian = new SpinChainHamiltonian(Chain, NbrSpins, JValues, JzValues, HzValues, Manager.GetBoolean("use-periodic"));
			}
		    }
		  else
		    {
		      if (HzValues == 0)
			{
			  Hamiltonian = new SpinChainJz2Hamiltonian(Chain, NbrSpins, JValues, JzValues, Jz2Values, Jxy4Values, Manager.GetBoolean("use-periodic"));
			}
		      else
			{
			  Hamiltonian = new SpinChainJz2Hamiltonian(Chain, NbrSpins, JValues, JzValues, Jz2Values, Jxy4Values, HzValues, Manager.GetBoolean("use-periodic"));
			}
		    }
		  char* TmpSzString = new char[64];
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
		  sprintf (TmpSzString, "%d %d", InitalSzValue, Mirror);
		  sprintf (TmpEigenstateString, "%s_%s_sz_%d_invsym_%d", OutputFileName, OutputParameterFileName, InitalSzValue, ((2 * Mirror) - 1));

		  if (Manager.GetBoolean("export-charpolynomial"))
		    {
		      int* SU2Degeneracy = 0;
		      if (Manager.GetBoolean("export-sresolvedcharpolynomial"))
			{			   
			  SU2Degeneracy = new int [(((SpinValue * NbrSpins) - InitalSzValue) >> 1) + 1];
			  for (int k = InitalSzValue; k <= (SpinValue * NbrSpins); k += 2)
			    {
			      switch (SpinValue)
				{
				case 1 :
				  {
				    Spin1_2ChainMirrorSymmetry TmpHilbert(NbrSpins, k, Mirror, 1000000);
				    SU2Degeneracy[(k - InitalSzValue) >> 1] = TmpHilbert.GetHilbertSpaceDimension();
				  }
				  break;
				case 2 :
				  {
				    Spin1ChainWithInversionSymmetry TmpHilbert (NbrSpins, (2 * Mirror) - 1, k, 1000000);
				    SU2Degeneracy[(k - InitalSzValue) >> 1] = TmpHilbert.GetHilbertSpaceDimension();
				  }
				  break;
				}
			    }
			  for (int k = InitalSzValue; k < (SpinValue * NbrSpins); k += 2)
			    {
			      SU2Degeneracy[(k - InitalSzValue) >> 1] -= SU2Degeneracy[1 + ((k - InitalSzValue) >> 1)];
			    }
			}
		      SpinChainComputeCharacteristicPolynomial(Hamiltonian, Chain, TmpEigenstateString, Architecture.GetArchitecture(), Manager.GetBoolean("export-sresolvedcharpolynomial"), InitalSzValue, SpinValue * NbrSpins, SU2Degeneracy, Manager.GetBoolean("charpolynomial-nofour"));
		    }			      


		  
		  GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
					   FirstRun, TmpEigenstateString);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  FirstRun = false;
		  delete Hamiltonian;
		  delete[] TmpSzString;
		  delete[] TmpEigenstateString;
		}
	      delete Chain;
	    }
	}
    }
  else
    {
      for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
	{
	  AbstractSpinChain* Chain = 0;
	  switch (SpinValue)
	    {
	    case 1 :
	      {
		Chain = new Spin1_2Chain (NbrSpins, InitalSzValue, 1000000);
		//	    Chain = new Spin1_2ChainNew (NbrSpins, InitalSzValue, 1000000);
	      }
	      break;
	    case 2 :
	      Chain = new Spin1Chain (NbrSpins, InitalSzValue, 1000000);
	      break;
	    case 3 :
	      Chain = new Spin3_2Chain (NbrSpins, InitalSzValue, 1000000);
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
	  //	  SpinNonLocalHeisenbergOperator TmpOperator(Chain, NbrSpins, Manager.GetBoolean("use-periodic"));
	  //	  RealMatrix TmpOperatorMatrix (Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
	  //	  TmpOperator.GetOperator(TmpOperatorMatrix);
	  //	  cout << "O=" << endl;
	  //	  cout << TmpOperatorMatrix << endl;
	  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());		  
	  SpinChainHamiltonian* Hamiltonian = 0;
	  if (Jz2Values == 0)
	    {
	      if (HzValues == 0)
		{
		  Hamiltonian = new SpinChainHamiltonian(Chain, NbrSpins, JValues, JzValues, Manager.GetBoolean("use-periodic"));
		}
	      else
		{
		  Hamiltonian = new SpinChainHamiltonian(Chain, NbrSpins, JValues, JzValues, HzValues, Manager.GetBoolean("use-periodic"));
		}
	    }
	  else
	    {
	      if (HzValues == 0)
		{
		  Hamiltonian = new SpinChainJz2Hamiltonian(Chain, NbrSpins, JValues, JzValues, Jz2Values, Jxy4Values, Manager.GetBoolean("use-periodic"));
		}
	      else
		{
		  Hamiltonian = new SpinChainJz2Hamiltonian(Chain, NbrSpins, JValues, JzValues, Jz2Values, Jxy4Values, HzValues, Manager.GetBoolean("use-periodic"));
		}
	    }

	  // RealMatrix TmpHamiltonian (Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
	  // Hamiltonian->GetHamiltonian(TmpHamiltonian);

	  //	  cout << "H=" << endl;
	  //	  cout << TmpHamiltonian << endl;	  
	  
	  // RealMatrix Tmp1 = TmpHamiltonian * TmpOperatorMatrix;
	  //	  cout << Tmp1;
	  //	  cout << endl;
	  // RealMatrix Tmp2 = TmpOperatorMatrix * TmpHamiltonian;
	  //	  cout << Tmp2;
	  // cout << "[H,O]=" << endl;
	  // Tmp1 -= Tmp2;
	  //	  cout << Tmp1;
	  
	  // SpinS2Operator S2Operator(Chain, Chain->GetSpinChainLength());
	  // RealMatrix TmpS2Operator (Chain->GetHilbertSpaceDimension(), Chain->GetHilbertSpaceDimension(), true);
	  // S2Operator.GetOperator(TmpS2Operator);
	  // Tmp1 = TmpS2Operator * TmpOperatorMatrix;
	  // Tmp2 = TmpOperatorMatrix * TmpS2Operator;
	  // cout << "[S^2,O]=" << endl;
	  // Tmp1 -= Tmp2;
	  //	  cout << Tmp1;
	  
	  char* TmpSzString = new char[64];
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
	  sprintf (TmpSzString, "%d", InitalSzValue);
	  sprintf (TmpEigenstateString, "%s_%s_sz_%d", OutputFileName, OutputParameterFileName, InitalSzValue);

	  
	  GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName, FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;

	  // RealDiagonalMatrix TmpSpectrum = Task.GetEigenvalues();
	  // for (int k = 0; k < Chain->GetHilbertSpaceDimension();)
	  //   {
	  //     int l = k + 1;
	  //     while ((l < Chain->GetHilbertSpaceDimension()) && (fabs(TmpSpectrum[k] - TmpSpectrum[l]) < 1e-10))
	  // 	{
	  // 	  l++;
	  // 	}
	  //     if ((l - k) > 1)
	  // 	{
	  // 	  RealVector* TmpVecs1 = new RealVector[l - k];		    
	  // 	  RealVector* TmpVecs2 = new RealVector[l - k];		    
	  // 	  for (int m = k; m < l; ++m)
	  // 	    {
	  // 	      char* TmpEigenFile = new char[strlen(TmpEigenstateString) + 64];
	  // 	      sprintf (TmpEigenFile, "%s.%d.vec", TmpEigenstateString, m);
	  // 	      TmpVecs1[m - k].ReadVector(TmpEigenFile);
	  // 	      TmpVecs2[m - k] = RealVector(Chain->GetHilbertSpaceDimension());
	  // 	      TmpOperator.LowLevelMultiply(TmpVecs1[m - k], TmpVecs2[m - k], 0, Chain->GetHilbertSpaceDimension());
	  // 	      delete[] TmpEigenFile;
	  // 	    }
	  // 	  RealSymmetricMatrix TmpMatrix (l - k);
	  // 	  for (int m = k; m < l; ++m)
	  // 	    {
	  // 	      for (int n = m; n < l; ++n)
	  // 		{
	  // 		  double Tmp = (TmpVecs1[m - k] * TmpVecs2[n - k]);
	  // 		  TmpMatrix.SetMatrixElement (m - k, n - k, Tmp);
	  // 		}
	  // 	    }
	  // 	  RealDiagonalMatrix TmpMatrixDiag (l - k);
	  // 	  TmpMatrix.LapackDiagonalize(TmpMatrixDiag);
	  // 	  for (int m = k; m < l; ++m)
	  // 	    {
	  // 	      cout << m << ": " << TmpMatrixDiag[m - k] << " " << TmpSpectrum[m] << endl;
	  // 	    }
	  // 	  k = l;
	  // 	}
	  //     else
	  // 	{
	  // 	  char* TmpEigenFile = new char[strlen(TmpEigenstateString) + 64];
	  // 	  sprintf (TmpEigenFile, "%s.%d.vec", TmpEigenstateString, k);
	  // 	  RealVector TmpVec1 (Chain->GetHilbertSpaceDimension(), true);
	  // 	  RealVector TmpVec2 (Chain->GetHilbertSpaceDimension(), true);
	  // 	  TmpVec1.ReadVector(TmpEigenFile);	      
	  // 	  TmpOperator.LowLevelMultiply(TmpVec1, TmpVec2, 0, Chain->GetHilbertSpaceDimension());
	  // 	  cout << k << ": " << (TmpVec1 * TmpVec2) << " " << ((TmpVec2 * TmpVec2) - ((TmpVec1 * TmpVec2) * (TmpVec1 * TmpVec2))) << " " << (TmpVec1 * TmpVec1) << " " << TmpSpectrum[k] << endl;
	  // 	  ++k;
	  // 	  delete[] TmpEigenFile;
	  // 	}
	  //   }

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


// compute the characteristic polynomial for the real hamiltonians
//
// hamiltonian = pointer to the hamiltonian
// chain = pointer to the Hilbert space
// outputFileName = file name prefix for the characteristic polynomial
// architecture = pointer to the architecture
// sResolvedFlag = true if the characteristic polynomial should be computed per S sector
// minSValue = twice the minimum S value to consider (when sResolvedFlag is true)
// minSValue = twice the maximum S value to consider (when sResolvedFlag is true)
// sU2Degeneracy = array providing the Hilbert space dimension per S value (when sResolvedFlag is true)
// discardFourFactor = discard a global four factor used to ensure integer numbers

void SpinChainComputeCharacteristicPolynomial(SpinChainHamiltonian* hamiltonian, AbstractSpinChain* chain, char* outputFileName, AbstractArchitecture* architecture, bool sResolvedFlag, int minSValue, int maxSValue, int* sU2Degeneracy, bool discardFourFactor)
{
#ifdef __GMP__
  if (sResolvedFlag == false)
    {
      cout << "Computing the hamiltonian" << endl;
      RealMatrix TmpRawMatrix(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      hamiltonian->GetHamiltonian(TmpRawMatrix);
      double* TmpNormalizationFactors = chain->GetBasisNormalization();
      for (int i = 0; i < chain->GetHilbertSpaceDimension(); ++i)
	{
	  for (int j = 0; j < chain->GetHilbertSpaceDimension(); ++j)
	    {
	      double Tmp;
	      TmpRawMatrix.GetMatrixElement(i, j, Tmp);
	      Tmp *= TmpNormalizationFactors[i];
	      Tmp /= TmpNormalizationFactors[j];
	      TmpRawMatrix.SetMatrixElement(i, j, Tmp);
	    }
	}
      cout << "Converting to integer matrix" << endl;
      LongIntegerMatrix TmpMatrix;
      if (discardFourFactor == false)
	{
	  TmpMatrix = LongIntegerMatrix(TmpRawMatrix, 4.0);
	}
      else
	{
	  TmpMatrix = LongIntegerMatrix(TmpRawMatrix);
	}
      cout << "Start computing characteristic polynomial (degree " << chain->GetHilbertSpaceDimension() << ")" << endl;
      mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomial(architecture);
      char* PolynomialOutputFileName = new char[strlen(outputFileName) + 256];
      if (discardFourFactor == false)
	{
	  sprintf (PolynomialOutputFileName, "%s.charpol", outputFileName);
	}
      else
	{
	  sprintf (PolynomialOutputFileName, "%s.no4.charpol", outputFileName);
	}
      ofstream OutputFile;
      OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
      OutputFile << CharacteristicPolynomial[0];
      for (int i = 1; i <= chain->GetHilbertSpaceDimension(); ++i)
	{
	  OutputFile << "," << CharacteristicPolynomial[i];
	}
      OutputFile << endl;
      OutputFile.close();
    }
  else
    {
      cout << "Computing the hamiltonian and the S^2 operator" << endl;
      SpinS2Operator S2Operator(chain, chain->GetSpinChainLength());
      RealMatrix TmpRawMatrix(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      hamiltonian->GetHamiltonian(TmpRawMatrix);
      RealMatrix TmpRawMatrixS2(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      S2Operator.GetOperator(TmpRawMatrixS2);
      LongIntegerMatrix TmpMatrix(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      LongIntegerMatrix TmpMatrix2(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      LongIntegerMatrix TmpProjector(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      LongIntegerMatrix TmpMatrixS2(chain->GetHilbertSpaceDimension(), chain->GetHilbertSpaceDimension(), true);
      double* TmpNormalizationFactors = chain->GetBasisNormalization();
      for (int i = 0; i < chain->GetHilbertSpaceDimension(); ++i)
	{
	  for (int j = 0; j < chain->GetHilbertSpaceDimension(); ++j)
	    {
	      double Tmp;
	      TmpRawMatrix.GetMatrixElement(i, j, Tmp);
	      Tmp *= TmpNormalizationFactors[i];
	      Tmp /= TmpNormalizationFactors[j];
	      TmpRawMatrix.SetMatrixElement(i, j, Tmp);
	      TmpRawMatrixS2.GetMatrixElement(i, j, Tmp);
	      Tmp *= TmpNormalizationFactors[i];
	      Tmp /= TmpNormalizationFactors[j];
	      TmpRawMatrixS2.SetMatrixElement(i, j, Tmp);
	    }
	}

      for (int s = minSValue; s <= maxSValue; s += 2)
	{
	  if (sU2Degeneracy[(s - minSValue) >> 1] > 0)
	    {
	      cout << "computing 2S=" << s << " sector" << endl;
	      cout << "  building projected hamiltonian" << endl;
	      TmpProjector.SetToIdentity();
	      long TmpNormalisation = 1l;
	      for (int j = minSValue; j <= maxSValue; j += 2)
		{
		  if (j != s)
		    {
		      long TmpShift;
		      if ((minSValue & 1) == 0)
			{
			  TmpMatrixS2 = LongIntegerMatrix(TmpRawMatrixS2);
			  TmpShift = ((long) (-j * (j + 2))) >> 2;
			}
		      else
			{
			  TmpMatrixS2 = LongIntegerMatrix(TmpRawMatrixS2, 4.0);
			  TmpShift = ((long) (-j * (j + 2)));
			}
		      for (int k = 0; k < chain->GetHilbertSpaceDimension(); ++k)
			{
			  TmpMatrixS2.AddToMatrixElement(k, k, TmpShift);
			}
		      TmpProjector.Multiply(TmpMatrixS2);
		      if ((minSValue & 1) == 0)
			{
			  TmpNormalisation *= ((long) ((s * (s + 2)) - (j * (j + 2)))) >> 2;
			}
		      else
			{
			  TmpNormalisation *= ((long) ((s * (s + 2)) - (j * (j + 2))));
			}
		    }
		}
	      if (discardFourFactor == false)
		{
		  TmpMatrix = LongIntegerMatrix(TmpRawMatrix, 4.0);
		}
	      else
		{
		  TmpMatrix = LongIntegerMatrix(TmpRawMatrix);
		}
	      TmpMatrix2 = TmpMatrix * TmpProjector;
	      TmpMatrix = TmpProjector * TmpMatrix2;
	      TmpMatrix /= (TmpNormalisation);

	      cout << "Start computing characteristic polynomial (degree " << chain->GetHilbertSpaceDimension() << ")" << endl;
	      mpz_t TmpNormalisation2;
	      mpz_init_set_si(TmpNormalisation2, TmpNormalisation);
	      mpz_t* CharacteristicPolynomial = TmpMatrix.CharacteristicPolynomial(architecture);
	      for (int i = 1; i <= sU2Degeneracy[(s - minSValue) >> 1]; ++i)
		{		  
		  for (int j = 0; j < i; ++j)
		    {
		      mpz_divexact(CharacteristicPolynomial[chain->GetHilbertSpaceDimension() - i], CharacteristicPolynomial[chain->GetHilbertSpaceDimension() - i], TmpNormalisation2);
		    }
		}
	      cout << "  testing polynomial (checking if compatible with 2S=" << s << " with dim=" << sU2Degeneracy[(s - minSValue) >> 1] << ")" << endl;
	      int PolynomialTestingFlag = 0;
	      for (int i = 0; i < (chain->GetHilbertSpaceDimension() - sU2Degeneracy[(s - minSValue) >> 1]); ++i)
		{
		  if (mpz_sgn( CharacteristicPolynomial[i]) != 0)
		    {
		      PolynomialTestingFlag++;
		    }
		}
	      if (PolynomialTestingFlag == 0)
		{
		  cout << "  all clear" << endl;
		}
	      else
		{
		  cout << "  " << PolynomialTestingFlag << " error(s) detected" << endl;
		}
	      
	      char* PolynomialOutputFileName = new char[strlen(outputFileName) + 256];
	      if (discardFourFactor == false)
		{
		  sprintf (PolynomialOutputFileName, "%s_s_%d.charpol", outputFileName, s);
		}
	      else
		{
		  sprintf (PolynomialOutputFileName, "%s_s_%d.no4.charpol", outputFileName, s);
		}		
	      ofstream OutputFile;
	      OutputFile.open(PolynomialOutputFileName, ios::binary | ios::out);
	      OutputFile << CharacteristicPolynomial[chain->GetHilbertSpaceDimension() - sU2Degeneracy[(s - minSValue) >> 1]];
	      for (int i = chain->GetHilbertSpaceDimension() - sU2Degeneracy[(s - minSValue) >> 1] + 1; i <= chain->GetHilbertSpaceDimension(); ++i)
		{
		  OutputFile << "," << CharacteristicPolynomial[i];
		}
	      OutputFile << endl;
	      OutputFile.close();
	    }
	}
    }
#else
  cout << "GMP library is required for characteristic polynomials" << endl;
#endif	       
}
