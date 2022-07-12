#include "Operator/PottsNonLocalOccupationNumberOperator.h"

#include "HilbertSpace/Spin1_2ChainFixedParity.h"
#include "HilbertSpace/Potts3Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/SpinFiles/SpinFileTools.h"

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
  OptionManager Manager ("PottsChainNonLocalOccupationNumber" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('n', "zn", "Zn value of the Potts model", 3);
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PottsChainNonLocalOccupationNumber -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int ZnValue = Manager.GetInteger("zn");
  int NbrSpins = 0;
  int SzValue = 0;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  bool SVDFlag = Manager.GetBoolean("use-svd");
  bool RealVectorFlag =  Manager.GetBoolean("real-vectors");

  ComplexVector* GroundStates = 0;
  RealVector* RealGroundStates = 0;
  double* Weights =0;
  bool WeightFlag = false;
  int NbrSpaces = 1;
  char** GroundStateFiles = 0;
  AbstractSpinChain** Spaces = 0;
  int* TotalSz = 0;
  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalSz = new int[1];
      Weights = new double[1];
      Weights[0] = 1.0;
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));      
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
      NbrSpaces = DegeneratedFile.GetNbrLines();
      GroundStateFiles = new char* [NbrSpaces];
      TotalSz = new int[NbrSpaces];
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	  strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	}
      if (DegeneratedFile.GetNbrColumns() > 1)
	{
	  Weights = DegeneratedFile.GetAsDoubleArray(1);
	  WeightFlag = true;
	}
      else
	{
	  Weights = new double[NbrSpaces];
	  for (int i = 0; i < NbrSpaces; ++i)
	    Weights[i] = 1.0;
	}
    }
  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalSz[i] = 0;
      if (PottsFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrSpins, TotalSz[i]) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
    }
  
  if (RealVectorFlag == false)
    {
      GroundStates = new ComplexVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }
	}
    }
  else
    {
      RealGroundStates = new RealVector [NbrSpaces];  
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (RealGroundStates[i].ReadVector(GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }
	}
    }
  
  
  ZnValue = Manager.GetInteger("zn");
  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      cout << "Filename: " << GroundStateFiles[i] << " N= " << NbrSpins << " Q= " << TotalSz[i] << endl;
    }
  
  
  Spaces = new AbstractSpinChain* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      switch (ZnValue)
	{
 	case 2 :
 	  Spaces[i] = new Spin1_2ChainFixedParity (NbrSpins, TotalSz[i]);
 	  break;
	case 3 :
	  Spaces[i] = new Potts3Chain (NbrSpins, TotalSz[i], 1000000);
	  break;
	default :
	  {
	    cout << "Zn with n=" << ZnValue  << " is not available" << endl;
	    return -1;
	  }
       }
    }


  for (int i = 0; i < NbrSpaces; ++i)
    {
      switch (ZnValue)
	{
	case 2 :
	  {
	    PottsNonLocalOccupationNumberOperator Operator((Spin1_2Chain*) Spaces[i], NbrSpins);
	    Complex Tmp1 = Operator.MatrixElement(RealGroundStates[i], RealGroundStates[i]);
	  } 
	  break; 
	case 3 :
	  {
	    //	    FermionParityOperator Operator(Chain, NbrSpinsPeriodic);
	  } 
	  break;
	} 
    }
  return 0;
}
