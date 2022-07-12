#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"

#include "Hamiltonian/SpinChainAKLTStabilizerHamiltonian.h"

#include "HilbertSpace/Spin1_2ChainFull.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

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
  OptionManager Manager ("SpinChainAKLTAsStabilizer" , "0.01");
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
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  BooleanOption ('\n', "use-mirror", "use the mirror symmetry");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
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
      cout << "see man page for option syntax or type SpinChainAKLTAsStabilizer -h" << endl;
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
      sprintf (OutputFileName, "spin_%d_%schain_aklt_stabilizer_n_%d", (SpinValue / 2), BoundaryName, NbrSpins);
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
      sprintf (OutputFileName, "spin_%d_2_%schain_aklt_stabilizer_n_%d", SpinValue, BoundaryName, NbrSpins);
      if (Manager.GetBoolean("use-mirror") == true)
	{
	  sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz P_mirror", BoundaryName, SpinValue, NbrSpins);
	}
      else
	{
	  sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz", BoundaryName, SpinValue, NbrSpins);
	}
    }
    
  char* FullOutputFileName = new char [strlen(OutputFileName) + 64];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

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
      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
      sprintf (TmpEigenstateString, "%s", OutputFileName);
      char* TmpString = new char[1];
      TmpString[0] = '\0';
      Lanczos.SetRealAlgorithms();
      SpinChainAKLTStabilizerHamiltonian* Hamiltonian = new SpinChainAKLTStabilizerHamiltonian(Chain, NbrSpins, Manager.GetBoolean("use-periodic"));
      GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
			       FirstRun, TmpEigenstateString);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      delete[] TmpString;
      delete[] TmpEigenstateString;
    }

  // test code for the MPS of the stabilizer code 

//   ComplexMatrix* MPSMatrices = new ComplexMatrix [2];
//   MPSMatrices[0] = ComplexMatrix(2, 2);
//   MPSMatrices[1] = ComplexMatrix(2, 2);

//   MPSMatrices[0].SetMatrixElement(0 ,0, Phase(M_PI * 0.75) / M_SQRT2);
//   MPSMatrices[0].SetMatrixElement(1 ,0, 1.0 / M_SQRT2);
//   MPSMatrices[0].SetMatrixElement(0 ,1, 1.0 / M_SQRT2);
//   MPSMatrices[0].SetMatrixElement(1 ,1, Phase(-M_PI * 0.75) / M_SQRT2);

//   MPSMatrices[1].SetMatrixElement(0 ,0, Phase(-M_PI * 0.25) / M_SQRT2);
//   MPSMatrices[1].SetMatrixElement(1 ,0, 1.0 / M_SQRT2);
//   MPSMatrices[1].SetMatrixElement(0 ,1, 1.0 / M_SQRT2);
//   MPSMatrices[1].SetMatrixElement(1 ,1, Phase(M_PI * 0.25) / M_SQRT2);


//   cout << MPSMatrices[0] << endl;
//   cout << MPSMatrices[1] << endl;


//   ComplexVector TmpState (Chain->GetHilbertSpaceDimension(), true);
//   Chain->CreateStateFromMPSDescription(MPSMatrices, TmpState, -1, -1);
//   TmpState.Normalize();
//   TmpState.WriteVector("vector_pbc.vec");
//   Chain->CreateStateFromMPSDescription(MPSMatrices, TmpState, 0, 0);
//   TmpState.Normalize();
//   TmpState.WriteVector("vector_obc_0_0.vec");
//   Chain->CreateStateFromMPSDescription(MPSMatrices, TmpState, 0, 1);
//   TmpState.Normalize();
//   TmpState.WriteVector("vector_obc_0_1.vec");
//   Chain->CreateStateFromMPSDescription(MPSMatrices, TmpState, 1, 0);
//   TmpState.Normalize();
//   TmpState.WriteVector("vector_obc_1_0.vec");
//   Chain->CreateStateFromMPSDescription(MPSMatrices, TmpState, 1, 1);
//   TmpState.Normalize();
//   TmpState.WriteVector("vector_obc_1_1.vec");
//   for (int i = 0; i < Chain->GetHilbertSpaceDimension(); ++i)
//     Chain->PrintState (cout, i) << " : " << TmpState[i] << endl;

  delete Chain;
  delete[] OutputFileName;
  delete[] CommentLine;
  delete[] FullOutputFileName;


  

  return 0;
}
