#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnTorusWithSpin.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "GeneralTools/ListIterator.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/Options.h"

#include "MainTask/FQHEOnTorusMainTask.h"

#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);
    
  // some running options and help
  OptionManager Manager ("FQHETorusFermionsWithSpinSValue" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the file that contains the state whose average L value has to be evaluated");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusFermionsWithSpinSValue -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int TotalSz = 0;
  int NbrParticles = 0;
  int KyMax = 0;
  int TotalKy = 0;
  double XRatio = 0;
  bool Statistics;

  if (FQHEOnTorusWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("state"),
							  NbrParticles, KyMax, TotalKy, TotalSz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state")) == false)
    {
      cout << "state " << Manager.GetString("state") << " does not exist or can't be opened" << endl;
      return -1;           
    }
  RealVector State;
  if (State.ReadVector(Manager.GetString("state")) == false)
    {
      cout << "error while reading " << Manager.GetString("state") << endl;
      return -1;
    }


  FermionOnTorusWithSpin Space (NbrParticles, KyMax, TotalSz, TotalKy);	
   if (Space.GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space.GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
 
  Architecture.GetArchitecture()->SetDimension(Space.GetHilbertSpaceDimension());	
  AbstractQHEHamiltonian* Hamiltonian = new ParticleOnSphereWithSpinS2Hamiltonian(&Space, NbrParticles, KyMax - 1, TotalKy, TotalSz,
										  Architecture.GetArchitecture(), 1.0, 0);
  RealVector TmpState(Space.GetHilbertSpaceDimension());
  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
  Operation.ApplyOperation(Architecture.GetArchitecture());
  double S2Value = TmpState * State;
  double SValue = 0.5 * (sqrt ((4.0 * S2Value) + 1.0) - 1.0);
  cout << "<S^2> = " << S2Value << endl
       << "<S> = " << SValue << endl;

  delete Hamiltonian;
  return 0;
}
