#include "Vector/RealVector.h"

#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"

#include "Options/Options.h"
#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;

int main(int argc, char** argv)
{
  cout.precision(14);
  
  OptionManager Manager ("FQHESphereBosonsFuseParticles" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  //ArchitectureManager Architecture;
  
  Manager += SystemGroup;
  //Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  Manager += OutputGroup;
  
  (*SystemGroup) += new SingleStringOption ('\0', "state", "name of the vector files which basis will be fused");
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the vector decomposition into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the vector decomposition into a text file");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsFuseParticles -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if(Manager.GetString("state") == 0)
    {
      cout << "no input state" << endl << "see man page for option syntax or type FQHESphereLValue -h" << endl;
      return -1;
    }
  
  
  int NbrParticles = 0;
  int TotalLz = 0;
  int LzMax = 0;
  char* OutputFileName = Manager.GetString("bin-output");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  bool FermionFlag = false;
  
  
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMax, TotalLz, FermionFlag) == false)
    {
      return -1;
    }
  
  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }
  
  char* StateFileName = Manager.GetString("state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;
    }
  
  RealVector GroundState;
  if (GroundState.ReadVector (StateFileName) == false)
    {
      cout << "can't open vector file " << StateFileName << endl;
      return -1;      
    }
  
  if((NbrParticles & 1) == 1)
    {
      cout <<"the number of particles must be even"<<endl;
      return -1;
    }
  
  BosonOnSphereShort * OutputBasis=0;
  OutputBasis = new BosonOnSphereShort(NbrParticles, TotalLz, LzMax);
  
  BosonOnSphereShort * FinalSpace=0;
  
  FinalSpace = new BosonOnSphereShort(NbrParticles / 2, TotalLz, 2 * LzMax);
  
  
  if (OutputBasis->GetHilbertSpaceDimension() != GroundState.GetVectorDimension())
    {
      cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
      return -1;
    }
  
  RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
  
  OutputBasis->FuseParticlesInState(GroundState, OutputVector, FinalSpace, 0, OutputBasis->GetHilbertSpaceDimension());
  OutputVector.WriteVector(OutputFileName);
  if(OutputTxtFileName != 0)
    {
      ofstream File;
      File.open(OutputTxtFileName, ios::binary | ios::out);
      File.precision(14);
      for (long i = 0; i < FinalSpace->GetLargeHilbertSpaceDimension(); ++i)
	{
	  File << OutputVector[i] << " ";
	  FinalSpace->PrintStateMonomial(File, i) << endl;
	}
      File.close();
    }
  return 0;
}
