#include "config.h"

#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereDroplet.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereConvertHaldaneBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "input state file name");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "number of flux quanta", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system", -1);
 
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-fluxes1", "number of fluxes in a droplet", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-particles1", "max number of particles in a droplet", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-holes1", "max number of holes in a droplet", 0);

  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-fluxes2", "secondary condition for number of fluxes in a droplet", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-particles2", "secondary condition for  max number of particles in a droplet", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-holes2", "secondary condition for max number of holes in a droplet", 0);

  (*SystemGroup) += new BooleanOption  ('f', "fermion", "use fermionic statistics", false);
  (*SystemGroup) += new BooleanOption  ('b', "boson", "use bosonic statistics", false);
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support (only available when both the source and target spaces are squeezed basis)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (removing any occurence of haldane_)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereConvertDropletToNBodyBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-file") == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHESphereConvertDropletToNBodyBasis -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("input-file")) == false)
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int LzMax = Manager.GetInteger("lzmax"); 
  int TotalLz  = Manager.GetInteger("total-lz");

  int NbrFluxes1 = Manager.GetInteger("nbr-fluxes1");
  int MaxNbrParticles1 = Manager.GetInteger("max-particles1");
  int MaxNbrHoles1 = Manager.GetInteger("max-holes1");

  int NbrFluxes2 = Manager.GetInteger("nbr-fluxes2");
  int MaxNbrParticles2 = Manager.GetInteger("max-particles2");
  int MaxNbrHoles2 = Manager.GetInteger("max-holes2");

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  if (Manager.GetBoolean("boson") == true)
   {
     cout << "Bosons currently not supported."<<endl;
     exit(1);
   }

  if (((NbrParticles * LzMax) & 1) != (TotalLz & 1))
    {
      cout << "incompatible values for nbr-particles, nbr-flux and total-lz" << endl;
      return -1;
    }

  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      cout << "Error: output file needed." << endl; 
    }


  FermionOnSphereDroplet* InputSpace = new FermionOnSphereDroplet(NbrParticles, TotalLz, LzMax, NbrFluxes1, MaxNbrParticles1, MaxNbrHoles1, NbrFluxes2, MaxNbrParticles2, MaxNbrHoles2, Memory);

  FermionOnSphere* OutputSpace = new FermionOnSphere(NbrParticles, TotalLz, LzMax, Memory);

  RealVector State;
  if (State.ReadVector (Manager.GetString("input-file")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-file") << endl;
      return -1;      
    }
  if (InputSpace->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension())
   { 
     cout << "dimension mismatch between Hilbert space and input state" << endl;
     return -1;      
   }

  RealVector OutputState;
  OutputState = InputSpace->ConvertToNbodyBasis(State, *OutputSpace);
  if (OutputState.WriteVector(OutputFileName) == false)
   {
     cout << "error while writing output state " << OutputFileName << endl;
     return -1;
   }
}
