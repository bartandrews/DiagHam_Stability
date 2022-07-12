#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnDisk.h"
#include "HilbertSpace/FermionOnDisk.h"
#include "HilbertSpace/FermionOnDiskUnlimited.h"
#include "HilbertSpace/ParticleOnSphere.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEForgeEigenstate" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (overriding the one found in the vector file name if greater than 0)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "momentum", "single particle momentum (overriding the one found in the vector file name if greater than 0)", 0, true, 0);
  (*SystemGroup) += new BooleanOption  ('\n', "bosons", "use boson statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "fermions", "use fermion statistics");
  (*SystemGroup) += new SingleStringOption  ('o', "output", "name of the vector output file", "output.vec");
  (*SystemGroup) += new SingleStringOption  ('i', "input", "name of the file that contains the eigenstate description");
  (*SystemGroup) += new SingleStringOption  ('\n', "additional-state", "name of the file that contains an optional state that will be added to the one describde by the input file");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "additional-coefficient", "multiplicative coefficient that will be used on the additional state", 0.0);
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEForgeEigenstate -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (((SingleStringOption*) Manager["input"])->GetString() == 0)
    {
      cout << "QHEBosonsCorrelation requires a state" << endl;
      return -1;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["momentum"])->GetInteger();
  bool BosonFlag = ((BooleanOption*) Manager["bosons"])->GetBoolean();
  bool FermionFlag = ((BooleanOption*) Manager["fermions"])->GetBoolean();
  char* OutputName = ((SingleStringOption*) Manager["output"])->GetString();
  char* InputName = ((SingleStringOption*) Manager["input"])->GetString();
  
  RealVector State;
  if (((SingleStringOption*) Manager["additional-state"])->GetString() != 0)
    {
      State.ReadVector(((SingleStringOption*) Manager["additional-state"])->GetString());
      State *= ((SingleDoubleOption*) Manager["additional-coefficient"])->GetDouble();
    }
  ParticleOnSphere* Space = 0;
  if (BosonFlag == true)
    {
      Space = new BosonOnDisk(NbrParticles, Lz);
    }
  else
    {
#ifdef __64_BITS__
      if ((Lz - (((NbrParticles - 1) * (NbrParticles - 2)) / 2)) < 63)      
#else
      if ((Lz - (((NbrParticles - 1) * (NbrParticles - 2)) / 2)) < 31)
#endif
	Space = new FermionOnDisk (NbrParticles, Lz);
      else
	Space = new FermionOnDiskUnlimited (NbrParticles, Lz);      
    }
  ((BosonOnDisk*) Space)->ForgeEigenstate(InputName, State);
  if (State.WriteVector(OutputName) == false)
    {
      cout << "can't open vector file " << OutputName << endl;
      return -1;      
    }

  delete Space;
}
