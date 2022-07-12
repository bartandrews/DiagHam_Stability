#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHEDiskQuasiholePropagatorOperation.h"

#include "Options/Options.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHEDiskQuasiholePropagator" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;


  (*SystemGroup) += new SingleStringOption  ('\n', "base-state", "vector file that describes the base state ");
  (*SystemGroup) += new SingleStringOption  ('\n', "base-reference", "use a file as the definition of the reference state of the base state");
  (*SystemGroup) += new BooleanOption  ('\n', "basehuge-basis", "use huge Hilbert space support for the base state");
  (*SystemGroup) += new SingleStringOption  ('\n', "excited-state", "vector file that describes the excited state");
  (*SystemGroup) += new SingleStringOption  ('\n', "excited-reference", "use a file as the definition of the reference state of the excited state");
  (*SystemGroup) += new BooleanOption  ('\n', "excitedhuge-basis", "use huge Hilbert space support for the excited state");
  (*SystemGroup) += new BooleanOption  ('\n', "no-base", "do not compute the contribution fron the base state");  
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "consider fermions instead of bosons");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "baseload-hilbert", "load Hilbert space description from the base state",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "excitedload-hilbert", "load Hilbert space description from the excited file",0);  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskQuasiholePropagator -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  bool Statistics = Manager.GetBoolean("fermion");
  int ExcitedNbrParticles = 0;
  int ExcitedTotalLz = 0;
  int ExcitedLzMax = 0;
  int BaseNbrParticles = 0;
  int BaseTotalLz = 0;
  int BaseLzMax = 0;
  int DiffLzMax = 0;
  double Base = 0.0;

  if (Manager.GetBoolean("no-base") == false)
    {

      RealVector BaseVector;
      if (BaseVector.ReadVector (Manager.GetString("base-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("base-state") << endl;
	  return -1;      
	}
      
      ParticleOnSphere* BaseBasis = 0;
      if (Statistics == false)
	{
	  if (Manager.GetBoolean("basehuge-basis") == true)
	    {
	      if (Manager.GetString("baseload-hilbert") == 0)
		{
		  cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		  return -1;
		}
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("base-reference"), BaseNbrParticles, BaseLzMax, ReferenceState) == false)
		return -1;
	      for (int i = 0; i <= BaseLzMax; ++i)
		DiffLzMax -= i * ReferenceState[i];
	      BaseBasis = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("baseload-hilbert"), Manager.GetInteger("memory"));
	    }
	  else
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("base-reference"), BaseNbrParticles, BaseLzMax, ReferenceState) == false)
		return -1;
	      for (int i = 0; i <= BaseLzMax; ++i)
		DiffLzMax -= i * ReferenceState[i];
	      if (Manager.GetString("baseload-hilbert") != 0)
		BaseBasis = new BosonOnSphereHaldaneBasisShort(Manager.GetString("baseload-hilbert"));	  
	      else
		BaseBasis = new BosonOnSphereHaldaneBasisShort(BaseNbrParticles, BaseTotalLz, BaseLzMax, ReferenceState);	  
	    }
	}
      else
	{
	  if (Manager.GetBoolean("basehuge-basis") == true)
	    {
	      if (Manager.GetString("baseload-hilbert") == 0)
		{
		  cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		  return -1;
		}
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("base-reference"), BaseNbrParticles, BaseLzMax, ReferenceState) == false)
		return -1;
	      for (int i = 0; i <= BaseLzMax; ++i)
		DiffLzMax -= i * ReferenceState[i];
	      BaseBasis = new FermionOnSphereHaldaneHugeBasis (Manager.GetString("baseload-hilbert"), Manager.GetInteger("memory"));
	    }
	  else
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("base-reference"), BaseNbrParticles, BaseLzMax, ReferenceState) == false)
		return -1;
	      for (int i = 0; i <= BaseLzMax; ++i)
		DiffLzMax -= i * ReferenceState[i];
	      if (Manager.GetString("baseload-hilbert") != 0)
		BaseBasis = new FermionOnSphereHaldaneBasis(Manager.GetString("baseload-hilbert"));	  
	      else
		BaseBasis = new FermionOnSphereHaldaneBasis(BaseNbrParticles, BaseTotalLz, BaseLzMax, ReferenceState);
	    }
	}
    
      FQHEDiskQuasiholePropagatorOperation Operation(BaseBasis, &BaseVector);
      Operation.ApplyOperation(Architecture.GetArchitecture());
      Base = Operation.GetScalar().Re;
      delete BaseBasis;
      BaseVector = RealVector();
    }

  RealVector ExcitedVector;
  if (ExcitedVector.ReadVector (Manager.GetString("excited-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("excited-state") << endl;
      return -1;      
    }

  ParticleOnSphere* ExcitedBasis = 0;
  if (Statistics == false)
    {
      if (Manager.GetBoolean("excitedhuge-basis") == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("excited-reference"), ExcitedNbrParticles, ExcitedLzMax, ReferenceState) == false)
	    return -1;
	  for (int i = 0; i <= ExcitedLzMax; ++i)
	    DiffLzMax += i * ReferenceState[i];
	  if (Manager.GetString("excitedload-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  ExcitedBasis = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("excitedload-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("excited-reference"), ExcitedNbrParticles, ExcitedLzMax, ReferenceState) == false)
	    return -1;
	  for (int i = 0; i <= ExcitedLzMax; ++i)
	    DiffLzMax += i * ReferenceState[i];
	  if (Manager.GetString("excitedload-hilbert") != 0)
	    ExcitedBasis = new BosonOnSphereHaldaneBasisShort(Manager.GetString("excitedload-hilbert"));	  
	  else
	    ExcitedBasis = new BosonOnSphereHaldaneBasisShort(ExcitedNbrParticles, ExcitedTotalLz, ExcitedLzMax, ReferenceState);	  
	}
    }
  else
    {
      if (Manager.GetBoolean("excitedhuge-basis") == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("excited-reference"), ExcitedNbrParticles, ExcitedLzMax, ReferenceState) == false)
	    return -1;
	  for (int i = 0; i <= ExcitedLzMax; ++i)
	    DiffLzMax += i * ReferenceState[i];
	  if (Manager.GetString("excitedload-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  ExcitedBasis = new FermionOnSphereHaldaneHugeBasis (Manager.GetString("excitedload-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("excited-reference"), ExcitedNbrParticles, ExcitedLzMax, ReferenceState) == false)
	    return -1;
	  for (int i = 0; i <= ExcitedLzMax; ++i)
	    DiffLzMax += i * ReferenceState[i];
	  if (Manager.GetString("excitedload-hilbert") != 0)
	    ExcitedBasis = new FermionOnSphereHaldaneBasis(Manager.GetString("excitedload-hilbert"));	  
	  else
	    ExcitedBasis = new FermionOnSphereHaldaneBasis(ExcitedNbrParticles, ExcitedTotalLz, ExcitedLzMax, ReferenceState);
	}
    }

  
  FQHEDiskQuasiholePropagatorOperation Operation2(ExcitedBasis, &ExcitedVector);
  Operation2.ApplyOperation(Architecture.GetArchitecture());
  double Excited = Operation2.GetScalar().Re;
  if (Manager.GetBoolean("no-base") == false)
    {
      cout << "2^" << DiffLzMax <<  " * " << Excited << " / " << Base << " = " << (pow(2.0, (double) DiffLzMax) * Excited / Base) << endl;
    }
  else
    {
      cout << Excited << endl;
    }

  delete ExcitedBasis;
  return 0;
}
