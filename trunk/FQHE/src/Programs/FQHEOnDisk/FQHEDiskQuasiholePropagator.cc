#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneLargeBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHEDiskQuasiholePropagatorOperation.h"

#include "Options/Options.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;


// get the Hilbert  space
//
// statistics = true if fermions have to be used
// hugeFlag = true if huge basis has to be used
// referenceFileName = reference file name for the squeezed basis
// loadHilbertFileName = name of the file that contains the Hilbert space description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum angular momentum
// totalLz = reference on twice the total Lz projection
// totalSz = reference on twice the total Sz projection
// diffLzMax = reference on the difference between LzMax, TotalLz will be added to the current value
// su2Flag = true for spinful state
// memory = amount of memory allowed for huge mode
// return value = pointer to the Hilbert space
ParticleOnSphere* FQHEDiskQuasiholePropagatorGetHilbertSpace(bool statistics, bool hugeFlag, char* referenceFileName, char* loadHilbertFileName, int& nbrParticles, int& lzMax, int& totalLz, int& totalSz, int& diffLzMax, bool su2Flag, long memory);


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
  (*SystemGroup) += new SingleStringOption  ('\n', "excited-2ndstate", "vector file that describes the second excited state if a scalar product has to be computed (state should be in the same squeezed basis than excited-state)");
  (*SystemGroup) += new BooleanOption  ('\n', "no-base", "do not compute the contribution fron the base state");  
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "consider fermions instead of bosons");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*SystemGroup) += new BooleanOption  ('\n', "rational" , "use rational numbers instead of double precision floating point numbers");
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
  int TotalSz = 0;
  double Base = 0.0;
  LongRational RationalBase = 0l;
  bool SU2Flag = false;
  if (strstr(Manager.GetString("excited-state"), "su2") != 0)
    SU2Flag = true;

  if (Manager.GetBoolean("no-base") == false)
    {      
      ParticleOnSphere* BaseBasis = FQHEDiskQuasiholePropagatorGetHilbertSpace(Statistics, Manager.GetBoolean("basehuge-basis"), Manager.GetString("base-reference"), Manager.GetString("baseload-hilbert"), 
									       BaseNbrParticles, BaseLzMax, BaseTotalLz, TotalSz, DiffLzMax, SU2Flag, Manager.GetInteger("memory"));
      if (BaseBasis == 0)
	return -1;
      DiffLzMax *= -1;	

      if (Manager.GetBoolean("rational") == true)
	{
	  LongRationalVector BaseVector;
	  if (BaseVector.ReadVector (Manager.GetString("base-state")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("base-state") << endl;
	      return -1;      
	    }
	  FQHEDiskQuasiholePropagatorOperation Operation(BaseBasis, &BaseVector);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  RationalBase = Operation.GetLongRationalScalar();
	  BaseVector = LongRationalVector();
	}
      else
	{
	  RealVector BaseVector;
	  if (BaseVector.ReadVector (Manager.GetString("base-state")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("base-state") << endl;
	      return -1;      
	    }
	  FQHEDiskQuasiholePropagatorOperation Operation(BaseBasis, &BaseVector);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  Base = Operation.GetScalar().Re;
	  BaseVector = RealVector();
	}
      delete BaseBasis;
    }

  ParticleOnSphere* ExcitedBasis = FQHEDiskQuasiholePropagatorGetHilbertSpace(Statistics, Manager.GetBoolean("excitedhuge-basis"), Manager.GetString("excited-reference"), Manager.GetString("excitedload-hilbert"), 
									      ExcitedNbrParticles, ExcitedLzMax, ExcitedTotalLz, TotalSz, DiffLzMax, SU2Flag, Manager.GetInteger("memory"));
  if (ExcitedBasis == 0)
    return -1;
  if (Manager.GetBoolean("rational") == true)
    {
      LongRational Excited = 0l;
      LongRationalVector ExcitedVector;
      if (ExcitedVector.ReadVector (Manager.GetString("excited-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("excited-state") << endl;
	  return -1;      
	}
      if (Manager.GetString("excited-2ndstate") == 0)
	{
	  FQHEDiskQuasiholePropagatorOperation Operation2(ExcitedBasis, &ExcitedVector);
	  Operation2.ApplyOperation(Architecture.GetArchitecture());
	  Excited = Operation2.GetLongRationalScalar();
	}
      else
	{
	  LongRationalVector ExcitedVector2;
	  if (ExcitedVector2.ReadVector (Manager.GetString("excited-2ndstate")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("excited-2ndstate") << endl;
	      return -1;      
	    }
	  FQHEDiskQuasiholePropagatorOperation Operation2(ExcitedBasis, &ExcitedVector, &ExcitedVector2);
	  Operation2.ApplyOperation(Architecture.GetArchitecture());
	  Excited = Operation2.GetLongRationalScalar();
	}
      if (Manager.GetBoolean("no-base") == false)
	{
	  LongRational Tmp = Excited / RationalBase;
	  Tmp.Power2Multiply(DiffLzMax);
	  LongRational Tmp2 = Tmp;
	  Tmp2.Power2Multiply(DiffLzMax);
	  cout << "2^" << DiffLzMax <<  " * " << Excited << " / " << RationalBase << " = " << Tmp2 << " = " << Tmp2.GetNumericalValue() << endl;
	}
      else
	{
	  cout << Excited << " = " << Excited.GetNumericalValue() << endl;
	}
    }
  else
    {
      double Excited = 0.0;
      RealVector ExcitedVector;
      if (ExcitedVector.ReadVector (Manager.GetString("excited-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("excited-state") << endl;
	  return -1;      
	}
      if (Manager.GetString("excited-2ndstate") == 0)
	{
	  FQHEDiskQuasiholePropagatorOperation Operation2(ExcitedBasis, &ExcitedVector);
	  Operation2.ApplyOperation(Architecture.GetArchitecture());
	  Excited = Operation2.GetScalar().Re;
	}
      else
	{
	  RealVector ExcitedVector2;
	  if (ExcitedVector2.ReadVector (Manager.GetString("excited-2ndstate")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("excited-2ndstate") << endl;
	      return -1;      
	    }
	  FQHEDiskQuasiholePropagatorOperation Operation2(ExcitedBasis, &ExcitedVector, &ExcitedVector2);
	  Operation2.ApplyOperation(Architecture.GetArchitecture());
	  Excited = Operation2.GetScalar().Re;
	}
      if (Manager.GetBoolean("no-base") == false)
	{
	  cout << "2^" << DiffLzMax <<  " * " << Excited << " / " << Base << " = " << (pow(2.0, (double) DiffLzMax) * Excited / Base) << endl;
	}
      else
	{
	  cout << Excited << endl;
	}
    }

  delete ExcitedBasis;
  return 0;
}


// get the Hilbert  space
//
// statistics = true if fermions have to be used
// hugeFlag = true if huge basis has to be used
// referenceFileName = reference file name for the squeezed basis
// loadHilbertFileName = name of the file that contains the Hilbert space description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum angular momentum
// totalLz = reference on twice the total Lz projection
// totalSz = reference on twice the total Sz projection
// diffLzMax = reference on the difference between LzMax, TotalLz will be added to the current value
// su2Flag = true for spinful state
// memory = amount of memory allowed for huge mode
// return value = pointer to the Hilbert space

ParticleOnSphere* FQHEDiskQuasiholePropagatorGetHilbertSpace(bool statistics, bool hugeFlag, char* referenceFileName, char* loadHilbertFileName, int& nbrParticles, int& lzMax, int& totalLz, int& totalSz, int& diffLzMax, bool su2Flag, long memory)
{
  ParticleOnSphere* Basis = 0;
  if (su2Flag == false)
    {
      if (statistics == false)
	{
	  if (hugeFlag == true)
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(referenceFileName, nbrParticles, lzMax, ReferenceState) == false)
		return 0;
	      for (int i = 0; i <= lzMax; ++i)
		diffLzMax += i * ReferenceState[i];
	      if (loadHilbertFileName == 0)
		{
		  cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		  return 0;
		}
	      Basis = new BosonOnSphereHaldaneHugeBasisShort (loadHilbertFileName, memory);
	    }
	  else
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(referenceFileName, nbrParticles, lzMax, ReferenceState) == false)
		return 0;
	      for (int i = 0; i <= lzMax; ++i)
		diffLzMax += i * ReferenceState[i];
	      if (loadHilbertFileName != 0)
		Basis = new BosonOnSphereHaldaneBasisShort(loadHilbertFileName);	  
	      else
		Basis = new BosonOnSphereHaldaneBasisShort(nbrParticles, totalLz, lzMax, ReferenceState);	  
	    }
	}
      else
	{
	  if (hugeFlag == true)
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(referenceFileName, nbrParticles, lzMax, ReferenceState) == false)
		return 0;
	      for (int i = 0; i <= lzMax; ++i)
		diffLzMax += i * ReferenceState[i];
	      if (loadHilbertFileName == 0)
		{
		  cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		  return 0;
		}
	      Basis = new FermionOnSphereHaldaneHugeBasis (loadHilbertFileName, memory);
	    }
	  else
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(referenceFileName, nbrParticles, lzMax, ReferenceState) == false)
		return 0;
	      for (int i = 0; i <= lzMax; ++i)
		diffLzMax += i * ReferenceState[i];
	      if (loadHilbertFileName != 0)
		Basis = new FermionOnSphereHaldaneBasis(loadHilbertFileName);	  
	      else
		Basis = new FermionOnSphereHaldaneBasis(nbrParticles, totalLz, lzMax, ReferenceState);
	    }
	}
    }
  else
    {
      int** ReferenceStates = 0;
      int NbrReferenceStates;
      int totalSz = 0;
      if (referenceFileName == 0)
	{
	  cout << "error, a reference file is needed for fermions in Haldane basis" << endl;
	  return 0;
	}
      if (FQHEGetRootPartitionSU2(referenceFileName, nbrParticles, lzMax, ReferenceStates, NbrReferenceStates) == false)
	{
	  cout << "error while parsing " << referenceFileName << endl;
	  return 0;
	}
      if (hugeFlag == true)
        {
          if (loadHilbertFileName != 0)
            {
              Basis = new FermionOnSphereWithSpinHaldaneLargeBasis(loadHilbertFileName);
            }
          else
            {
              Basis = new FermionOnSphereWithSpinHaldaneLargeBasis(nbrParticles, totalLz, lzMax, totalSz, ReferenceStates, NbrReferenceStates); 
            }
        }
      else
        {	
 #ifdef __64_BITS__
          if (lzMax <= 31)
#else
	  if (lzMax <= 15)
#endif
	    {
	      Basis = new FermionOnSphereWithSpinHaldaneBasis(nbrParticles, totalLz, lzMax, totalSz, ReferenceStates, NbrReferenceStates);
	    }
	  else
	    {
#ifdef __128_BIT_LONGLONG__
	      if (lzMax <= 63)
#else
	        if (lzMax <= 31)
#endif
		  {
		    Basis = new FermionOnSphereWithSpinHaldaneBasisLong (nbrParticles, totalLz, lzMax, totalSz, ReferenceStates, NbrReferenceStates);
		  }
	        else
	   	  {
		    cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		    return 0;
		  }
	    }
      }
    }
  return Basis;
}
