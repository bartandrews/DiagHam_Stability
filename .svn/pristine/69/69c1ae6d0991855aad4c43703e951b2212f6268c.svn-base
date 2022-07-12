#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"

#include "MathTools/ClebschGordanCoefficients.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereASCIIToState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 20);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (only useful in su(2)/su(4) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "su3-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-tz", "twice the quantum number of the system associated to the Tz generator (only useful in su(3) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-y", "three time the quantum number of the system associated to the Y generator (only useful in su(3) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "su4-spin", "consider particles with SU(4) spin");
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin of the system (only usefull in su(4) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "2-ll", "consider particles within two Landau levels");
  (*SystemGroup) += new SingleStringOption ('\0', "ascii-state", "name of the input ASCII description of the state (should use the same convention than FQHESphereShowBasis output)");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of trying to replace .txt extension with .vec or appending .vec extension");
  (*SystemGroup) += new BooleanOption  ('\n', "no-normalization", "do not normalize the final state");
  (*SystemGroup) += new BooleanOption  ('\n', "rational", "use rational numbers instead of double precision floating point numbers");
  (*SystemGroup) += new BooleanOption  ('c', "complex", "create a vector with complex components");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereASCIIToState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if (Manager.GetString("ascii-state") == 0)
    {
      cout << "An input file has to be provided" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux"); 
  int TotalLz = Manager.GetInteger("lz-value");
  int TotalTz = Manager.GetInteger("total-tz");
  int TotalY = Manager.GetInteger("total-y");
  bool SU2SpinFlag = Manager.GetBoolean("su2-spin");
  bool SU3SpinFlag = Manager.GetBoolean("su3-spin");
  bool SU4SpinFlag = Manager.GetBoolean("su4-spin");
  bool TwoLLFlag = Manager.GetBoolean("2-ll");
  int TotalSz = Manager.GetInteger("total-sz");
    
  if (((NbrParticles * NbrFluxQuanta) & 1) != (TotalLz & 1)) 
    {
      cout << "incompatible values for the number of particles, the number of flux quanta and twice the total lz value (nbr-particles * nbr-flux and lz-value should have the same parity)" << endl;
      return -1;
    }

  MultiColumnASCIIFile InputFile(':');
  if (InputFile.Parse(Manager.GetString("ascii-state")) == false)
    {
      InputFile.DumpErrors(cout) << endl;
      return -1;
    }
  double* Coefficients = InputFile.GetAsDoubleArray(1);
  if (Coefficients == 0)
    {
      InputFile.DumpErrors(cout) << endl;
      return -1;
    }

  ParticleOnSphere* Space = 0;
  if (Manager.GetBoolean("boson") == true)
    {
      if (Manager.GetBoolean("haldane") == false)
	{
	  Space = new BosonOnSphereShort(NbrParticles, TotalLz, NbrFluxQuanta);
	}
      else
	{
	  int* ReferenceState = 0;
	  if (Manager.GetString("reference-file") == 0)
	    {
	      cout << "error, a reference file is needed" << endl;
	      return 0;
	    }
	  ConfigurationParser ReferenceStateDefinition;
	  if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
	    {
	      ReferenceStateDefinition.DumpErrors(cout) << endl;
	      return 0;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	    {
	      cout << "NbrParticles is not defined or as a wrong value" << endl;
	      return 0;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", NbrFluxQuanta) == false) || (NbrFluxQuanta < 0))
	    {
	      cout << "LzMax is not defined or as a wrong value" << endl;
	      return 0;
	    }
	  int MaxNbrLz;
	  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	    {
	      cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
	      return 0;     
	    }
	  if (MaxNbrLz != (NbrFluxQuanta + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return 0;     
	    }
	   Space = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);	  
	}
    }
  else
    {
      if ((SU2SpinFlag == false) && (SU3SpinFlag == false) && (SU4SpinFlag == false) && (TwoLLFlag == false))
	{
	  if (Manager.GetBoolean("haldane") == false)
	    {
#ifdef __64_BITS__
	      if (NbrFluxQuanta <= 63)
		Space = new FermionOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, NbrFluxQuanta);
#else
	      if (NbrFluxQuanta <= 31)
		Space = new FermionOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, NbrFluxQuanta);
#endif
	    }
	  else
	    {
	      int* ReferenceState = 0;
	      if (Manager.GetString("reference-file") == 0)
		{
		  cout << "error, a reference file is needed" << endl;
	      return 0;
		}
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
		{
		  ReferenceStateDefinition.DumpErrors(cout) << endl;
		  return 0;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
		{
		  cout << "NbrParticles is not defined or as a wrong value" << endl;
		  return 0;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", NbrFluxQuanta) == false) || (NbrFluxQuanta < 0))
		{
		  cout << "LzMax is not defined or as a wrong value" << endl;
		  return 0;
		}
	      int MaxNbrLz;
	      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		{
		  cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
		  return 0;     
		}
	      if (MaxNbrLz != (NbrFluxQuanta + 1))
		{
		  cout << "wrong LzMax value in ReferenceState" << endl;
		  return 0;     
		}
	      Space = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);	  
	    }
	}
      else
	{
	  if (SU2SpinFlag == true)
	    Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
	  else
	    if (SU3SpinFlag == true)
	      Space = new FermionOnSphereWithSU3Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalTz, TotalY);	  
	    else
	      if (TwoLLFlag == true)
		Space = new FermionOnSphereTwoLandauLevels(NbrParticles, TotalLz, NbrFluxQuanta + 2, NbrFluxQuanta);	  
	}
    }
  char* OutputFile = Manager.GetString("output-file");
  if (OutputFile == 0)
    {
      OutputFile = ReplaceExtensionToFileName(Manager.GetString("ascii-state"), "txt", "vec");
      if (OutputFile == 0)
	OutputFile = AddExtensionToFileName(Manager.GetString("ascii-state"), "vec");
    }
  if ((Manager.GetBoolean("rational") == false) && (Manager.GetBoolean("complex") == false))
  {
    RealVector State (Space->GetHilbertSpaceDimension(), true);

    for (int i = 0; i < InputFile.GetNbrLines(); ++i)
      {
	char* TmpString =  InputFile(0, i);
	CleanLine(TmpString);
	int TmpIndex = Space->FindStateIndex(TmpString);
	if (TmpIndex != -1)
	  {
	    State[TmpIndex] = Coefficients[i];
	  }
	else
	  {
	    cout << "warning , invalid state |" <<  TmpString << ">" << endl;
	  }
      }
    if (Manager.GetBoolean("no-normalization") == false)
      State /= State.Norm();
  State.WriteVector(OutputFile);
  }
  else
  {
    if (Manager.GetBoolean("complex") == false)
    {
      LongRationalVector State (Space->GetHilbertSpaceDimension(), true);
      LongRational* Coefficients = InputFile.GetAsLongRationalArray(1);
      if (Coefficients == 0)
	{
	  InputFile.DumpErrors(cout) << endl;
	  return -1;
	}
      for (int i = 0; i < InputFile.GetNbrLines(); ++i)
	{
	  char* TmpString =  InputFile(0, i);
	  CleanLine(TmpString);
	  int TmpIndex = Space->FindStateIndex(TmpString);
	  if (TmpIndex != -1)
	    {
	      State[TmpIndex] = Coefficients[i];
	    }
	  else
	    {
	      cout << "warning , invalid state |" <<  TmpString << ">" << endl;
	    }
	}
      State.WriteVector(OutputFile);
    }
    else
    {
      ComplexVector State (Space->GetHilbertSpaceDimension(), true);
      double* CoefficientsRe = InputFile.GetAsDoubleArray(1);
      double* CoefficientsIm = InputFile.GetAsDoubleArray(2);
      Complex* Coefficients = new Complex[InputFile.GetNbrLines()];
      if (CoefficientsIm!=NULL)
	{
	  for (int i = 0; i < InputFile.GetNbrLines(); ++i)
	    Coefficients[i]=Complex(CoefficientsRe[i],CoefficientsIm[i]);
	}
      else
	{
	  for (int i = 0; i < InputFile.GetNbrLines(); ++i)
	    Coefficients[i]=Complex(CoefficientsRe[i]);
	 }
      delete[] CoefficientsRe;
      delete[] CoefficientsIm;
      if (Coefficients == 0)
	{
	  InputFile.DumpErrors(cout) << endl;
	  return -1;
	}
      for (int i = 0; i < InputFile.GetNbrLines(); ++i)
	{
	  char* TmpString =  InputFile(0, i);
	  CleanLine(TmpString);
	  int TmpIndex = Space->FindStateIndex(TmpString);
	  if (TmpIndex != -1)
	    {
	      State[TmpIndex] = Coefficients[i];
	    }
	  else
	    {
	      cout << "warning , invalid state |" <<  TmpString << ">" << endl;
	    }
	}
      if (Manager.GetBoolean("no-normalization") == false)
	State /= State.Norm();
      State.WriteVector(OutputFile);      
    }
  }

  return 0;
}

