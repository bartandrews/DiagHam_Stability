#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

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
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (only useful in su(2)/su(4) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "su3-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-tz", "twice the quantum number of the system associated to the Tz generator (only useful in su(3) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-y", "three time the quantum number of the system associated to the Y generator (only useful in su(3) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "su4-spin", "consider particles with SU(4) spin");
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin of the system (only usefull in su(4) mode)", 0);
  (*SystemGroup) += new SingleStringOption ('\0', "ascii-state", "name of the input ASCII description of the state (should use the same convention than FQHESphereShowBasis output)");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of trying to replace .txt extension with .vec or appending .vec extension");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereASCIIToState -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if (((SingleStringOption*) Manager["ascii-state"])->GetString() == 0)
    {
      cout << "An input file has to be provided" << endl;
      return -1;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrFluxQuanta = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger(); 
  int TotalLz = ((SingleIntegerOption*) Manager["lz-value"])->GetInteger();
  int TotalTz = ((SingleIntegerOption*) Manager["total-tz"])->GetInteger();
  int TotalY = ((SingleIntegerOption*) Manager["total-y"])->GetInteger();
  bool SU2SpinFlag = ((BooleanOption*) Manager["su2-spin"])->GetBoolean();
  bool SU3SpinFlag = ((BooleanOption*) Manager["su3-spin"])->GetBoolean();
  bool SU4SpinFlag = ((BooleanOption*) Manager["su4-spin"])->GetBoolean();
  int TotalSz = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
    
  if (((NbrParticles * NbrFluxQuanta) & 1) != (TotalLz & 1)) 
    {
      cout << "incompatible values for the number of particles, the number of flux quanta and twice the total lz value (nbr-particles * nbr-flux and lz-value should have the same parity)" << endl;
      return -1;
    }

  MultiColumnASCIIFile InputFile(':');
  if (InputFile.Parse(((SingleStringOption*) Manager["ascii-state"])->GetString()) == false)
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
  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
    {
      Space = new BosonOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
    }
  else
    {
      if ((SU2SpinFlag == false) && (SU3SpinFlag == false) && (SU4SpinFlag == false))
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
 	if (SU2SpinFlag == true)
	  Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
	else
 	if (SU3SpinFlag == true)
	  Space = new FermionOnSphereWithSU3Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalTz, TotalY);	  
    }

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
  State /= State.Norm();

  char* OutputFile = ((SingleStringOption*) Manager["output-file"])->GetString();
  if (OutputFile == 0)
    {
      OutputFile = ReplaceExtensionToFileName(((SingleStringOption*) Manager["ascii-state"])->GetString(), "txt", "vec");
      if (OutputFile == 0)
	OutputFile = AddExtensionToFileName(((SingleStringOption*) Manager["ascii-state"])->GetString(), "vec");
    }
  State.WriteVector(OutputFile);

  return 0;
}

