#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Operator/ParticleOnSphereLMinusOperator.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "GeneralTools/ConfigurationParser.h"

#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"


#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("LMinusHaldane" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  
  Manager += SystemGroup;
  Manager += DataGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "twice the total lz value of the system for the initial state", 0);
  (*SystemGroup) += new SingleStringOption  ('s', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");
  (*SystemGroup) += new SingleStringOption  ('\n', "input-reference", "use a haldane basis with the given reference file for the input file");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-reference", "use a haldane basis with the given reference file for the output file (if omitted, full space is used)");


  (*DataGroup) += new SingleStringOption  ('i', "input-file", "input vector file name");
  (*DataGroup) += new SingleStringOption  ('o', "output-file", "output vector file name");

  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert-input", "save Hilbert space description for the input space in the indicated file and exit (only available for the non-symmetric Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert-input", "load Hilbert space description for the input space from the indicated file (only available for the non-symmetric Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert-output", "save Hilbert space description for the output space in the indicated file and exit (only available for the non-symmetric Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert-output", "load Hilbert space description for the output space from the indicated file (only available for the non-symmetric Haldane basis)",0);


  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type LMinus -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int TotalLz = ((SingleIntegerOption*) Manager["lz"])->GetInteger();
  bool FermionFlag = false;
  if (((SingleStringOption*) Manager["statistics"])->GetString() == 0)
    FermionFlag = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(((SingleStringOption*) Manager["input-file"])->GetString(), NbrParticles, LzMax, TotalLz, FermionFlag) == false)
    {
      return -1;
    }
  if ((((SingleStringOption*) Manager["statistics"])->GetString()) != 0)
    {
      if ((strcmp ("fermions", ((SingleStringOption*) Manager["statistics"])->GetString()) == 0))
	{
	  FermionFlag = true;
	}
      else
	if ((strcmp ("fermions", ((SingleStringOption*) Manager["statistics"])->GetString()) == 0))
	  {
	    FermionFlag = false;
	  }
	else
	  {
	    cout << ((SingleStringOption*) Manager["statistics"])->GetString() << " is an undefined statistics" << endl;
	  }  
    }
  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the same parity" << endl;
      return -1;           
    }

  RealVector InitialVector; 
  RealVector TargetVector; 
  if (InitialVector.ReadVector(Manager.GetString("input-file")) == false)
    {
      cout << "error while reading " << Manager.GetString("input-file") << endl;
      return -1;
    }
	
  long MemorySpace = 9l << 20;
  ParticleOnSphere* InputSpace=0;
  ParticleOnSphere* TargetSpace=0;
  cout << "Creating input space...";
  if (Manager.GetString("input-reference")!=NULL)
    {
      int* ReferenceState = 0;
      ConfigurationParser ReferenceStateDefinition;
      if (ReferenceStateDefinition.Parse(Manager.GetString("input-reference")) == false)
	{
	  ReferenceStateDefinition.DumpErrors(cout) << endl;
	  return 0;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	{
	  cout << "NbrParticles is not defined or as a wrong value" << endl;
	  return 0;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
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
      if (MaxNbrLz != (LzMax + 1))
	{
	  cout << "wrong LzMax value in ReferenceState" << endl;
	  return 0;     
	}
      if (FermionFlag == true)
	{
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	    if (LzMax <= 30)
#endif
	      {
		if (Manager.GetString("load-hilbert-input") != 0)
		  InputSpace = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert-input"), MemorySpace);
		else
		  InputSpace = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
		if (Manager.GetString("save-hilbert-input") != 0)
		  {
		    ((FermionOnSphereHaldaneBasis*) InputSpace)->WriteHilbertSpace(Manager.GetString("save-hilbert-input"));
		    return 0;
		  }
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	      if (LzMax <= 126)
#else
		if (LzMax <= 62)
#endif
		  {
		    if (Manager.GetString("load-hilbert-input") != 0)
		      InputSpace = new FermionOnSphereHaldaneBasisLong(Manager.GetString("load-hilbert-input"), MemorySpace);
		    else
		      InputSpace = new FermionOnSphereHaldaneBasisLong(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
		    if (Manager.GetString("save-hilbert-input") != 0)
		      {
			((FermionOnSphereHaldaneBasisLong*) InputSpace)->WriteHilbertSpace(Manager.GetString("save-hilbert-input"));
			return 0;
		      }
		  }
	}
      else
	{
#ifdef  __64_BITS__
	  if ((LzMax + NbrParticles - 1) < 63)
#else
	    if ((LzMax + NbrParticles - 1) < 31)	
#endif
	      InputSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, LzMax, ReferenceState);	  
	}
    }
  else
    {
      cout << "An input reference file is required!"<<endl;
      exit(1);
    }
  
  cout << "Parameters: d="<<InputSpace->GetHilbertSpaceDimension()<<", Lz="<<TotalLz<<"/2"<<endl;
    
    
  int TargetLz=TotalLz-2;
  cout << "Creating target space ...";
  if (Manager.GetString("output-reference")!=NULL)
    {
      int* ReferenceState = 0;
      ConfigurationParser ReferenceStateDefinition;
      if (ReferenceStateDefinition.Parse(Manager.GetString("output-reference")) == false)
	{
	  ReferenceStateDefinition.DumpErrors(cout) << endl;
	  return 0;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
	{
	  cout << "NbrParticles is not defined or as a wrong value" << endl;
	  return 0;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
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
      if (MaxNbrLz != (LzMax + 1))
	{
	  cout << "wrong LzMax value in ReferenceState" << endl;
	  return 0;     
	}
      if (FermionFlag == true)
	{
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	    if (LzMax <= 30)
#endif
	      {
		if (Manager.GetString("load-hilbert-output") != 0)
		  TargetSpace = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert-output"), MemorySpace);
		else
		  TargetSpace = new FermionOnSphereHaldaneBasis(NbrParticles, TargetLz, LzMax, ReferenceState, MemorySpace);
		if (Manager.GetString("save-hilbert-output") != 0)
		  {
		    ((FermionOnSphereHaldaneBasis*) TargetSpace)->WriteHilbertSpace(Manager.GetString("save-hilbert-output"));
		    return 0;
		  }
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	      if (LzMax <= 126)
#else
		if (LzMax <= 62)
#endif
		  {
		    if (Manager.GetString("load-hilbert-output") != 0)
		      TargetSpace = new FermionOnSphereHaldaneBasisLong(Manager.GetString("load-hilbert-output"), MemorySpace);
		    else
		      TargetSpace = new FermionOnSphereHaldaneBasisLong(NbrParticles, TargetLz, LzMax, ReferenceState, MemorySpace);
		    if (Manager.GetString("save-hilbert-output") != 0)
		      {
			((FermionOnSphereHaldaneBasisLong*) TargetSpace)->WriteHilbertSpace(Manager.GetString("save-hilbert-output"));
			return 0;
		      }
		  }
	}
      else
	{
#ifdef  __64_BITS__
	  if ((LzMax + NbrParticles - 1) < 63)
#else
	    if ((LzMax + NbrParticles - 1) < 31)	
#endif
	      TargetSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles, TargetLz, LzMax, ReferenceState);	  
	}
      cout << "Parameters: d="<<TargetSpace->GetHilbertSpaceDimension()<<", Lz="<<TargetLz<<"/2"<<endl;
      if (TargetLz!=TotalLz-2)
	{
	  cout << "Error: Wrong target angular momentum"<<endl;
	  exit(1);
	}
    }
  else
    {
      if (FermionFlag == true)
	{
#ifdef __64_BITS__
	  if (LzMax <= 63)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles, TargetLz, LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles, TargetLz, LzMax, MemorySpace);
	    }
#else
	  if (LzMax <= 31)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles, TargetLz, LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles, TargetLz, LzMax, MemorySpace);
	    }
#endif
	}
      else
	{
	  TargetSpace = new BosonOnSphereShort(NbrParticles, TargetLz, LzMax);
	}
      cout << "Full Hilbert space used: d="<< TargetSpace->GetHilbertSpaceDimension() << " Lz="<<TargetLz<<"/2"<<endl;
    }

  
  InputSpace->SetTargetSpace(TargetSpace);
  TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension(), true);
  if (TargetSpace->GetHilbertSpaceDimension()!=InputSpace->GetTargetHilbertSpaceDimension())
    {
      cout << "Problem with setting target space"<<endl;
      exit(-1);
    }
  ParticleOnSphereLMinusOperator LMinus(InputSpace, TotalLz, LzMax);
  LMinus.Multiply(InitialVector, TargetVector);
//   double Norm =TargetVector.Norm();
//   cout << "Target norm = "<<Norm <<endl;
  TargetVector/=TargetVector.Norm();

  char *OutputName;
  if (Manager.GetString("output-file")!=NULL)
    {
      OutputName = new char[strlen(Manager.GetString("output-file"))+1];
      strcpy(OutputName,Manager.GetString("output-file"));
    }
  else
    {
      OutputName = new char[strlen(Manager.GetString("input-file"))+10];
      sprintf(OutputName,"%s_L-",Manager.GetString("input-file"));
    } 
  if (TargetVector.WriteVector(OutputName) == false)
    {
      cout << "error while writing " << OutputName << endl;
      return -1;
    }
  delete InputSpace;
  delete TargetSpace;
  delete [] OutputName;
  return 0;
}

