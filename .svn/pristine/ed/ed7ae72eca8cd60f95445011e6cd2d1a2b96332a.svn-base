#include "Vector/RealVector.h"

#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"

#include "Options/Options.h"
#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"                      
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"


#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  
  OptionManager Manager ("FQHESphere2LLBosonicStateTimePolarizedSlaters" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  ArchitectureManager Architecture;
	
  Manager += SystemGroup;
  Manager += MiscGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += OutputGroup;
	
  (*SystemGroup) += new SingleStringOption ('\0', "state", "name of the vector file in the 2LL");
  (*SystemGroup) += new BooleanOption  ('\n', "2-ll", "consider particles within two Landau levels");    
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the Jack polynomial decomposition into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the Jack polynomial decomposition into a text file");
  (*SystemGroup) += new SingleIntegerOption  ('s', "sz", "value of sz in the final space",0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "value of lz in the final space",0);

  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files, default is nbody)");
	// (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  // (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  // (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleStringOption  ('\n', "resume-file", "use this file as the partial vector to resume from");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "resume-idx", "use this file as the partial vector to resume from", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");  
  
  (*SystemGroup) += new BooleanOption ('\n', "2-ll-lz", "use lz symmetry to reduce number of 2-ll boson configs taken",false);
  (*SystemGroup) += new BooleanOption ('\n', "2-ll-sz", "use sz symmetry to reduce number of 2-ll boson configs taken",false);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "mpi-stages", "the number of stages divide into when using MPI  (default is 20)", 20);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "smp-stages", "the number of stages divide into when using SMP  (default is 20)", 20);
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "the output vector will be normalize on the factory",false);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsTimesFermions -h" << endl;
      return -1;
    }
	
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
 
  bool LL2 = Manager.GetBoolean("2-ll"); 
  bool LzSym = Manager.GetBoolean("2-ll-lz");   
  bool SzSym = Manager.GetBoolean("2-ll-sz");  
	
  RealVector InitialState;
  if (InitialState.ReadVector(Manager.GetString("state")) == false)
    {
      cout << "error while reading " << Manager.GetString("state") << endl;
      return -1;
    }
 
  int NbrParticles = 0;
  int LzMax = 0;
  int TotalLz = 0;
  int FinalLz = Manager.GetInteger("lz");
  int TotalSz = Manager.GetInteger("sz");
  bool FermionFlag = false;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMax, TotalLz, FermionFlag) == false)
   {
     return -1;
   }
 
   ParticleOnSphere * InitialSpace;
   if ( LL2 == true)
     InitialSpace = new BosonOnSphereTwoLandauLevels (NbrParticles,TotalLz,LzMax+ 2,LzMax);
   else
     InitialSpace = new BosonOnSphereShort (NbrParticles,TotalLz,LzMax);
	
  if (InitialSpace->GetHilbertSpaceDimension() != InitialState.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << InitialState.GetVectorDimension() << ") and the Hilbert space (" << InitialSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
    
  if ( LL2 == true)
    {
      int OldDimension = InitialSpace->GetHilbertSpaceDimension();           
      int NewDimension = ((BosonOnSphereTwoLandauLevels*)InitialSpace)->RemoveZeros(InitialState, LzSym);      
      cout << "Removing zero valued elements: " << OldDimension << " -> " << NewDimension << endl;
    }
    
    FermionOnSphere * SlaterSpace = 0;
    FermionOnSphere * SlaterSpaceUp = 0;
    FermionOnSphere * SlaterSpaceDown = 0;
    int NbrFermionsUp = (NbrParticles + TotalSz)/2;
    int NbrFermionsDown = (NbrParticles - TotalSz)/2;
    
    if (TotalSz == 0)
      {
	SlaterSpace = new FermionOnSphere (NbrFermionsUp,0,NbrFermionsUp-1);
      }
    else
      {
	SlaterSpaceUp = new FermionOnSphere (NbrFermionsUp,0,NbrFermionsUp-1);
	SlaterSpaceDown = new FermionOnSphere (NbrFermionsDown,FinalLz - TotalLz, NbrFermionsUp - 1);
      }

	
  FermionOnSphereWithSpin * FinalSpace;
  if ( Manager.GetBoolean("haldane") == false)
    {
       FinalSpace = new FermionOnSphereWithSpin (NbrParticles,FinalLz, LzMax + NbrFermionsUp - 1,TotalSz);
    }
  else if (Manager.GetString("reference-file") != 0 )
    {
      int** ReferenceStates = 0;
      int NbrReferenceStates;		    
      bool TexturelessFlag;
      if (FQHEGetRootPartitionSU2(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceStates, NbrReferenceStates, TexturelessFlag) == false)
	{
	  cout << "error while parsing " << Manager.GetString("reference-file") << endl;	      
	  return 0;
	}      		  		         
      if (TexturelessFlag == false ) 
	{	
	  if ( Manager.GetBoolean("lzsymmetrized-basis") && Manager.GetBoolean("szsymmetrized-basis") )
	    {		      
	      bool LzMinusParity = ((BooleanOption*) Manager["minus-lzparity"])->GetBoolean();      
	      bool SzMinusParity = ((BooleanOption*) Manager["minus-szparity"])->GetBoolean();	  	  
	      FinalSpace = new FermionOnSphereWithSpinHaldaneLzSzSymmetry(NbrParticles, LzMax, SzMinusParity, LzMinusParity, ReferenceStates, NbrReferenceStates);
	    }
	  else
	    {
	      FinalSpace = new FermionOnSphereWithSpinHaldaneBasis(NbrParticles, TotalLz, LzMax, TotalSz, ReferenceStates, NbrReferenceStates);
	    }	    
	}
      else
	{			  			    
	  int **TexturelessReferenceState = new int*[NbrReferenceStates];
	  for ( int j = 0 ; j < NbrReferenceStates ; j++ ) 
	    {
	      TexturelessReferenceState[j] = new int[LzMax+1];
	      for ( int i = 0 ; i < (LzMax + 1) ; i++ )
		{
		  if ( ReferenceStates[j][i] == 3 ) 
		    {
			TexturelessReferenceState[j][i] = 2;
		    }
		  else if ( (ReferenceStates[j][i] == 1) || (ReferenceStates[j][i] == 2) ) 
		    {
			TexturelessReferenceState[j][i] = 1;
		    }
		  else
		    {
			TexturelessReferenceState[j][i] = 0;
		    }
		}
	    }	
	  if ( Manager.GetBoolean("lzsymmetrized-basis") && Manager.GetBoolean("szsymmetrized-basis") )
	    {		      
	      bool LzMinusParity = ((BooleanOption*) Manager["minus-lzparity"])->GetBoolean();      
	      bool SzMinusParity = ((BooleanOption*) Manager["minus-szparity"])->GetBoolean();
	      FinalSpace = new FermionOnSphereWithSpinHaldaneLzSzSymmetry(NbrParticles, LzMax, SzMinusParity, LzMinusParity, TexturelessReferenceState, NbrReferenceStates, true);
	    }
	  else
	    {
	      FinalSpace = new FermionOnSphereWithSpinHaldaneBasis(NbrParticles, TotalLz, LzMax, TotalSz, TexturelessReferenceState, NbrReferenceStates, true);
	    }
	}
    }    
  else
    {
      cout << "error, no reference file." << endl;
      return 0;
    }
    
  RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
	int ResumeIdx = Manager.GetInteger("resume-idx");
		
	if (Manager.GetString("resume-file") != 0 )
		{
			if ( Architecture.GetArchitecture()->ReadVector(OutputVector, Manager.GetString("resume-file")) == false )
			{
				cout << "error while reading " << Manager.GetString("resume-file") << endl;
				return -1;
			}		
		}
		
	if(TotalSz == 0)
	{
	
		
		FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation Operation(InitialSpace, SlaterSpace, FinalSpace, &InitialState, &OutputVector,LL2,LzSym, SzSym,Manager.GetInteger("mpi-stages"), Manager.GetInteger("smp-stages"), ResumeIdx);	
		Operation.ApplyOperation(Architecture.GetArchitecture());    
  

		if(Manager.GetBoolean("normalize"))
			FinalSpace->ConvertFromUnnormalizedMonomial(OutputVector,0l,true);
		
		Architecture.GetArchitecture()->WriteVector(OutputVector, Manager.GetString("bin-output"));
	}
	else
	{
	  char* EigenvectorName = 0;
	  if (((SingleStringOption*) Manager["interaction-name"])->GetString() == 0)
	    {
	      EigenvectorName = new char [256];
	      sprintf (EigenvectorName, "fermions_su2_product_n_%d_2s_%d_sz_%d_lz_%d",  NbrParticles, LzMax + NbrFermionsUp - 1 ,TotalSz, FinalLz);
	    }
	  else
	    {
	      EigenvectorName = new char [256 + strlen(((SingleStringOption*) Manager["interaction-name"])->GetString())];
	      sprintf (EigenvectorName, "fermions_su2_%s_n_%d_2s_%d_sz_%d_lz_%d", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax +NbrFermionsUp - 1 ,TotalSz, FinalLz);
	    }
	    
	    char *	TmpVectorName = new char [strlen(EigenvectorName) + 16];
		int NbrState = 0;
			 for (int IndexUp = 0; IndexUp <  SlaterSpaceUp->GetHilbertSpaceDimension(); IndexUp++)
 {
			for (int IndexDown = 0; IndexDown <  SlaterSpaceDown->GetHilbertSpaceDimension(); IndexDown++)
 {
		sprintf (TmpVectorName, "%s.%d.vec", EigenvectorName, NbrState);
		
		FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation Operation(InitialSpace, SlaterSpaceUp , SlaterSpaceDown, FinalSpace, &InitialState, &OutputVector,IndexUp, IndexDown , true, false, false, Manager.GetInteger("mpi-stages"), Manager.GetInteger("smp-stages"), ResumeIdx);
		
		Operation.ApplyOperation(Architecture.GetArchitecture()); 
		if(Manager.GetBoolean("normalize"))
			FinalSpace->ConvertFromUnnormalizedMonomial(OutputVector,0l,true);
		Architecture.GetArchitecture()->WriteVector(OutputVector, TmpVectorName);
		OutputVector.ClearVector();
		NbrState++;
 }
}
	}	
}


