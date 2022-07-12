#include "Vector/RealVector.h"

#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
#include "HilbertSpace/BosonOnSphereShort.h"


#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

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
  OptionManager Manager ("FQHESphere2LLTimes2LLProjection" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax1", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('z', "totallz1", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
	(*SystemGroup) += new SingleIntegerOption  ('\n', "lzmax2", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "totallz2", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "express the projected states in the normalized basis");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphere2LLTimes2LLProjection -h" << endl;
      return -1;
    }	
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  
  int NbrFermion = Manager.GetInteger("nbr-particles");
  int LzMax1 = Manager.GetInteger("lzmax1");
  int TotalLz1 = Manager.GetInteger("totallz1");
	int LzMax2 = Manager.GetInteger("lzmax2");
  int TotalLz2 = Manager.GetInteger("totallz2");
	
  
  int Parity = TotalLz1 & 1;
  if (Parity != ((NbrFermion * LzMax1) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;
    }
    
    int Parity2 = TotalLz2 & 1;
  if (Parity2 != ((NbrFermion * LzMax2) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;
    }
  

  
  FermionOnSphereTwoLandauLevels * SpaceLL1 = new FermionOnSphereTwoLandauLevels (NbrFermion, TotalLz1, LzMax1 + 2, LzMax1);
  FermionOnSphereTwoLandauLevels * SpaceLL2 = new FermionOnSphereTwoLandauLevels (NbrFermion, TotalLz2, LzMax2 + 2, LzMax2);
  BosonOnSphereShort * FinalSpace = new BosonOnSphereShort (NbrFermion, TotalLz1 + TotalLz2, LzMax1 + LzMax2);
  
	      
	
  char * PreffixFileName = new char [50];
  char * OutputFileName = new char [60];
	sprintf (PreffixFileName, "bosons_%s_n_%d_2s_%d_lz_%d",Manager.GetString("interaction-name"), NbrFermion,LzMax1 + LzMax2, TotalLz1 + TotalLz2);
	
	RealVector * FermionState1 = new RealVector(SpaceLL1->GetHilbertSpaceDimension(),true);
	RealVector * FermionState2 = new RealVector(SpaceLL2->GetHilbertSpaceDimension(),true);
  RealVector * OutputVector =  new RealVector(FinalSpace->GetHilbertSpaceDimension(),true);
   int Index = 0;   
  for (int i = 0; i< SpaceLL1->GetHilbertSpaceDimension(); 	i++)
	{
		(*FermionState1)[i]=1;
		for (int j = 0; j< SpaceLL2->GetHilbertSpaceDimension(); 	j++)
	{
		(*FermionState2)[j]=1;
		
		SpaceLL1->FermionicStateTimeFermionicState( (*FermionState1), (*FermionState2), (*OutputVector), SpaceLL2 , FinalSpace, 0,SpaceLL2->GetHilbertSpaceDimension());
		FinalSpace->ConvertFromUnnormalizedMonomial((*OutputVector),0,true);
		
	  sprintf (OutputFileName,"%s.%d.vec",PreffixFileName,Index);
	  (*OutputVector).WriteVector(OutputFileName);
		Index++;
		(*FermionState2)[j]=0;
	}
	(*FermionState1)[i]=0;
	}
	  return 0;
}
	
      