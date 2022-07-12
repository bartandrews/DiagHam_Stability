#include "config.h"
#include "Vector/RealVector.h"
#include "GeneralTools/Endian.h"

#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"
#include "Options/Options.h"


#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereJastrowMultiplicationOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereJastrowDivisionOperation.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <cstring>
#include <stdio.h>



using std::cout;
using std::endl;

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereBosonsFermionsConverter" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  ArchitectureManager Architecture;
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  Manager += OutputGroup;
 
  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "name of the file that contains the bosonic(fermionic) state that will be multiplied (divided) by a Jastrow Factor");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new BooleanOption  ('\n', "p-truncated", "use a p-truncated basis instead of the full n-body basis");
  (*SystemGroup) += new SingleIntegerOption ('\n', "p-truncation", "p-truncation for the p-truncated basis (if --p-truncated is used)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "boson-truncation", "maximum occupation for a given orbital for the p-truncated basis", 1);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state (should be the one of the bosonic state)");
  (*SystemGroup) += new SingleStringOption  ('\n', "vector-file", "single column file describing a list of state belonging to the same hilbert space to be converted");
  (*OutputGroup) += new MultipleStringOption ('o', "output-states", "output file name (if none, guess it from the input file name)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsFermionsConverter -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  
  int NbrVectors;
  char** VectorFiles = Manager.GetStrings("states", NbrVectors);
  
  if((NbrVectors == 0)&&(Manager.GetString("vector-file")==0))
    {
      cout << "no input state" << endl << "see man page for option syntax or type FQHESphereBosonsFermionsConverter -h" << endl;
      return -1;
    }
  
  if((Manager.GetString("vector-file")!=0)&&(IsFile(Manager.GetString("vector-file"))==false))
    {
      cout << "can't open file " <<  Manager.GetString("vector-file")<< endl;
      return -1;
    }
  
  int NbrSpaces=0;
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  char** StateVectorFiles=0;
  char** StateVector=0;
  
  if (Manager.GetString("vector-file") != 0)
    {
      MultiColumnASCIIFile VectorFile;
      if (VectorFile.Parse(Manager.GetString("vector-file")) == false)
	{
	  VectorFile.DumpErrors(cout);
	  return -1;
	}
      NbrSpaces = VectorFile.GetNbrLines();
      StateVectorFiles = new char* [NbrSpaces];
      
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  StateVectorFiles[i] = new char [strlen(VectorFile(0, i)) + 1];
	  strcpy (StateVectorFiles[i], VectorFile(0, i));   
	}
    }
  StateVector = new char * [NbrSpaces+NbrVectors];
  
  for(int i = 0;i<NbrVectors;i++)
    {
      StateVector[i] = new char [strlen(VectorFiles[i]) + 1];
      strcpy (StateVector[i], VectorFiles[i]);
    }
  for(int i = 0;i<NbrSpaces;i++)
    {
      StateVector[i+NbrVectors]=new char [strlen(StateVectorFiles[i]) + 1];
      strcpy (StateVector[i+NbrVectors], StateVectorFiles[i]);
    }
  
  NbrVectors += NbrSpaces;
  cout <<"NbrVectors = "<< NbrVectors<<endl;
  
  int* NbrParticles = new int [NbrVectors];
  int* LzMax = new int [NbrVectors];
  int* TotalLz = new int [NbrVectors];
  bool* FermionFlag = new bool [NbrVectors];
  
  for(int i = 0; i < NbrVectors; i++)
    {
      NbrParticles[i] = 0;
      LzMax[i] = 0;
      TotalLz[i] = 0;
      FermionFlag[i] = true;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(StateVector[i], NbrParticles[i], LzMax[i], TotalLz[i], FermionFlag[i]) == false)
	{
	  return -1;
	}
      if (FermionFlag[i] == true)
	LzMax[i] -= NbrParticles[i] - 1;
    }
  
  for(int i = 1; i < NbrVectors; i++)
    {
      if((NbrParticles[i] != NbrParticles[i - 1]) || (LzMax[i] != LzMax[i - 1]) || (TotalLz[i] != TotalLz[i - 1]) || 
	 (FermionFlag[i] != FermionFlag[i - 1]))
	{ 
	  cout << "Vectors are not all in the same HilbertSpace." << endl;
	  return -1;
	}
    }
  
  int Parity = TotalLz[0] & 1;
  if (Parity != ((NbrParticles[0] * LzMax[0]) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;
    }
  
  int NbrOutput;
  char** VectorFilesOut = Manager.GetStrings("output-states", NbrOutput);
  if (NbrOutput != 0)
    {
      if (NbrVectors != NbrOutput)
	{
	  cout << "The number of vectors in entry must be the same that the number of Output vectors." << endl;
	  return -1;
	}
    }
  else
    {
      VectorFilesOut = new char*[NbrVectors];
      for(int i = 0; i < NbrVectors; i++)
	{
	  VectorFilesOut[i] = new char [strlen (StateVector[i])+ 32];
	  char* TmpPos = VectorFilesOut[i];
	  if (FermionFlag[i] == true)
	    {
	      char* TmpPos2 = strstr(StateVector[i], "_n_");
	      sprintf (TmpPos, "bosons_");
	      TmpPos += 7;
	      strncpy (TmpPos, StateVector[i] + 9, TmpPos2 - StateVector[i] - 9);
	      TmpPos += TmpPos2 - StateVector[i] - 9;
	      TmpPos2 = strstr(StateVector[i], "_lz_");
	      sprintf (TmpPos, "_n_%d_2s_%d%s", NbrParticles[i], LzMax[i], TmpPos2);
	    }
	  else
	    {
	      char* TmpPos2 = strstr(StateVector[i], "_n_");
	      sprintf (TmpPos, "fermions_");
	      TmpPos += 9;
	      strncpy (TmpPos, StateVector[i] + 7, TmpPos2 - StateVector[i] - 7);
	      TmpPos += TmpPos2 - StateVector[i] - 7;
	      TmpPos2 = strstr(StateVector[i], "_lz_");
	      sprintf (TmpPos, "_n_%d_2s_%d%s", NbrParticles[i], LzMax[i] + NbrParticles[i] - 1, TmpPos2);
	    }
	}
    }

  RealVector * InputState = new RealVector[NbrVectors];
  
  for(int i = 0; i < NbrVectors; i++)
    {
      if (IsFile(StateVector[i]) == false)
	{
	  cout << "state " << StateVector[i] << " does not exist or can't be opened" << endl;
	  return -1;
	}
      if(InputState[i].ReadVector(StateVector[i]) == false)
	{
	  cout << "error while reading " << StateVector[i] << endl;
	  return -1;
	}
    }
  
  BosonOnSphereShort* Space = 0;
  
  if ((HaldaneBasisFlag == false) && (Manager.GetBoolean("p-truncated") == false))
    {
#ifdef  __64_BITS__
      if ((LzMax[0] + NbrParticles[0] - 1) < 63)
#else
	if ((LzMax[0] + NbrParticles[0] - 1) < 31)	
#endif
	  {
	    Space = new BosonOnSphereShort (NbrParticles[0], TotalLz[0], LzMax[0]);
	  }
	else
	  {
		cout << "error, the needed space needed class BosonOnSphere" <<endl;
	  }
    }
  else
    {
      int* ReferenceState = 0;
      if (Manager.GetString("reference-file") == 0)
	{
	  cout << "error, a reference file is needed for bosons in Haldane basis" << endl;
	  return -1;
	}
      ConfigurationParser ReferenceStateDefinition;
      if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
	{
	  ReferenceStateDefinition.DumpErrors(cout) << endl;
 	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles[0]) == false) || (NbrParticles[0] <= 0))
	{
	  cout << "NbrParticles is not defined or as a wrong value" << endl;
	  return -1;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax[0]) == false) || (LzMax[0] < 0))
	{
	  cout << "LzMax is not defined or as a wrong value" << endl;
	  return 0;
	}
      int MaxNbrLz;
      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	{
	  cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
	  return -1;     
	}
      if (MaxNbrLz != (LzMax[0] + 1))
	{
	  cout << "wrong LzMax value in ReferenceState" << endl;
	  return -1;
	}
#ifdef  __64_BITS__
      if (LzMax[0]  < 63)
#else
	if (LzMax[0]  < 31)	
#endif
	  {
	    if (Manager.GetBoolean("p-truncated") == true)
	      {
		Space = new BosonOnSpherePTruncated(NbrParticles[0], TotalLz[0], LzMax[0], Manager.GetBoolean("p-truncation"), (int) Manager.GetInteger("boson-truncation"), ReferenceState);	  
	      }
	    else
	      Space = new BosonOnSphereHaldaneBasisShort(NbrParticles[0], TotalLz[0], LzMax[0], ReferenceState);
	  }
    }
  
  for(int i = 0; i < NbrVectors; i++)
    {
      if (Space->GetHilbertSpaceDimension() != InputState[i].GetVectorDimension())
	{
	  cout << "dimension mismatch between the state (" << InputState[i].GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}
    }
  
  
  RealVector* OutputVector = new RealVector[NbrVectors];
  for (int i = 0; i < NbrVectors; i++)
    OutputVector[i] = RealVector(Space->GetHilbertSpaceDimension(), true);
  
  if (FermionFlag[0] == false)
    {
      FQHESphereJastrowMultiplicationOperation Operation(Space, InputState, OutputVector, NbrVectors);
      Operation.ApplyOperation(Architecture.GetArchitecture());
    }
  else
    {
      FQHESphereJastrowDivisionOperation Operation(Space, InputState, OutputVector, NbrVectors);
      Operation.ApplyOperation(Architecture.GetArchitecture());
    }


  for(int k = 0; k < NbrVectors; k++)
    {
      cout << "writing " << VectorFilesOut[k] << endl;
      OutputVector[k].WriteVector(VectorFilesOut[k]);
    }
 
  return 0;
}
