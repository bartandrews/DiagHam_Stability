#include "HilbertSpace/ParticleOnSphereManager.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereImportVector" , "0.01");

  ParticleOnSphereManager ParticleManager(true, false, 1);
  ParticleManager.AddOptionGroup(&Manager);
  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption ('\0', "state", "vector to import in ascii format");
  (*SystemGroup) += new SingleStringOption ('\n', "raw-state", "vector to import in FORTRAN binary format");
  (*SystemGroup) += new SingleStringOption ('b', "basis", "description of basis in which vector is formatted");
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*SystemGroup) += new SingleStringOption ('o', "output-state", "use this name for the output vector state instead of standard terminology");
  (*SystemGroup) += new BooleanOption  ('r', "reverse", "reverse order of states");
  (*SystemGroup) += new BooleanOption  ('e', "export", "export the into ascii format");
  (*MiscGroup) += new BooleanOption  ('v', "verbose", "give a lot of output");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int Lz = Manager.GetInteger("lz-value");
  bool Verbose = Manager.GetBoolean("verbose");

  ifstream InputFile;
  if ((Manager.GetString("state")!=0)&&(!Manager.GetBoolean("export")))
    {
      InputFile.open(Manager.GetString("state"), ios::in);
      if (!InputFile.is_open())
	{
	  cout << "Could not open text-file with vector description!"<<endl;
	  exit(-1);	    
	}
    }
  else
    {
      if (Manager.GetBoolean("export"))
	{
	  InputFile.open(Manager.GetString("state"), ios::binary | ios::in);
	  if (!InputFile.is_open())
	    {
	      cout << "Could not open input vector!"<<endl;
	      exit(-1);	    
	    }
	}
      else
	{
	  if (Manager.GetString("raw-state")!=0)
	    {
	      InputFile.open(Manager.GetString("raw-state"), ios::binary | ios::in);
	      if (!InputFile.is_open())
		{
		  cout << "Could not open binary-file with vector description!"<<endl;
		  exit(-1);	    
		}      
	    }
	  else
	    {
	      cout << "Require a state vector in ascii or raw format!"<<endl;
	      exit(1);
	    }
	}
    }


  ParticleOnSphere *Space = ParticleManager.GetHilbertSpace(Lz);

  cout << "Hilbert-space dimension: "<<Space->GetHilbertSpaceDimension() << endl;

  int *Map = new int[Space->GetHilbertSpaceDimension()];
  int *Signs = NULL; 

  if (Manager.GetString("basis") != 0)
    {
      Signs = new int[Space->GetHilbertSpaceDimension()];
      ifstream BasisFile;
      BasisFile.open(Manager.GetString("basis"), ios::in);
      if (!BasisFile.is_open())
	{
	  cout << "Could not open text-file with basis description!"<<endl;
	  exit(-1);	    
	}
      int Index, LastIndex=-1;
      char Descriptor[256];
      while (!BasisFile.eof())
	{
	  BasisFile >> Index >> Descriptor;
	  --Index; // convert to range 0..dim-1
	  if ((Index>LastIndex)&&(Index<Space->GetHilbertSpaceDimension()))
	    {
	      
	      if ((strlen(Descriptor)==(unsigned)LzMax+1))
		{
		  // recognized format used by Edward Rezayi for bilayer basis:
		  int pos=0, NToBeFound=NbrParticles;
		  unsigned long State=0x0l;
		  double Coefficient=1.0, TmpCoeff;
		  if (Verbose)
		    cout << Descriptor << " : " << endl;
		  while ((pos<LzMax+1)&&(NToBeFound>0))
		    {
		      if (Descriptor[pos]=='1')
			{			  
			  State = ((FermionOnSphere*)Space)->Ad(State, pos, TmpCoeff);
			  Coefficient*=TmpCoeff;
			  --NToBeFound;
			}
		      ++pos;
		    }
		  if (NToBeFound!=0)
		    {
		      cout << "Wrong number of particles - check your entries!"<<endl;
		      exit(-1);
		    }
		  Map[Index]=((FermionOnSphere*)Space)->CarefulFindStateIndex(State,-1);
		  Signs[Index]=Coefficient;
		}
	      else
		{
		  // could define other bases here...
		  cout << "LzMax does not correspond to given basis - please check your entry!"<<endl;
		  exit(-1);
		}
	      if (Verbose)
		{
		  cout << Index << "\t" << Descriptor << "\t->"<<Signs[Index]<<"*'"<<Map[Index]<<"' : ";
		  Space->PrintState(cout, Map[Index]) << endl;
		}
	    }	  
	  LastIndex=Index;	  
	}
      if (LastIndex!=Space->GetHilbertSpaceDimension()-1)
	{
	  cout << "Attention: dimensions do not match!"<<endl;
	  exit(1);
	}
      else 
	{
	  cout << "Basis read successfully!"<<endl;
	}
      BasisFile.close();
      for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
	if (Map[i]==-1)
	  {
	    cout << "Basis not complete!"<<endl;
	    exit(-1);
	  }
    }
  else
    {
      int Dim=Space->GetHilbertSpaceDimension();
      if (Manager.GetBoolean("reverse"))
	for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i) Map[i]=Dim-1-i;
      else
	for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i) Map[i]=i;
    }

  if (Manager.GetBoolean("export"))
    {
      InputFile.close();
      RealVector InputVector;
      if (InputVector.ReadVector(Manager.GetString("state"))==false)
	{
	  cout << "Could not read input vector"<<endl;
	  exit(1);
	}
      char *OutputName = Manager.GetString("output-state");
      if (OutputName==NULL)
	{
	  OutputName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "txt");
	  if (OutputName==NULL)
	    OutputName = AddExtensionToFileName(Manager.GetString("state"), "txt");
	}
      ofstream OutputFile;
      OutputFile.open(OutputName,ios::out);
      OutputFile.precision(15);
      for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
	OutputFile<<InputVector[Map[i]]<<endl;
      OutputFile.close();
      delete [] OutputName;
    }
  else // regular import
    {
      RealVector ResultingVector(Space->GetLargeHilbertSpaceDimension());

      char *OutputName = Manager.GetString("output-state");
  
      if (OutputName==NULL)
	{
	  OutputName = new char[512];
	  if (Manager.GetString("state")!=NULL)
	    {
	      char *Path, *FileName;
	      ExtractPathAndFileName (Manager.GetString("state"), Path, FileName);
	      sprintf(OutputName,"fermions_sphere_import-%s_n_%d_2s_%d_lz_%d.vec",
		      FileName, NbrParticles, LzMax, Lz);
	      delete [] Path;
	      delete [] FileName;
	    }
	  else if (Manager.GetString("raw-state")!=NULL)
	    {
	      char *Path, *FileName;
	      ExtractPathAndFileName (Manager.GetString("raw-state"), Path, FileName);
	      sprintf(OutputName,"fermions_sphere_import-%s_n_%d_2s_%d_lz_%d.vec",
		      FileName, NbrParticles, LzMax, Lz);
	      delete [] Path;
	      delete [] FileName;
	    }
	}
      else
	{
	  OutputName = new char[strlen(Manager.GetString("output-state")+1)];
	  strcpy(OutputName,Manager.GetString("output-state"));
	}

      if (Manager.GetString("state")!=0)
	{
	  // int VIndex, VLastIndex=-1;
	  char NextLine[256];
	  double Value, Sign=1.0;
	  int Index=0;
	  while (InputFile.getline(NextLine,256))
	    {
	      if (Signs)
		Sign = Signs[Index];
	      std::istringstream LineStream(NextLine);
	      LineStream >> Value;
	      if (Verbose)
		cout << Index << "\t" << Value << endl;
	      ResultingVector[Map[Index]]=Sign*Value;
	      Index++;
	    }
	}
      else if (Manager.GetString("raw-state")!=0)
	{
	  double *InputVector = new double[Space->GetHilbertSpaceDimension()];
	  int header;
	  InputFile.read ((char*)&header, sizeof(int));      
	  if (header/sizeof(double)!=Space->GetHilbertSpaceDimension())
	    {
	      cout<<"Apparent problem with reading input vector"<<endl;
	      exit(1);
	    }	  
	  InputFile.read ((char*)InputVector, Space->GetHilbertSpaceDimension()*sizeof(double));
	  for (int i=0; i<10; ++i) cout << "Input["<<i<<"]="<<InputVector[i]<<endl;
	  if (Signs!=NULL)
	    for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
	      ResultingVector[Map[i]]=Signs[i]*InputVector[i];
	  else
	    for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
	      ResultingVector[Map[i]]=InputVector[i];
	}
      else
	{
	  cout << "An error occurred"<<endl;
	  exit(1);
	}
      cout << "Writing vector "<<OutputName<<endl;
      ResultingVector.WriteVector(OutputName);
      delete [] OutputName;
    }
}
