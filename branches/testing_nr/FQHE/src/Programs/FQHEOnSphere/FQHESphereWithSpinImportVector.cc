#include "HilbertSpace/ParticleOnSphereManager.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

#include "Options/Options.h"

#include "Vector/RealVector.h"

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <bitset>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::bitset;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereWithSpinImportVector" , "0.01");

  ParticleOnSphereManager ParticleManager(true, false, 2);
  ParticleManager.AddOptionGroup(&Manager);
  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = Manager.GetOptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption ('\0', "state", "vector to import in ascii format");
  (*SystemGroup) += new SingleStringOption ('b', "basis", "description of basis in which vector is formatted (none = standard basis)");
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-state", "use this name for the output vector state instead of standard terminology");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int Lz = Manager.GetInteger("lz-value");
  int Sz = Manager.GetInteger("total-sz");

  ifstream InputFile;
  InputFile.open(Manager.GetString("state"), ios::in);
  if (!InputFile.is_open())
    {
      cout << "Could not open text-file with vector description!"<<endl;
      exit(-1);	    
    }  

  ParticleOnSphere *Space = ParticleManager.GetHilbertSpace(Lz);

  cout << "Hilbert-space dimension: "<<Space->GetHilbertSpaceDimension() << endl;   

  int *Map = new int[Space->GetHilbertSpaceDimension()];
  double *Signs = new double[Space->GetHilbertSpaceDimension()];
  if (Manager.GetString("basis") != 0)
    {
      for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
      {
	Map[i]=-1;
      }
      ifstream BasisFile;
      BasisFile.open(Manager.GetString("basis"), ios::in);
      if (!BasisFile.is_open())
	{
	  cout << "Could not open text-file with basis description!"<<endl;
	  exit(-1);	    
	}
      int Index, LastIndex=-1;
      char Descriptor[256];
      bitset<20> b;
      while (BasisFile.getline(Descriptor,256))
	{
	  std::istringstream LineStream(Descriptor);
	  LineStream >> Index >> Descriptor;	  
	  --Index; // convert to range 0..dim-1
	  if ((Index>LastIndex)&&(Index<Space->GetHilbertSpaceDimension()))
	    {
	      
	      if ((strlen(Descriptor)==2*(unsigned)LzMax+3)&&(Descriptor[LzMax+1]==','))
		{
		  // recognized format used by Edward Rezayi for bilayer basis:
		  int pos=0, Nup=(NbrParticles+Sz)/2, Ndown=(NbrParticles-Sz)/2;
		  unsigned int State=0x0l;
		  double Coefficient=1.0, TmpCoeff;
		  cout << Descriptor << " : " << endl;
		  while ((pos<LzMax+1)&&(Nup>0))
		    {
		      if (Descriptor[pos]=='1')
			{			  
			  State = ((ParticleOnSphereWithSpin*)Space)->Ad(State, pos, 1, TmpCoeff);
			  b=State;
			  cout << "up " << Nup << " " << b << endl;
			  Coefficient*=TmpCoeff;
			  --Nup;
			}
		      ++pos;
		    }
		  if (Nup!=0)
		    {
		      cout << "Wrong number of particles, or total spin - check your entries!"<<endl;
		      exit(-1);
		    }
		  pos=0;
		  while ((pos<LzMax+1)&&(Ndown>0))
		    {
		      if (Descriptor[LzMax+2+pos]=='1')
			{
			  State = ((ParticleOnSphereWithSpin*)Space)->Ad(State, pos, 0, TmpCoeff);
			  b=State;
			  cout << "down " << Ndown << " " << b << endl;
			  Coefficient*=TmpCoeff;
			  --Ndown;
			}
		      ++pos;
		    }
		  if (Ndown!=0)
		    {
		      cout << "Wrong number of particles, or total spin - check your entries!"<<endl;
		      exit(-1);
		    }
		  Map[Index]=((ParticleOnSphereWithSpin*)Space)->CarefulFindStateIndex(State,-1);
		  Signs[Index]=Coefficient;
		}
	      else
		{
		  // could define other bases here...
		  cout << "LzMax does not correspond to given basis - please check your entry!"<<endl;
		  exit(-1);
		}
	      cout << Index << "\t" << Descriptor << "\t->"<<Signs[Index]<<"*'"<<Map[Index]<<"' : ";
	      Space->PrintState(cout, Map[Index]) << endl;
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
    for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
      {
	Map[i]=i;
	Signs[i]=1.0;
      }

  RealVector ResultingVector(Space->GetHilbertSpaceDimension());

  char *OutputName = Manager.GetString("output-state");

  if (OutputName==NULL)
    {
      OutputName = new char[512];
      sprintf(OutputName,"fermions_sphere_spin_import-%s_n_%d_2S_%d_Sz_%d_lz_%d.vec",
	      Manager.GetString("state"), NbrParticles, LzMax, Sz, Lz);      
    }
  
  int VIndex, VLastIndex=-1;
  char NextLine[256];
  double Value;
  int Index;
  while (InputFile.getline(NextLine,256))
    {
      std::istringstream LineStream(NextLine);
      LineStream >> Index >> Value;
      cout << Index << "\t" << Value << endl;
      ResultingVector[Map[Index-1]]=Signs[Index-1]*Value;
    }

  ResultingVector.WriteVector(OutputName);

  delete [] OutputName;
}
