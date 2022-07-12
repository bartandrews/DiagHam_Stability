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

void calcArekN12_2S22(int*Map, double *Signs, ParticleOnSphereWithSpin* Space, int N=12, int N_phi=22,int Lz_total=0, int Dim=198472577);

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
  (*SystemGroup) += new SingleStringOption ('s', "state", "vector to import in ascii format");
  (*SystemGroup) += new SingleStringOption ('r', "raw-state", "vector to import in FORTRAN binary format");
  (*SystemGroup) += new SingleStringOption ('b', "basis", "description of basis in which vector is formatted (none = standard basis)");
  (*SystemGroup) += new SingleIntegerOption ('c', "coded-basis", "use internally coded basis (1 = A. Wojs 12@22)");
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-state", "use this name for the output vector state instead of standard terminology");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int Lz = Manager.GetInteger("lz-value");
  int Sz = Manager.GetInteger("total-sz");

  ifstream InputFile;
  if (Manager.GetString("state")!=0)
    {
      InputFile.open(Manager.GetString("state"), ios::in);
      if (!InputFile.is_open())
	{
	  cout << "Could not open text-file with vector description!"<<endl;
	  exit(-1);	    
	}
    }
  else if (Manager.GetString("raw-state")!=0)
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
  else if (Manager.GetInteger("coded-basis")!=0)
    {
      switch (Manager.GetInteger("coded-basis"))
	{
	case 1:
	  {
	    cout << "Using fixed basis for N=12, 2S=22"<<endl;
	    if (NbrParticles!=12)
	      {
		cout << "Need to have N=12 particles for this basis"<<endl;
		exit(1);
	      }
	    if (LzMax!=22)
	      {
		cout << "Need to have 2S=22 flux for this basis"<<endl;
		exit(1);
	      }
	    if ((Sz!=0)||(Lz!=0))
	      {
		cout << "Need to have Sz=Lz=0 for this basis"<<endl;
		exit(1);
	      }
	    int Dimension = 198472577;
	    if (Space->GetHilbertSpaceDimension()!=Dimension)
	      {
		cout << "Problem with dimension!"<<endl;
		exit(1);
	      }
	    calcArekN12_2S22(Map, Signs, (ParticleOnSphereWithSpin*)Space, NbrParticles, LzMax, Lz, Dimension);
	    break;
	  }
	default:
	  {
	    cout << "Undefined basis!"<<endl;
	    exit(1);
	    break;
	  }
	}
    }
  else
    {
      for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
	{
	  Map[i]=i;
	  Signs[i]=1.0;
	}
    }

  RealVector ResultingVector(Space->GetHilbertSpaceDimension());

  char *OutputName = Manager.GetString("output-state");

  if (OutputName==NULL)
    {
      OutputName = new char[512];
      sprintf(OutputName,"fermions_sphere_spin_import-%s_n_%d_2S_%d_Sz_%d_lz_%d.vec",
	      Manager.GetString("state"), NbrParticles, LzMax, Sz, Lz);      
    }

  if (Manager.GetString("state")!=0)
    {
      // int VIndex, VLastIndex=-1;
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
      for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
      	ResultingVector[Map[i]]=Signs[i]*InputVector[i];
    }
  else
    {
      cout << "An error occurred"<<endl;
      exit(1);
    }

  ResultingVector.WriteVector(OutputName);

  delete [] OutputName;
}


/**********************************************************************
c PARAMETERS
c     N=12
c     lT=22
c     MT=0
c     numd=198472577
**********************************************************************/

// one-off basis from Arek

void calcArekN12_2S22(int*Map, double *Signs, ParticleOnSphereWithSpin* Space, int N, int N_phi,int Lz_total, int Dim)
{
  int i=0;
  unsigned long StateDesc;
  double Coeff,TmpCoeff;

  for (int k06u=0; k06u<=N_phi; ++k06u)
    for (int k05u=0; k05u<k06u; ++k05u)
      for (int k04u=0; k04u<k05u; ++k04u)
	for (int k03u=0; k03u<k04u; ++k03u)
	  for (int k02u=0; k02u<k03u; ++k02u)
	    for (int k01u=0; k01u<k02u; ++k01u)
	      
	      for (int k06d=0; k06d<=N_phi; ++k06d)
		for (int k05d=0; k05d<k06d; ++k05d)
		  for (int k04d=0; k04d<k05d; ++k04d)
		    for (int k03d=0; k03d<k04d; ++k03d)
		      for (int k02d=0;  k02d<k03d; ++k02d)
			for (int k01d=0; k01d<k02d; ++k01d)
			  {
			    int ksum=k01u+k02u+k03u+k04u+k05u+k06u
			      +k01d+k02d+k03d+k04d+k05d+k06d;

			    if(2*ksum-N*N_phi == Lz_total)
			      {
				Coeff=1.0;
				StateDesc=0x0ul;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k01u, 1, TmpCoeff);
				Coeff*=TmpCoeff;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k02u, 1, TmpCoeff);
				Coeff*=TmpCoeff;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k03u, 1, TmpCoeff);
				Coeff*=TmpCoeff;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k04u, 1, TmpCoeff);
				Coeff*=TmpCoeff;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k05u, 1, TmpCoeff);
				Coeff*=TmpCoeff;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k06u, 1, TmpCoeff);
				Coeff*=TmpCoeff;

				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k01d, 0, TmpCoeff);
				Coeff*=TmpCoeff;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k02d, 0, TmpCoeff);
				Coeff*=TmpCoeff;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k03d, 0, TmpCoeff);
				Coeff*=TmpCoeff;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k04d, 0, TmpCoeff);
				Coeff*=TmpCoeff;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k05d, 0, TmpCoeff);
				Coeff*=TmpCoeff;
				StateDesc = ((ParticleOnSphereWithSpin*)Space)->Ad(StateDesc, k06d, 0, TmpCoeff);
				Coeff*=TmpCoeff;
				if(i >= Dim)
				  {
				    cout << "Error generating Map: size overflow!"<<endl;
				    exit(1);
				  }
				Map[i]=((ParticleOnSphereWithSpin*)Space)->CarefulFindStateIndex(StateDesc,-1);
				if ((Map[i]<0)||(Map[i]>Space->GetHilbertSpaceDimension()))
				  {
				    cout << "Problem with state ["<<i<<"] mapped to "<<Map[i]<<endl;
				    Map[i]=0;
				  }				
				Signs[i]=Coeff;
				++i;
			      }
			  }
    
  if(i != Dim)
    {
      printf("wrong numd");
      exit(1);
    }
}

/************************************************************************/
