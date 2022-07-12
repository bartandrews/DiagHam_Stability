#include "HilbertSpace/BosonOnSquareLatticeWannierSpace.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/ParticleOnSquareLatticeWannierInterface.h"

#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETopInsulatorWannierBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption ('\n', "state", "input file in a Wannier basis representation");
  (*SystemGroup) += new SingleStringOption ('\n', "second-state", "second input file in a Wannier basis representation if ones want to compute the overlap between two FCI states");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ky", "initial momentum along the y direction", 0);

  (*SystemGroup) += new SingleIntegerOption  ('\n', "target-ky", "target momentum on torus along the y direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "no-autodetect", "do not autdetect system parameter from state file name");

  (*SystemGroup) += new BooleanOption  ('\n', "optimize", "optimize overlap over phases");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-iter", "maximum number of iteration for optimization",1000);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tolerance", "tolerance on phases for optimization",0.001);
  (*SystemGroup) += new BooleanOption  ('\n', "random", "random initial value of phases for optimization");

  (*SystemGroup) += new BooleanOption  ('\n', "block-diagonal", "the FQHE eigenstate is already in the torus hilbert space, no need for projection");

  (*SystemGroup) += new SingleStringOption ('\n', "torus", "a vector on the torus, with which to evaluate the overlap");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name to save the resulting vector, instead of replacement by torus conventions");
  (*OutputGroup) += new BooleanOption('\n',"no-save","do not save output vector");
  (*MiscGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (not available for the non-periodic momentum space or the case with spin)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSiteX = Manager.GetInteger("nbr-sitex"); 
  int NbrSiteY = Manager.GetInteger("nbr-sitey"); 
  int TotalKy = Manager.GetInteger("ky");
  int TargetKy = Manager.GetInteger("target-ky");
  int TargetX = TargetKy / NbrSiteY;

  int MaxIter =  Manager.GetInteger("max-iter");
  double Tolerance =  Manager.GetDouble("tolerance");

  if (Manager.GetString("state") == 0)
    {
      cout << "an input state is required: use --state"<<endl;
      return -1;
    }
  if  (Manager.GetBoolean("no-autodetect") == false)
    {
      char * Filename = Manager.GetString("state");
      char * StrWannier = strstr(Filename, "Wannier");
      if (StrWannier == 0)
	StrWannier = strstr(Filename, "wannier");
      if (StrWannier==0)
	{
	  cout << "this does not appear to be a state in a wannier basis - use --no-autodetect to override"<<endl;
	  return -1;
	}
      double Mass = 0.0;
      int TmpKx=0; // not used for Wannier...
      bool Statistics = false;
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(Manager.GetString("state"),
							      NbrParticles, NbrSiteX, NbrSiteY, TmpKx, TotalKy, Mass, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
	  return -1;
	}
    }

  if (Manager.GetString("second-state") != 0)
    {
      cout << "Calculation of overlap between two FCI states"<<endl;
      ComplexVector InputState;
      if (InputState.ReadVector(Manager.GetString("state")) == false)
	{
	  cout << "error while reading " << Manager.GetString("state") << endl;
	  return -1;
	}
      ComplexVector InputStateSecond;
      if (InputStateSecond.ReadVector(Manager.GetString("second-state")) == false)
	{
	  cout << "error while reading " << Manager.GetString("second-state") << endl;
	  return -1;
	}
      cout << "Norm(State1) = " << InputState.Norm()  << endl;
      cout << "Norm(State2) = " << InputStateSecond.Norm()  << endl;
      Complex Ovl = InputState * InputStateSecond;
      cout << "Ovl = " << Ovl << endl;
      cout << "Norm(Ovl) = " << Norm(Ovl) << endl;
      return -1;
    }
 
  // if(Manager.GetBoolean("block-diagonal"))
  //   {
  //     if (Manager.GetString("torus")!=NULL)
  // 	{
  // 	  RealVector ReferenceState;
  // 	  if (ReferenceState.ReadVector(Manager.GetString("torus")) == false)
  // 	    {
  // 	      cout << "error while reading " << Manager.GetString("torus") << endl;
  // 	      return -1;
  // 	    }

  // 	  ComplexVector InputState;
  // 	  if (InputState.ReadVector(Manager.GetString("state")) == false)
  // 	    {
  // 	      cout << "error while reading " << Manager.GetString("state") << endl;
  // 	      return -1;
  // 	    }

  // 	  Complex Ovl=ReferenceState*InputState;
  // 	  cout << "Overlap with torus eigenstate:         "<<Ovl<<endl;
  // 	  cout << "               square overlap:         "<<SqrNorm(Ovl)<<endl;

  // 	  // Write the overlap to a output file
  // 	  ofstream outdata; // outdata is like cin

  // 	  outdata.open("overlap.dat"); // opens the file
  // 	  // if( !outdata ) { // file couldn't be opened
  // 	  // 	cerr << "Error: file could not be opened" << endl;
  // 	  // 	exit(1);
  // 	  // }

  // 	  outdata << Norm(Ovl) << endl;
  // 	  outdata.close();
  // 	}
  //   }
  // else
    // {
      AbstractQHEParticle* Space;
      ParticleOnSquareLatticeWannierInterface *SpacePtr;
  
      if (Manager.GetBoolean("boson") == false)
	{
	  cout << "Currently, no fermionic Wannier states are implemented!"<<endl;
	  return -1;
	}
      else
	{
	  if (Manager.GetInteger("nbr-subbands") == 1)
	    {
	       if(Manager.GetBoolean("block-diagonal"))
		 {
		   Space = new BosonOnSquareLatticeWannierSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKy, TargetX);
		 }
	       else
		 {
		   Space = new BosonOnSquareLatticeWannierSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKy);
		 }
	      // two-step typecast to accommodate double inheritance properties
	      SpacePtr = (ParticleOnSquareLatticeWannierInterface*)((BosonOnSquareLatticeWannierSpace*)Space);
	      if (Manager.GetString("save-hilbert") != 0)
		{
		  Space->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		  return 0;
		}
	    }
	  else
	    {
	      cout << "Wannier states not yet implemented for multiple subbands." << endl;
	      return -1;
	    }
	}

      ParticleOnTorus *TargetSpace = NULL;

      if (Manager.GetBoolean("boson"))
	TargetSpace = new BosonOnTorusShort(NbrParticles, NbrSiteX*NbrSiteY, TargetKy);
      else
	TargetSpace = new FermionOnTorus(NbrParticles, NbrSiteX*NbrSiteY, TargetKy);

      cout << "Dimension of Wannier space: "<<Space->GetHilbertSpaceDimension()<<endl;
      cout << "Dimension of torus space:   "<<TargetSpace->GetHilbertSpaceDimension()<<endl;
  
      int NbrTorusComponents = 0;
      double TorusProjection = 0.0;
      double WeightOtherComponents = 0.0;
      ComplexVector InputState;
      ComplexVector TorusState(TargetSpace->GetHilbertSpaceDimension(), true);
      if (InputState.ReadVector(Manager.GetString("state")) == false)
	{
	  cout << "error while reading " << Manager.GetString("state") << endl;
	  return -1;
	}
      if (Space->GetHilbertSpaceDimension() != InputState.GetVectorDimension())
	{
	  cout << "dimension mismatch between the state (" << InputState.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}

      int Index;
      int FirstTimeFlag=1;
      double PhaseToMakeVectorReal;
      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	{

	  if (SpacePtr->GetLinearizedMomentum(i)==TargetKy)
	    {
	      if(FirstTimeFlag==1 && Manager.GetBoolean("block-diagonal") )
		{
		  PhaseToMakeVectorReal = Arg(InputState[i]);
		  FirstTimeFlag = 0;
		}
	      //cout << InputState[i] << endl;
	      Index = SpacePtr->ProjectToTorus(TargetSpace, i);
	      TorusState[Index] = InputState[i];
	      //if( Manager.GetBoolean("block-diagonal"))  
	      //	TorusState[Index] *= Phase(-PhaseToMakeVectorReal);
	      //cout << "Component inside torus hilbert space : " << 	TorusState[Index] << endl;
	      TorusProjection+=SqrNorm(InputState[i]);
	      ++NbrTorusComponents;
	    }
	  else
	    {
	      WeightOtherComponents+=SqrNorm(InputState[i]);
	      //cout << "SqrNorm of Component outside torus hilbert space : " << SqrNorm(InputState[i]) << endl;
	    }
	  
	}
      if (NbrTorusComponents!=TargetSpace->GetHilbertSpaceDimension())
	{
	  cout << "Remark: Not all torus states found: "<<NbrTorusComponents<<" vs "<<TargetSpace->GetHilbertSpaceDimension()<<endl;
	}

      cout << "Weight of projection onto torus space: "<<TorusProjection<<endl;
      cout << "Weight outside of torus space:         "<<WeightOtherComponents<<endl;

      if (Manager.GetString("torus")!=NULL)
	{
	  if (IsFile(Manager.GetString("torus")) == false)
	    {
	      cout << "Error: Torus vector not found"<<endl;
	      exit(1);
	    }
	  RealVector RealReferenceState;
	  if (RealReferenceState.ReadVectorTest(Manager.GetString("torus")) == false)
	    {
	      ComplexVector ReferenceState;
	      if (ReferenceState.ReadVectorTest(Manager.GetString("torus")) == false)
		{
		  cout << "Could not open torus vector"<<endl;
		  exit(1);
		}
	      if (ReferenceState.ReadVector(Manager.GetString("torus")) == false)
		{
		  cout << "error while reading " << Manager.GetString("torus") << endl;
		  return -1;
		}
	      if (TargetSpace->GetHilbertSpaceDimension() != ReferenceState.GetVectorDimension() && !(Manager.GetBoolean("block-diagonal")))
		{
		  cout << "dimension mismatch between the state (" << ReferenceState.GetVectorDimension() << ") and the Hilbert space (" << TargetSpace->GetHilbertSpaceDimension() << ")" << endl;
		  return -1;
		}
	      Complex Ovl=ReferenceState*TorusState;
	      cout << "Overlap with torus eigenstate:         "<<Ovl<<endl;
	      cout << "               square overlap:         "<<SqrNorm(Ovl)<<endl;
	      
	      // Write the overlap to a output file
	      ofstream outdata; // outdata is like cin
	      
	      outdata.open("overlap.dat"); // opens the file
	      // if( !outdata ) { // file couldn't be opened
	      // 	cerr << "Error: file could not be opened" << endl;
	      // 	exit(1);
	      // }
	      
	      outdata << SqrNorm(Ovl) << endl;
	      outdata << TorusProjection << endl;
	      outdata.close();

	    }
	  else
	    {
	      RealVector ReferenceState;
	      if (ReferenceState.ReadVector(Manager.GetString("torus")) == false)
		{
		  cout << "error while reading " << Manager.GetString("torus") << endl;
		  return -1;
		}
	      if (TargetSpace->GetHilbertSpaceDimension() != ReferenceState.GetVectorDimension() && !(Manager.GetBoolean("block-diagonal")))
		{
		  cout << "dimension mismatch between the state (" << ReferenceState.GetVectorDimension() << ") and the Hilbert space (" << TargetSpace->GetHilbertSpaceDimension() << ")" << endl;
		  return -1;
		}
	      Complex Ovl=ReferenceState*TorusState;
	      cout << "Overlap with torus eigenstate:         "<<Ovl<<endl;
	      cout << "               square overlap:         "<<SqrNorm(Ovl)<<endl;
	      
	      // Write the overlap to a output file
	      ofstream outdata; // outdata is like cin
	      
	      outdata.open("overlap.dat"); // opens the file
	      // if( !outdata ) { // file couldn't be opened
	      // 	cerr << "Error: file could not be opened" << endl;
	      // 	exit(1);
	      // }
	      
	      outdata << SqrNorm(Ovl) << endl;
	      outdata << TorusProjection << endl;
	      outdata.close();
	    }

	}

      if (!Manager.GetBoolean("no-save"))
	{
	  char *Output = Manager.GetString("output-file");
	  if (Output==NULL)
	    {      
	      Output = AddSegmentInFileName(Manager.GetString("state"), "to_torus_", "annier_", false);
	      if (Output == NULL)
		Output = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "torus_vec");
	      if (Output == NULL)
		Output = AddExtensionToFileName(Manager.GetString("state"), "torus_vec");
	    }
	  TorusState.WriteVector(Output);
	}
  
      //delete Output;
      delete Space;
      delete TargetSpace;
    
  
  return 0;
}

