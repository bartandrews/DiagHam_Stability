#include "HilbertSpace/BosonOnTorusWithSpin.h"
#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "Options/Options.h"
#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <cstdlib>
#include <climits>
#include <cmath>
#include <cstring>
#include <sys/time.h>
#include <stdio.h>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;



int main(int argc, char** argv)
{
  cout.precision(14);
  
  // some running options and help
  OptionManager Manager ("FQHETorusTwoLandauLevelProjection" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* SystemGroup = new OptionGroup("system options");
  
  Manager += MiscGroup;
  Manager += SystemGroup;
  Manager += OutputGroup;
  
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the file corresponding to the state to be projected");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particule (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption ('l', "max-momentum", "maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption ('y', "ky-momentum", "constraint on the total momentum modulo the maximum momentum (negative if none)", 0);
  (*SystemGroup) += new BooleanOption  ('u', "unnormalized", "leave the vector unormalized at the end");
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the result of the projection into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the result of the projection into a text file");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusTwoLandauLevelProjection -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrFluxQuanta = Manager.GetInteger("max-momentum"); 
  int Ky = Manager.GetInteger("ky-momentum") % NbrFluxQuanta;
  
  bool FermionFlag = false;
  char * StateName= Manager.GetString("state");
  char* OutputFileName = ((SingleStringOption*) Manager["bin-output"])->GetString();
  char* OutputTxtFileName = ((SingleStringOption*) Manager["txt-output"])->GetString();
  bool UnNormalize = Manager.GetBoolean("unnormalized"); 
  bool ComplexFlag;  
  
  if (StateName == 0)
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereTwoLandauLevelProjection -h" << endl;
      return -1;
    }
  
 
  if (IsFile(StateName) == false)
    {
      cout << "state " << StateName << " does not exist or can't be opened" << endl;
      return -1;
    }
  
 if (FQHEOnTorusFindSystemInfoFromVectorFileName(StateName,NbrParticles, NbrFluxQuanta, Ky, FermionFlag) == false)
    {
      cout << "error while retrieving system parameters from file name " <<StateName  << endl;
      return -1;
    } 	
 cout <<"Using system parameters: N = "<<NbrParticles <<" Kmax = " <<NbrFluxQuanta<< " Ky = " <<Ky<<endl;	      
 RealVector State;
 ComplexVector ComplexState;
 if (State.ReadVectorTest(StateName) == true)
   {
     if (State.ReadVector ( StateName) == false)
       {
	 cout << "can't open vector file " << StateName << endl;
	 return -1;      
       }
     ComplexFlag = false;
   }   
 else
   {
     if (ComplexState.ReadVector (StateName) == false)
       {
	 cout << "can't open vector file " << StateName << endl;
	 return -1;
       }
     ComplexFlag = true;
   }
  
  if (FermionFlag) 
    {
      
      /*ParticleOnTorusWithSpin* Space = new FermionOnSphereTwoLandauLevels (NbrParticles, TotalLz, LzMaxUp, LzMaxDown);
      
      if (Space->GetHilbertSpaceDimension() != GroundState.GetVectorDimension())
	{
	  cout <<Space->GetHilbertSpaceDimension()<<" "<<GroundState.GetVectorDimension()<<endl;
	  cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
	  return -1;
	}
      
      FermionOnSphere * FinalSpace=new FermionOnSphere(NbrParticles,TotalLz,LzMaxDown);
      
      RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
      
      ((FermionOnSphereTwoLandauLevels*) Space)->ProjectionInTheLowestLevel(GroundState,OutputVector,FinalSpace);
      
      if ( !UnNormalize ) 
	OutputVector.Normalize();
      
      OutputVector.WriteVector(OutputFileName);	
      ofstream File;
      if(OutputTxtFileName!=0)
	{
	  File.open(OutputTxtFileName, ios::binary | ios::out);
	  File.precision(14);
	  for (long i = 0; i < FinalSpace->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << OutputVector[i] << " ";
	      FinalSpace->PrintStateMonomial(File, i) << endl;
	    }
	}
      File.close();
      */
    }
  else 
    {
      BosonOnTorusWithSpin* Space = new BosonOnTorusWithSpin(NbrParticles, NbrFluxQuanta, Ky);
      
      BosonOnTorusShort* FinalSpace = 0;

#ifdef  __64_BITS__
      if ( (NbrFluxQuanta + NbrParticles - 1) < 63)
#else
	if ((NbrFluxQuanta + NbrParticles - 1) < 31)
#endif
	  {
	    FinalSpace = new BosonOnTorusShort(NbrParticles, NbrFluxQuanta, Ky);    
	  }
	else
	  {
	    cout <<"This Hilbert space is not supported yet"<<endl;
	    return -1;
	    //	    FinalSpace = new BosonOnTorus(NbrParticles, MaxMomentum, Momentum);
	  }
      
      if(ComplexFlag == false)
	{
	  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
	    {
	      cout <<Space->GetHilbertSpaceDimension()<<" "<<State.GetVectorDimension()<<endl;
	      cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
	      return -1;
	    } 
	  
	  RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
	  
	  Space->ProjectionInTheLowestLevel(State,OutputVector,FinalSpace);
	  
	  if ( !UnNormalize ) 
	    OutputVector.Normalize();
	  
	  OutputVector.WriteVector(OutputFileName);	
	  ofstream File;
	  if(OutputTxtFileName!=0)
	    {
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);
	      for (long i = 0; i < FinalSpace->GetLargeHilbertSpaceDimension(); ++i)
		{
		  File << OutputVector[i] << " ";
		  FinalSpace->PrintStateMonomial(File, i) << endl;
		}
	    }
	  File.close();
	}
      else
	{
	  if (Space->GetHilbertSpaceDimension() != ComplexState.GetVectorDimension())
	    {
	      cout <<Space->GetHilbertSpaceDimension()<<" "<<ComplexState.GetVectorDimension()<<endl;
	      cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
	      return -1;
	    } 
	  
	  ComplexVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
	  
	  Space->ProjectionInTheLowestLevel(ComplexState,OutputVector,FinalSpace);
	  
	  if ( !UnNormalize ) 
	    OutputVector.Normalize();
	  
	  OutputVector.WriteVector(OutputFileName);	
	  ofstream File;
	  if(OutputTxtFileName!=0)
	    {
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);
	      for (long i = 0; i < FinalSpace->GetLargeHilbertSpaceDimension(); ++i)
		{
		  File << OutputVector[i] << " ";
		  FinalSpace->PrintStateMonomial(File, i) << endl;
		}
	    }
	  File.close(); 
	}
    }
  return 0;
}
