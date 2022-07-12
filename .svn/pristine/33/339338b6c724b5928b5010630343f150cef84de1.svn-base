#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphere.h"

#include "Operator/ParticleOnSphereNBodyOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
  OptionManager Manager ("QHEBosonsNBodyCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (overriding the one found in the vector file name)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (overriding the one found in the vector file name)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the lz value corresponding to the eigenvector (overriding the one found in the vector file name)", -1);
  (*SystemGroup) += new SingleStringOption  ('s', "state", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-nbody", "number of n-body term (aka n of the <Psi+^n(z) Psi^n(z)>", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "integrate", "integrate ovez the z value of the mean value instead taking z=0", false);
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsNBodyCorrelation -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["lz-value"])->GetInteger();
  int NbrNBody = ((SingleIntegerOption*) Manager["nbr-nbody"])->GetInteger();
  if (((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "QHEBosonsCorrelation requires a state" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["state"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
      return -1;      
    }
  char* StrNbrParticles = 0;
  if (NbrParticles < 0)
    {
      StrNbrParticles = strstr(((SingleStringOption*) Manager["state"])->GetString(), "_n_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 3;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      NbrParticles = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess number of particles from file name " << ((SingleStringOption*) Manager["state"])->GetString() << endl
	       << "use --nbr-particles option" << endl;
	  return -1;            
	}
    }
  if (LzMax < 0)
    {
      StrNbrParticles = strstr(StrNbrParticles, "_2s_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 4;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '_') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '_') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      LzMax = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '_';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess maximum momentum from file name " << ((SingleStringOption*) Manager["state"])->GetString() << endl
	       << "use --lzmax option" << endl;
	  return -1;            
	}
    }
  if (Lz < 0)
    {
      StrNbrParticles = strstr(StrNbrParticles, "_lz_");
      if (StrNbrParticles != 0)
	{
	  StrNbrParticles += 4;
	  int SizeString = 0;
	  while ((StrNbrParticles[SizeString] != '\0') && (StrNbrParticles[SizeString] != '.') && (StrNbrParticles[SizeString] >= '0') 
		 && (StrNbrParticles[SizeString] <= '9'))
	    ++SizeString;
	  if ((StrNbrParticles[SizeString] == '.') && (SizeString != 0))
	    {
	      StrNbrParticles[SizeString] = '\0';
	      Lz = atoi(StrNbrParticles);
	      StrNbrParticles[SizeString] = '.';
	      StrNbrParticles += SizeString;
	    }
	  else
	    StrNbrParticles = 0;
	}
      if (StrNbrParticles == 0)
	{
	  cout << "can't guess total momentum from file name " << ((SingleStringOption*) Manager["state"])->GetString() << endl
	       << "use --lz-value option" << endl;
	  return -1;            
	}
    }


  BosonOnSphere Space (NbrParticles, Lz, LzMax);
  ParticleOnSphereFunctionBasis Basis(LzMax);
  int* AnnihilationIndices = new int [NbrNBody];
  int* CreationIndices = new int [NbrNBody];
  if (((BooleanOption*) Manager["integrate"])->GetBoolean() == false)
    {
      double Factor = 1.0;
      double Tmp = ((double) (LzMax + 1)) / (4.0 * M_PI);
      for (int i = 0; i < NbrNBody; ++i)
	{
	  AnnihilationIndices[i] = 0;
	  CreationIndices[i] = 0;
	  Factor *= Tmp;
	}
      ParticleOnSphereNBodyOperator Operator (&Space, CreationIndices, AnnihilationIndices, NbrNBody);
      Factor *= Operator.MatrixElement(State, State).Re;
      cout << Factor << endl;
    }
  else
    {
    }
  delete[] AnnihilationIndices;
  delete[] CreationIndices;
  return 0;
}


