#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "Hamiltonian/ParticleOnSphereLaplacianDeltaHamiltonian.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


Complex LaughlinWaveFunction(RealVector& position, int nbrFermions);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEFermionsDeltaOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*SystemGroup) += new SingleStringOption  ('\n', "exact-state", "name of the file containing the vector obtained using exact diagonalization");
  (*SystemGroup) += new BooleanOption  ('\n', "with-timecoherence", "don't use time coherence between two successive evaluation of the wave function");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsDeltaOverlap -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();

  if (((SingleStringOption*) Manager["exact-state"])->GetString() == 0)
    {
      cout << "QHEFermionsDeltaOverlap requires an exact state" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["exact-state"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["exact-state"])->GetString() << endl;
      return -1;      
    }
  FermionOnSphere Space (NbrFermions, 0, LzMax);
  ParticleOnSphereFunctionBasis Basis(LzMax);
  RealVector Location(2 * NbrFermions, true);
  srand48(29457);
/*  for (int k = 0; k < 10; ++k)
    {
      for (int i = 0; i < NbrFermions; ++i)
	{
	  Location[i << 1] = M_PI * drand48();
	  Location[(i << 1) + 1] = 2.0 * M_PI * drand48();
	  cout << Location[i << 1] << " " << Location[(i << 1) + 1] << endl;
	}
      //  Location[4] = Location[0];
      //  Location[5] = Location[1];
      ParticleOnSphereFunctionBasis Basis(LzMax);
      QHEParticleWaveFunctionOperation Operation(&Space, &State, &Location, &Basis);
      Operation.ApplyOperation(Architecture.GetArchitecture());      
      Complex ValueExact (Operation.GetScalar());
      //      Complex ValueExact = Space.EvaluateWaveFunction(State, Location, Basis);
      Complex ValueLaughlin = LaughlinWaveFunction(Location, NbrFermions) * 0.36563112422012;
//      cout << ValueExact  << endl; 
      cout << ValueExact  << " " << ValueLaughlin << " " << (Norm(ValueExact) / Norm(ValueLaughlin)) << endl;        
      cout << "-------------------------------------" << endl;
    }
  return 0;*/
  double Factor = 1.0;
  for (int j = 0; j < NbrFermions; ++j)
    {
      Factor *= 4.0 * M_PI; //2.0 * M_PI * M_PI;
    }
  Complex Overlap;
  Complex ErrorOverlap;
  double Normalization = 0.0;
  double ErrorNormalization = 0.0;
  Complex Tmp;
  Complex Tmp3;
  double Tmp2;
  int NextCoordinates = 0;
  for (int j = 0; j < NbrFermions; ++j)
    {
      Location[j << 1] = acos (1.0- (2.0 * drand48()));
      Location[1 + (j << 1)] = 2.0 * M_PI * drand48();
    }
  for (int i = 0; i < NbrIter; ++i)
    {
      /*      for (int j = 0; j < NbrFermions; ++j)
	{
	  Location[NextCoordinates << 1] = acos (1.0- (2.0 * drand48()));
	  Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * drand48();
	  }*/
      Location[NextCoordinates << 1] = acos (1.0- (2.0 * drand48()));
      Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * drand48();
      NextCoordinates = (int) (((double) NbrFermions) * drand48());
      if (NextCoordinates == NbrFermions)
	--NextCoordinates;
      Tmp = LaughlinWaveFunction (Location, NbrFermions);
      int TimeCoherence = -1;
      if (((BooleanOption*) Manager["with-timecoherence"])->GetBoolean() == true)
	TimeCoherence = NextCoordinates;
      QHEParticleWaveFunctionOperation Operation(&Space, &State, &Location, &Basis, TimeCoherence);
      Operation.ApplyOperation(Architecture.GetArchitecture());      
      Complex ValueExact (Operation.GetScalar());
/*      if ((i == 168131))
	{
	  Space.Debug = 1;
	  cout << i << " " << ValueExact << " " << Tmp << endl;
	  QHEParticleWaveFunctionOperation Operation2(&Space, &State, &Location, &Basis, TimeCoherence);
	  Operation2.ApplyOperation(Architecture.GetArchitecture());      
	  Complex ValueExact (Operation.GetScalar());
	}*/
      Tmp2 = (Tmp.Re * Tmp.Re) + (Tmp.Im * Tmp.Im);
      Tmp3 = (Conj(Tmp) * ValueExact);
      Overlap += Tmp3;// * Factor;
      ErrorOverlap.Re += Tmp3.Re * Tmp3.Re;//  * Factor * Factor;
      ErrorOverlap.Im += Tmp3.Im * Tmp3.Im;//  * Factor * Factor;
//      ErrorOverlap += Tmp3 * Tmp3  * Factor * Factor;
      Normalization += Tmp2;// * Factor;
      ErrorNormalization += Tmp2 * Tmp2;// *  Factor * Factor;
      if ((i > 0) && ((i % 1000) == 0))
	{
	  cout << i << endl;
	  Complex Tmp4 = Overlap / ((double) i);
	  cout << (Tmp4 * Factor);
//	  Complex Tmp5 (sqrt((((ErrorOverlap.Re / ((double) (i))) - (Overlap.Re * Overlap.Re)) / ((double) (i))) / Factor),
//			sqrt((((ErrorOverlap.Im / ((double) (i))) - (Overlap.Im * Overlap.Im)) / ((double) (i))) / Factor));
	  Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
			sqrt( ((ErrorOverlap.Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
	  cout << " +/- " << (Tmp5 * Factor) << endl;
	  double Tmp6 = Normalization / ((double) i);
	  cout << Factor * Tmp6;
//	  Tmp6 = sqrt((((ErrorNormalization / ((double) (i))) - (Tmp6 * Tmp6)) / ((double) (i))) / Factor);
	  double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	  cout << " +/- " << (Tmp7  * Factor) << endl;	  
	  Tmp5.Re /= Tmp4.Re;
	  Tmp5.Im /= Tmp4.Im;
	  Tmp5.Re = fabs(Tmp5.Re);
	  Tmp5.Im = fabs(Tmp5.Im);
	  Tmp5.Re += (Tmp7 / Tmp6);
	  Tmp5.Im += (Tmp7 / Tmp6);
	  Tmp4 *= sqrt(Factor / Tmp6);	  
	  Tmp5.Re *= Tmp4.Re;
	  Tmp5.Im *= Tmp4.Im;
	  cout << Tmp4 << " " << Tmp5 << endl;
	  cout << "-----------------------------------------------" << endl;
	}
    } 
 return 0;
}


Complex LaughlinWaveFunction(RealVector& position, int nbrFermions)
{
  Complex Value (1.0, 0.0);
  Complex Tmp;
  for (int i = 0; i < (nbrFermions - 1); ++i)
    for (int j = i + 1; j < nbrFermions; ++j)
      {
	Tmp.Re = sin(0.5 * (position[j << 1] - position[i << 1])) * cos(0.5 * (position[1 + (i << 1)] - position[1 + (j << 1)]));
	Tmp.Im = sin(0.5 * (position[i << 1] + position[j << 1])) * sin(0.5 * (position[1 + (i << 1)] - position[1 + (j << 1)]));
	Value *= Tmp;
	Value *= Tmp;
	Value *= Tmp;
      }
  return Value;
}
