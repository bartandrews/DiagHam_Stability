#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "Options/Options.h"


#include "Vector/ComplexVector.h"

#include "GeneralTools/ConfigurationParser.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/Endian.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


void RandomUV (ComplexVector& uv, RealVector& positions, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator);

void RandomUVOneCoordinate(ComplexVector& uv, RealVector& positions, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);

void FlipCoordinates (ComplexVector& uv, int i, int j);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereSUKToU1MCOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
   Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

	(*SystemGroup) += new SingleStringOption  ('\0', "initial-vector", "name of the file corresponding to the ground state of the whole system");
	(*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 0);
	(*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 0);
	(*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
   (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereSUKToU1MCOverlap -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

	int NbrParticles = Manager.GetInteger("nbr-particles");
	int LzMax = Manager.GetInteger("lzmax");
	int TotalLz = Manager.GetInteger("total-lz");
  bool HaldaneBasisFlag = ((BooleanOption*) Manager["haldane"])->GetBoolean();
	bool FermionFlag = true;
	char* StateFileName = Manager.GetString("initial-vector");
  AbstractQHEParticle* ExactSpace = 0;
  RealVector* ExactState = 0;
  AbstractFunctionBasis* ExactBasis = 0;


	if (StateFileName == 0)
	{
		cout << "error, an initial state file should be provided. See man page for option syntax or type  FQHESpheretimeEvolution -h" << endl;
		return -1;
	}
	
	
	if (FQHEOnSphereFindSystemInfoFromVectorFileName(StateFileName, NbrParticles,LzMax, TotalLz, FermionFlag) == false)
	{
		return -1;
	}
	
	if (FermionFlag == false)
	{
		if(HaldaneBasisFlag == false)
		{
			#ifdef  __64_BITS__
			if((LzMax + NbrParticles - 1) < 63)
				#else
				if ((LzMax + NbrParticles - 1) < 31)
					#endif
					{	
						ExactSpace = new BosonOnSphereShort (NbrParticles, TotalLz, LzMax);
					}
					else
					{
						cout << "This Space requires BosonOnSphere class for which this kind of calculation is not available"<<endl;
						return -1;
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
			if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
			{
				cout << "NbrParticles is not defined or as a wrong value" << endl;
				return -1;
			}
			if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax < 0))
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
			if (MaxNbrLz != (LzMax + 1))
			{
				cout << "wrong LzMax value in ReferenceState" << endl;
				return -1;
			}
			#ifdef  __64_BITS__
			if (LzMax  < 63)
				#else
				if (LzMax  < 31)
					#endif
					ExactSpace = new BosonOnSphereHaldaneBasisShort (NbrParticles,TotalLz, LzMax, ReferenceState);
		}
	}
	else
	{
		if (HaldaneBasisFlag == false)
		{
			#ifdef __64_BITS__
			if (LzMax <= 63)
				#else
				if (LzMax <= 31)
					#endif
					{
						ExactSpace = new FermionOnSphere (NbrParticles,TotalLz,LzMax);
					}
					else
					{
						cout << "This Space requires FermionOnSphereLong class for which this kind of calculation is not available"<<endl;
						return -1;
					}
		}
		else
		{
			int * ReferenceState = 0;
			ConfigurationParser ReferenceStateDefinition;
			if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
			{
				ReferenceStateDefinition.DumpErrors(cout) << endl;
				return -1;
			}
			if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
			{
				cout << "NbrParticles is not defined or as a wrong value" << endl;
				return -1;
			}
			if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax < 0))
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
			if (MaxNbrLz != (LzMax + 1))
			{
				cout << "wrong LzMax value in ReferenceState" << endl;
				return -1;
			}
			ExactSpace = new FermionOnSphereHaldaneBasis (NbrParticles,TotalLz, LzMax, ReferenceState);
		}
	}
	
	ExactState = new RealVector [1];
	if (ExactState[0].ReadVector (StateFileName) == false)
	{
		cout << "can't open vector file " << StateFileName<< endl;
	  return -1;      
	}


  if (ExactSpace->GetHilbertSpaceDimension() != ExactState[0].GetVectorDimension())
    {
      cout << "dimension mismatch : hilbert space = " << ExactSpace->GetHilbertSpaceDimension() << ", exact state = " << ExactState[0].GetVectorDimension() << endl;
      return -1;
    }
    
  ExactBasis = new ParticleOnSphereFunctionBasis (LzMax);	
  
  AbstractRandomNumberGenerator* RandomNumber = 0;
  RandomNumber = new StdlibRandomNumberGenerator (29457);
  
	
  
  ComplexVector UV (NbrParticles * 2, true);
  RealVector TmpPositions (NbrParticles * 2, true);
  RandomUV (UV, TmpPositions, NbrParticles, RandomNumber);
  QHEParticleWaveFunctionOperation Operation(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
  Operation.ApplyOperation(Architecture.GetArchitecture());      
  Complex Tmp = Operation.GetScalar();
  cout << "wave function initial value : " << Tmp << endl;;
  cout << "test vanishing properties : " << endl;
  for (int i = 1; i < 7; ++i)
    {
      TmpPositions[(i << 1)] = TmpPositions[0];// * (1 + 0.0001 *RandomNumber->GetRealRandomNumber());
      TmpPositions[(i << 1) + 1] = TmpPositions[1];// * (1 + 0.0001 *RandomNumber->GetRealRandomNumber());
      QHEParticleWaveFunctionOperation Operation2(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
      Operation2.ApplyOperation(Architecture.GetArchitecture());      
      Tmp = Operation2.GetScalar();
      cout << (i  + 1) << " body cancellation : " << Norm(Tmp) << endl;
    }
  RandomUV (UV, TmpPositions, NbrParticles, RandomNumber);
	QHEParticleWaveFunctionOperation Operation3(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
	Operation3.ApplyOperation(Architecture.GetArchitecture());      
	Complex Tmp1 = Operation3.GetScalar();
	cout << "wave function initial value : " << Tmp1 << endl;;
	cout << "test vanishing properties with 2 groups: " << endl;
	for (int i = 1; i < NbrParticles; ++i)
	{
		int TmpI=abs((i+1)/2);
		if(i%2==1)
		{
			TmpPositions[TmpI << 1] = TmpPositions[0];// * (1 + 0.0001 *RandomNumber->GetRealRandomNumber());
			TmpPositions[(TmpI << 1) + 1] = TmpPositions[1];// * (1 + 0.0001 *RandomNumber->GetRealRandomNumber());
			QHEParticleWaveFunctionOperation Operation2(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
			Operation2.ApplyOperation(Architecture.GetArchitecture());      
			Tmp1 = Operation2.GetScalar();
			cout <<"first cluster of "<< (TmpI  + 1) << " particles cancellation : " << Norm(Tmp1) << endl;
		}
		else
		{
			TmpPositions[2*(NbrParticles  - (TmpI << 1)) ] = TmpPositions[2*(NbrParticles-1)];// * (1 + 0.0001 *RandomNumber->GetRealRandomNumber());
			TmpPositions[2*(NbrParticles  - (TmpI << 1)) + 1] = TmpPositions[2*(NbrParticles-1)+1];// * (1 + 0.0001 *RandomNumber->GetRealRandomNumber());
			QHEParticleWaveFunctionOperation Operation2(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
			Operation2.ApplyOperation(Architecture.GetArchitecture());      
			Tmp1 = Operation2.GetScalar();
			cout <<"second cluster of " << (TmpI  + 1) << " particles cancellation : " << Norm(Tmp1) << endl;
		}
	}
	
	RandomUV (UV, TmpPositions, NbrParticles, RandomNumber);
	QHEParticleWaveFunctionOperation Operation4(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
	Operation4.ApplyOperation(Architecture.GetArchitecture());      
	Complex Tmp2 = Operation3.GetScalar();
	cout << "wave function initial value : " << Tmp1 << endl;;
	cout << "test vanishing properties with 3 groups: " << endl;
	
	for (int i = 3; i < NbrParticles; i++)
	{
		int TmpI=abs(i/3);
		if(i%3 == 1)
		{
			TmpPositions[2*i] = TmpPositions[2];
			TmpPositions[2*i + 1] = TmpPositions[3];
			QHEParticleWaveFunctionOperation Operation2(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
			Operation2.ApplyOperation(Architecture.GetArchitecture());      
			Tmp1 = Operation2.GetScalar();
			cout <<"second cluster of "<< (TmpI  + 1) << " particles cancellation : " << Norm(Tmp1) << endl;
		}
		else
		{
			if (i%3 == 2)
			{
			TmpPositions[2*i] = TmpPositions[4];
			TmpPositions[2*i + 1] = TmpPositions[5];
			QHEParticleWaveFunctionOperation Operation2(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
			Operation2.ApplyOperation(Architecture.GetArchitecture());      
			Tmp1 = Operation2.GetScalar();
			cout <<"third cluster of " << (TmpI  + 1) << " particles cancellation : " << Norm(Tmp1) << endl;
			}
			else
			{
				TmpPositions[2*i] = TmpPositions[0];
				TmpPositions[2*i + 1] = TmpPositions[1];
				QHEParticleWaveFunctionOperation Operation2(ExactSpace, &(ExactState[0]), &TmpPositions, ExactBasis);
				Operation2.ApplyOperation(Architecture.GetArchitecture());      
				Tmp1 = Operation2.GetScalar();
				cout <<"first cluster of " << (TmpI  + 1) << " particles cancellation : " << Norm(Tmp1) << endl;
			}			
		}
	}
	
		
	
  return 0;
}

void RandomUV (ComplexVector& uv, RealVector& positions, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  for (int j = 0; j < nbrParticles; ++j)
    {
      double x = acos (1.0 - (2.0 * randomNumberGenerator->GetRealRandomNumber()));
      double y = 2.0 * M_PI * randomNumberGenerator->GetRealRandomNumber();
      positions[2 * j] = x;
      positions[(2 * j) + 1] = y;
      uv.Re(2 * j) = cos(0.5 * x);
      uv.Im(2 * j) = uv.Re(2 * j) * sin(0.5 * y);
      uv.Re(2 * j) *= cos(0.5 * y);
      uv.Re(2 * j + 1) = sin(0.5 * x);
      uv.Im(2 * j + 1) = - uv.Re(2 * j + 1) * sin(0.5 * y);
      uv.Re(2 * j + 1) *= cos(0.5 * y);      
    }
}

void RandomUVOneCoordinate(ComplexVector& uv, RealVector& positions, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  double x = acos (1.0 - (2.0 * randomNumberGenerator->GetRealRandomNumber()));
  double y = 2.0 * M_PI * randomNumberGenerator->GetRealRandomNumber();
  positions[coordinate] = x;
  uv.Re(coordinate) = cos(0.5 * x);
  uv.Im(coordinate) = uv.Re(coordinate) * sin(0.5 * y);
  uv.Re(coordinate) *= cos(0.5 * y);
  ++coordinate;
  positions[coordinate] = y;
  uv.Re(coordinate) = sin(0.5 * x);
  uv.Im(coordinate) = - uv.Re(coordinate) * sin(0.5 * y);
  uv.Re(coordinate) *= cos(0.5 * y);      
}


void FlipCoordinates (ComplexVector& uv, int i, int j)
{
  Complex Tmp = uv[2 * i];
  uv.Re(2 * i) = uv.Re(2 * j);
  uv.Re(2 * j) = Tmp.Re;
  uv.Im(2 * i) = uv.Im(2 * j);
  uv.Im(2 * j) = Tmp.Im;
  Tmp = uv[2 * i + 1];
  uv.Re(2 * i + 1) = uv.Re(2 * j + 1);
  uv.Re(2 * j + 1) = Tmp.Re;
  uv.Im(2 * i + 1) = uv.Im(2 * j + 1);
  uv.Im(2 * j + 1) = Tmp.Im;
}
