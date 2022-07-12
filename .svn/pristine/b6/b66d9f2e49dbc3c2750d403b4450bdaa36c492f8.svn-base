#include "Vector/RealVector.h"

#include "HilbertSpace/ParticleOnSphereManager.h"

#include "Operator/ParticleOnSphereWithSpinDensityDensityOperator.h"
#include "Operator/ParticleOnSphereWithSpinDensityOperator.h"
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

#include "GeneralTools/FilenameTools.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
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
  OptionManager Manager ("FQHESphereFermionsWithSpinCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  ParticleOnSphereManager ParticleManager(true, false, 2);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the lz value corresponding to the eigenvector (override autodetection from input file name if greater or equal to zero)", 0, true, 0);
  (*SystemGroup) += new SingleStringOption  ('e', "eigenstate", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleStringOption  ('i', "interaction-name", "name of the interaction (used for output file name)", "sphere_spin");
  (*SystemGroup) += new SingleStringOption ('a', "add-filename", "add a string with additional informations to the output file name(just before the .dat extension)");
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one or two body coefficients that are requested to evaluate the density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('b', "bilayer", "adjust normalization as a bilayer correlation function");
  (*SystemGroup) += new BooleanOption  ('\n', "shift", "calculate 'shift' as defined for bilayer states");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpinCorrelation -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int TotalLz = ((SingleIntegerOption*) Manager["total-lz"])->GetInteger();
  int SzTotal = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  int NbrPoints = ((SingleIntegerOption*) Manager["nbr-points"])->GetInteger();
  bool DensityFlag = Manager.GetBoolean("density");
  bool Statistics = true;

  if (Manager.GetString("eigenstate") == 0)
    {
      cout << "FQHESphereFermionsWithSpinCorrelation requires a state" << endl;
      return -1;
    }
  if (NbrParticles==0)
    if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("eigenstate"), NbrParticles, LzMax, TotalLz, 
							     SzTotal, Statistics) == false)
      {
	cout << "error while retrieving system informations from file name " << Manager.GetString("eigenstate") << endl;
	return -1;
      }
  cout << NbrParticles << " " << TotalLz << " " << SzTotal << " " << endl;

  RealVector State;
  if (State.ReadVector (Manager.GetString("eigenstate")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
      return -1;      
    }
  //  char* OutputNameCorr = new char [256 + strlen (((SingleStringOption*) Manager["interaction-name"])->GetString())];
  //   if (((SingleStringOption*) Manager["add-filename"])->GetString() == 0)
  //     {
  //       sprintf(OutputNameCorr,"%s",Manager.GetString("eigenstate"));
  //       OutputNameCorr[strlen(OutputNameCorr)-4]='\0';
  //       sprintf(OutputNameCorr,"%s.%s.dat",OutputNameCorr,typeStr);
  //     }
  //   else
  //     {
  //       sprintf (OutputNameCorr, "fermions_%s_n_%d_2s_%d_sz_%d_lz_%d_%s.%s.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax, SzTotal, TotalLz, 
  // 	       ((SingleStringOption*) Manager["add-filename"])->GetString(), typeStr);
  //     }

  ParticleOnSphereWithSpin* Space;
#ifdef __64_BITS__
  if (LzMax <= 31)
    {
      Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, SzTotal, 0);
    }
  else
    {
      cout << "States of this Hilbert space cannot be represented in a single word." << endl;
      return -1;
    }
#else
  if (LzMax <= 15)
    {
      Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, SzTotal, 0);
    }
  else
    {
      cout << "States of this Hilbert space cannot be represented in a single word." << endl;
      return -1;
    }
#endif

  cout << "dim = " << Space->GetHilbertSpaceDimension() << endl;

  ParticleOnSphereFunctionBasis Basis(LzMax);

  Complex* Sum = new Complex [3];
  Complex Sum2 (0.0, 0.0);
  Complex TmpValue;
  RealVector Value(2, true);
  double X = 0.0;
  double XInc = M_PI / ((double) NbrPoints);
  Complex** PrecalculatedValues = new Complex* [3];
  for (int i = 0; i < 3; ++i)
    PrecalculatedValues[i] = new Complex [LzMax + 1];
   
  if (DensityFlag == true)
    {
      for (int i = 0; i <= LzMax; ++i)
	{
	  ParticleOnSphereWithSpinDensityOperator Operator (Space, i, 0, i, 0);
	  PrecalculatedValues[0][i] =   Operator.MatrixElement(State, State);
	}
      for (int i = 0; i <= LzMax; ++i)
	{
	  ParticleOnSphereWithSpinDensityOperator Operator (Space, i, 1, i, 1);
	  PrecalculatedValues[1][i] =   Operator.MatrixElement(State, State);
	}
    }
  else
    {
      Basis.GetFunctionValue(Value, TmpValue, LzMax);
      for (int i = 0; i <= LzMax; ++i)
	{
	  ParticleOnSphereWithSpinDensityDensityOperator Operator (Space, i, 0, LzMax, 0, i, 0, LzMax, 0);
	  PrecalculatedValues[0][i] =   Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
	}
      for (int i = 0; i <= LzMax; ++i)
	{
	  ParticleOnSphereWithSpinDensityDensityOperator Operator (Space, i, 1, LzMax, 1, i, 1, LzMax, 1);
	  PrecalculatedValues[2][i] =  Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
	}
      for (int i = 0; i <= LzMax; ++i)
	{
	  ParticleOnSphereWithSpinDensityDensityOperator Operator1 (Space, i, 0, LzMax, 1, i, 0, LzMax, 1);
	  PrecalculatedValues[1][i] =  Operator1.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
	}
    }

  ofstream File;
  File.precision(14);
  if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
    File.open(((SingleStringOption*) Manager["output-file"])->GetString(), ios::binary | ios::out);
  else
    {
      char* TmpFileName = 0;
      if (DensityFlag == false)
	{
	  if (Manager.GetBoolean("coefficients-only"))
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rhorho-c");
	  else
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rhorho");
	}
      else
	{
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rho");
	}
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("eigenstate") << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  
  double Factor2 = sqrt (0.5 * LzMax );
  if (((BooleanOption*) Manager["radians"])->GetBoolean() == true) 
    Factor2 = 1.0;
  if (DensityFlag == true)
    {
      double Factor1 = 1.0;
      File << "# density coefficients for " << Manager.GetString("eigenstate") << endl;
      File << "#" << endl << "# (l+S) n_l^{u} n_l^{d}" << endl;
      double Sum2 = 0.0;
      for (int i = 0; i <= LzMax; ++i)
	{
	  File << "# " << i;
	  for (int j = 0; j < 2; ++j)
	    {
	      File << " " << PrecalculatedValues[j][i].Re;
	      Sum2 += PrecalculatedValues[j][i].Re;
	    }
	  File << endl;
	}
      File << "# sum = " << Sum2 << endl;
      File << "# dist (rad) rho_{u} rho_{d} rho_{u}+rho_{d}" << endl;
      for (int x = 0; x < NbrPoints; ++x)
	{
	  Value[0] = X;
	  for (int j = 0; j < 3; ++j)
	    Sum[j] = 0.0;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      Basis.GetFunctionValue(Value, TmpValue, i);
	      Sum[0] += PrecalculatedValues[0][i] * (Conj(TmpValue) * TmpValue);
	      Sum[1] += PrecalculatedValues[1][i] * (Conj(TmpValue) * TmpValue);
	    }
	  File << (X * Factor2);
	  File << " " << (Sum[0].Re  * Factor1) << " " << (Sum[1].Re  * Factor1) << " " << ((Sum[0].Re + Sum[1].Re) * Factor1) << endl;
	  X += XInc;
	}
      File.close();     
      return 0;
    }
  if (Manager.GetBoolean("coefficients-only"))
    {
      File << "# density-density correlation coefficients for " << Manager.GetString("eigenstate") << endl;
      File << "#" << endl << "# (l+S) n_l^{u,u} n_l^{u,d} n_l^{d,d}" << endl;
      for (int i = 0; i <= LzMax; ++i)
	{
	  File << i;
	  for (int j=0; j<3; ++j)
	    File << " " << PrecalculatedValues[j][i].Re;
	  File << endl;
	}
      File.close();
      for (int i = 0; i < 3; ++i)
	delete[] PrecalculatedValues[i];
      delete[] PrecalculatedValues;      
    }
  else
    {
      double Factor1;
      if (Manager.GetBoolean("bilayer"))
	Factor1 = (64.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
      else
	Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
      File << "# dist (rad) rho_{u,u} rho_{u,d} rho_{d,d}";
      for (int x = 0; x < NbrPoints; ++x)
	{
	  Value[0] = X;
	  for (int j = 0; j < 3; ++j)
	    Sum[j] = 0.0;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      Basis.GetFunctionValue(Value, TmpValue, i);
	      for (int j = 0; j < 3; ++j)	    
		Sum[j] += PrecalculatedValues[j][i] * (Conj(TmpValue) * TmpValue);
	    }
	  File << (X * Factor2);
	  for (int j = 0; j < 3; ++j)
	    File << " " << (Norm(Sum[j])  * Factor1);
	  File << endl;
	  X += XInc;
	}
      File.close();
      for (int i = 0; i < 3; ++i)
	delete[] PrecalculatedValues[i];
      delete[] PrecalculatedValues;
      delete[] Sum;

      // CALCULATION OF "SHIFT OPERATOR":
      if (Manager.GetBoolean("shift"))
	{
	  char* OutputNameCorr = new char [256 + strlen (((SingleStringOption*) Manager["interaction-name"])->GetString())];
	  if (((SingleStringOption*) Manager["add-filename"])->GetString() == 0)
	    {
	      sprintf(OutputNameCorr,"%s",Manager.GetString("eigenstate"));
	      OutputNameCorr[strlen(OutputNameCorr)-4]='\0';
	      sprintf(OutputNameCorr,"%s.shift.dat",OutputNameCorr);
	    }
	  else
	    {
	      sprintf (OutputNameCorr, "fermions_%s_n_%d_2s_%d_Sz_%d_lz_%d_%s.shift.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrParticles, LzMax, SzTotal, TotalLz, 
		       ((SingleStringOption*) Manager["add-filename"])->GetString());
	    }

	  Complex* PrecalculatedValues2 = new Complex [2*(LzMax + 1)];
	  Complex TmpValue2;      
	  int Pos=0;      
	  for (int m = 0; m <= LzMax; ++m)
	    {	    
	      ParticleOnSphereWithSpinDensityDensityOperator Operator (Space, m, 0, LzMax, 1, m, 0, LzMax, 1);
	      PrecalculatedValues2[Pos] = Operator.MatrixElement(State, State);
	      ++Pos;
	      ParticleOnSphereWithSpinDensityDensityOperator Operator2 (Space, LzMax, 0, m, 1, m, 1, LzMax, 0);
	      PrecalculatedValues2[Pos] = Operator2.MatrixElement(State, State);
	      ++Pos;
	    }

	  ofstream File;
	  File.precision(14);
	  File.open(OutputNameCorr, ios::binary | ios::out);
	  X=0.0;
	  for (int x = 0; x < NbrPoints; ++x)
	    {
	      Value[0] = X;
	      Value[1] = 0.0;
	      Pos = 0;
	      Complex Sum3 = 0.0;
	      Sum2 = 0.0;	  
	      for (int m = 0; m <= LzMax; ++m)
		{
		  Basis.GetFunctionValue(Value, TmpValue, m);
		  double CommonFactor=SqrNorm(TmpValue);      
		  Sum3 += CommonFactor*PrecalculatedValues2[Pos];
		  ++Pos;
		  Sum2 += CommonFactor*PrecalculatedValues2[Pos];
		  ++Pos;
		}
	      if (Sum2!=0.0)
		File << (X * Factor2) << " " << Real(Sum3/Sum2) << " " << Imag(Sum3/Sum2) << endl;
	      else
		File << (X * Factor2) << " nan" << endl;
	      if (x==NbrPoints-1)
		cout << "Shift = " << Real(Sum3/Sum2) << endl; 
	      X += XInc;
	    }
	  File.close();
	  delete[] PrecalculatedValues;
	  delete [] OutputNameCorr;
	}
    }  
  return 0;
}


