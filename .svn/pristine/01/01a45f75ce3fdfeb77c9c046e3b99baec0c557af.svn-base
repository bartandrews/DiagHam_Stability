#include "Vector/RealVector.h"
#include "Matrix/ComplexMatrix.h"

#include "HilbertSpace/BosonOnDisk.h"
#include "HilbertSpace/FermionOnDisk.h"
#include "HilbertSpace/FermionOnDiskUnlimited.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/FermionOnDiskHaldaneBasis.h"
#include "HilbertSpace/BosonOnDiskHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnDiskHaldaneHugeBasisShort.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "FunctionBasis/ParticleOnDiskFunctionBasis.h"

#include "Tools/FQHEFiles/FQHEOnDiskFileTools.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// get the Hilbert space and the vector state form the input file name
//
// inputState = input file name
// nbrParticles = reference on the number of particles 
// forceMaxMomentum = reference on the forced max momentum
// totalLz = reference on the total Lz value
// statistics = reference on the statistic flag
// space = reference on the pointer to the Hilbert space
// state = reference on the state vector
// referenceFile = name of reference file for the squeezed Hilbert space (0 if non squeezed basis)
// loadHilbert = name of the file where the Hilbert space is stored (0 if none) 
// hugeMode = flag for huge basis mode
// hugeMemory = amount of memory in huge mode
// return value = true if no error occured
bool FQHEDiskCorrelationGetHilbertSpace(char* inputState, int& nbrParticles, int& forceMaxMomentum, int& totalLz, bool& statistics, 
					ParticleOnSphere*& space, RealVector& state, char* referenceFile = 0, char* loadHilbert = 0, bool hugeMode = false, long hugeMemory = 0l);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHEDiskCorrelation" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PlotOptionGroup = new OptionGroup ("plot options");  
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += PlotOptionGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (overriding the one found in the vector file name if greater than 0)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "momentum", "single particle momentum (overriding the one found in the vector file name if greater than 0)", 0, true, 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "force-maxmomentum", "force the maximum single particle momentum to a particular value (overriding the one found in the vector file name if greater than 0)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use boson statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermion statistics");
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the vector file describing the state whose correlation has to be plotted");
  (*SystemGroup) += new SingleStringOption  ('\n', "input-states", "use a file to describe the state as a linear combination");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one body coefficients that are requested to evaluate the correlation profile", false);
  (*SystemGroup) += new BooleanOption ('\n', "correlation-profil", "only compute the correlation profil along x");  
  (*SystemGroup) += new BooleanOption ('\n', "force-xsymmetry", "assume the wave function is invariant under the x <->-x symmetry");  
  (*SystemGroup) += new BooleanOption ('\n', "force-ysymmetry", "assume the wave function is invariant under the y <->-y symmetry");  
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);

  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "use-precomputed", "use precomputed matrix elements to do the plot");
  (*PlotOptionGroup) += new SingleStringOption ('\n', "output", "output file name (default output name replace the .vec extension of the input file with .rhorho.dat)", 0);
  (*PlotOptionGroup) += new SingleIntegerOption ('\n', "nbr-samples", "number of samples in radial direction", 1000, true, 10);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskCorrelation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if ((Manager.GetString("state") == 0) && (Manager.GetString("input-states") == 0))
    {
      cout << "FQHEDiskCorrelation requires an input state" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int TotalLz = Manager.GetInteger("momentum");
  int ForceMaxMomentum = Manager.GetInteger("force-maxmomentum");
  bool CoefficientOnlyFlag = Manager.GetBoolean("coefficients-only");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  char* OutputName = Manager.GetString("output");
  int NbrSamples = Manager.GetInteger("nbr-samples");
  bool Statistics = true;
  bool Plot2DFlag = !(Manager.GetBoolean("correlation-profil"));

  if ((Manager.GetBoolean("boson") == true) || (Manager.GetBoolean("fermion") == true))
    {
      if (Manager.GetBoolean("boson") == true)
        Statistics = false;
      else
        Statistics = true;
    }
  
//  if (Manager.GetString("input-states") == 0)
    {
      RealVector State;
      ParticleOnSphere* Space = 0;
      
      if (HaldaneBasisFlag == true)
	{
	  if (FQHEDiskCorrelationGetHilbertSpace(Manager.GetString("state"), NbrParticles, ForceMaxMomentum, TotalLz, Statistics, 
					     Space, State, Manager.GetString("reference-file"), Manager.GetString("load-hilbert")) == false)
	    return -1;
	}
      else
	{
	  if (FQHEDiskCorrelationGetHilbertSpace(Manager.GetString("state"), NbrParticles, ForceMaxMomentum, TotalLz, Statistics, 
					     Space, State) == false)
	    return -1;
	}
      
      ParticleOnDiskFunctionBasis Basis(ForceMaxMomentum);
      
      double Scale = 2.0 * sqrt((double) ForceMaxMomentum);
      cout << "scale="  << Scale << endl;
      Complex* PrecalculatedValues = new Complex [ForceMaxMomentum + 1];
	  
      for (int m = 0; m <= ForceMaxMomentum; ++m)
	{
	  ParticleOnSphereDensityDensityOperator Operator (Space, m, 0, m, 0);
	  OperatorMatrixElementOperation Operation(&Operator, State, State);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  PrecalculatedValues[m] = Operation.GetScalar();
	}
      ofstream File;
      File.precision(14);
      if (OutputName == 0)
	OutputName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "rhorho.dat");
      File.open(OutputName, ios::binary | ios::out);
      if (CoefficientOnlyFlag == false)
	{
	  Complex Tmp (0.0);
	  double RInc = Scale / ((double) NbrSamples);
	  NbrSamples *= 2;      
	  RealVector Value(2, true);
	  Complex TmpValue;
	  double GaussianWeight;
	  Value[0] = -Scale;
	  
	  if (Plot2DFlag == false)
	    {
	      for (int i = 0; i <= NbrSamples; ++i)
		{
		  Tmp = 0.0;
		  GaussianWeight = exp (-0.5 * (Value[0] * Value[0]));
		  for (int j = 0; j <= ForceMaxMomentum; ++j)
		    {
		      Basis.GetFunctionValue(Value, TmpValue, j);
		      Tmp += PrecalculatedValues[j] *  GaussianWeight * SqrNorm(TmpValue);
		    }
		  File << Value[0] << " " << Tmp.Re << endl;
		  Value[0] += RInc;
		}
	    }
	  else
	    {
	      for (int i = 0; i <= NbrSamples; ++i)
		{
		  Value[1] = -Scale;
		  for (int j = 0; j <= NbrSamples; ++j)
		    {
		      Tmp = 0.0;
		      for (int k = 0; k <= ForceMaxMomentum; ++k)
			{
			  Basis.GetFunctionValue(Value, TmpValue, k);
			  Tmp += PrecalculatedValues[k]  * SqrNorm(TmpValue);
			}
		      Tmp *= exp (-0.5 * ((Value[0] * Value[0]) + (Value[1] * Value[1])));
		      File << Value[0] << " " << Value[1] << " " << Tmp.Re << endl;
		      Value[1] += RInc;
		    }
		  File << endl;
		  Value[0] += RInc;
		}
	    }
	}
      else
	{
	  File << "# correlation coefficients for " << Manager.GetString("state") << endl;
	  if (PrecalculatedValues[ForceMaxMomentum].Re != 0.0)
	    {
	      File << "#" << endl << "# m    n_m    (n_m / n_Nphi)" << endl;
	      for (int i = 0; i <= ForceMaxMomentum; ++i)
		File << i << " " << PrecalculatedValues[i].Re << " " << (PrecalculatedValues[i].Re / PrecalculatedValues[ForceMaxMomentum].Re) << endl;
	    }
	  else
	    {
	      File << "#" << endl << "# m    n_m" << endl;
	      for (int i = 0; i <= ForceMaxMomentum; ++i)
		File << i << " " << PrecalculatedValues[i].Re << endl;
	    }
	}
      File.close();

      delete[] PrecalculatedValues;
      delete Space;
    }
//   else
//     {
//       MultiColumnASCIIFile InputVectors;
//       if (InputVectors.Parse(Manager.GetString("input-states")) == false)
// 	{
// 	  InputVectors.DumpErrors(cout) << endl;
// 	  return -1;
// 	}
//       if (FQHEOnDiskFindSystemInfoFromFileName(InputVectors(0, 0), NbrParticles, ForceMaxMomentum, TotalLz, Statistics) == false)
// 	{
// 	  return -1;      
// 	}
//       for (int i = 1; i < InputVectors.GetNbrLines(); ++i)
// 	{
// 	  int TmpNbrParticles = 0;
// 	  int TmpForceMaxMomentum = 0;
// 	  int TmpTotalLz = 0;
// 	  if (FQHEOnDiskFindSystemInfoFromFileName(InputVectors(0, i), TmpNbrParticles, TmpForceMaxMomentum, TmpTotalLz, Statistics) == false)
// 	    {
// 	      return -1;      
// 	    }
// 	  if (NbrParticles != TmpNbrParticles)
// 	    {
// 	      cout << InputVectors(0, i) << " and " << InputVectors(0, 0) << " have different numbers of particles" << endl;
// 	      return -1;
// 	    }
// 	  if (ForceMaxMomentum != TmpForceMaxMomentum)
// 	    {
// 	      cout << InputVectors(0, i) << " and " << InputVectors(0, 0) << " have different forced max momenta" << endl;
// 	      return -1;
// 	    }
// 	}
//       double* Coefficients = 0;
//       if (InputVectors(1, 0) != 0)
// 	{
// 	  Coefficients = InputVectors.GetAsDoubleArray(1);
// 	}
//       else
// 	{
// 	  cout << "no coefficients defined in " << Manager.GetString("input-states") << endl;
// 	  return -1; 
// 	}

//       ParticleOnDiskFunctionBasis Basis(TotalLz);
      
//       double Scale = 2.0 * sqrt((double) TotalLz);
//       ComplexMatrix PrecalculatedValues(ForceMaxMomentum + 1, ForceMaxMomentum + 1, true);
//       ComplexMatrix** RawPrecalculatedValues = new ComplexMatrix*[InputVectors.GetNbrLines()];
//       for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
// 	{
// 	  RawPrecalculatedValues[i] = new ComplexMatrix[InputVectors.GetNbrLines()];
// 	  for (int j = i; j < InputVectors.GetNbrLines(); ++j)
// 	    RawPrecalculatedValues[i][j] = ComplexMatrix(ForceMaxMomentum + 1, ForceMaxMomentum + 1, true);
// 	}

//       if (Manager.GetString("use-precomputed") == 0)
// 	{	  
// 	  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
// 	    {
// 	      ParticleOnSphere* LeftSpace = 0;
// 	      RealVector LeftState;
// 	      int LeftTotalLz = 0;
// 	      if ((InputVectors.GetNbrColumns() <= 2) || (strcmp ("none", InputVectors(2, i)) == 0))
// 		{
// 		  if (FQHEDiskCorrelationGetHilbertSpace(InputVectors(0, i), NbrParticles, ForceMaxMomentum, LeftTotalLz, Statistics, 
// 						     LeftSpace, LeftState) == false)
// 		    return -1;
// 		}
// 	      else
// 		{
// 		  if ((InputVectors.GetNbrColumns() <= 3) || (strcmp ("none", InputVectors(3, i)) == 0))
// 		    {
// 		      if (FQHEDiskCorrelationGetHilbertSpace(InputVectors(0, i), NbrParticles, ForceMaxMomentum, LeftTotalLz, Statistics, 
// 							 LeftSpace, LeftState, InputVectors(2, i)) == false)
// 			return -1;
// 		    }
// 		  else
// 		    {
// 		      if (FQHEDiskCorrelationGetHilbertSpace(InputVectors(0, i), NbrParticles, ForceMaxMomentum, LeftTotalLz, Statistics, 
// 							 LeftSpace, LeftState, InputVectors(2, i), InputVectors(3, i)) == false)
// 			return -1;
// 		    }
// 		}
// 	      for (int m1 = 0; m1 <= ForceMaxMomentum; ++m1)
// 		{
// 		  for (int m2 = 0; m2 <= ForceMaxMomentum; ++m2)
// 		    {
// 		      for (int m3 = 0; m3 <= ForceMaxMomentum; ++m3)
// 			{
// 			  int m4 = m1 + m2 - m3;
// 			  if ((m4 >= 0) && (m4 <= ForceMaxMomentum))
// 			    {
// 			      ParticleOnSphereDensityDensityOperator Operator (LeftSpace, m1, m2, m3, m4);
// 			      OperatorMatrixElementOperation Operation(&Operator, LeftState, LeftState);
// 			      Operation.ApplyOperation(Architecture.GetArchitecture());
// 			      Complex Tmp = Operation.GetScalar();
// 			      RawPrecalculatedValues[i][i].AddToMatrixElement(m, m, Tmp);
// 			      Tmp *= (Coefficients[i] * Coefficients[i]);
// 			      PrecalculatedValues.AddToMatrixElement(m, m, Tmp);
// 			    }
// 			  for (int j = i + 1; j < InputVectors.GetNbrLines(); ++j)
// 			    {
// 			      ParticleOnSphere* RightSpace = 0;
// 			      RealVector RightState;
// 			      int RightTotalLz = 0;
// 			      if ((InputVectors.GetNbrColumns() <= 2) || (strcmp ("none", InputVectors(2, j)) == 0))
// 				{
// 				  if (FQHEDiskCorrelationGetHilbertSpace(InputVectors(0, j), NbrParticles, ForceMaxMomentum, RightTotalLz, Statistics, 
// 									 RightSpace, RightState) == false)
// 				    return -1;
// 				}
// 			      else
// 				{
// 				  if ((InputVectors.GetNbrColumns() <= 3) || (strcmp ("none", InputVectors(3, j)) == 0))
// 				    {
// 				      if (FQHEDiskCorrelationGetHilbertSpace(InputVectors(0, j), NbrParticles, ForceMaxMomentum, RightTotalLz, Statistics, 
// 									     RightSpace, RightState, InputVectors(2, j)) == false)
// 					return -1;
// 				    }
// 				  else
// 				    {
// 				      if (FQHEDiskCorrelationGetHilbertSpace(InputVectors(0, j), NbrParticles, ForceMaxMomentum, RightTotalLz, Statistics, 
// 									     RightSpace, RightState, InputVectors(2, j), InputVectors(3, j)) == false)
// 					return -1;
// 				    }
// 				}
// 		  for (int m = 0; m <= ForceMaxMomentum; ++m)
// 		    {
// 		      for (int n = 0; n <= ForceMaxMomentum; ++n)
// 			{
// 			  if ((RightTotalLz - n) == (LeftTotalLz - m))
// 			    {
// 			      RightSpace->SetTargetSpace(LeftSpace);
// 			      ParticleOnSphereDensityDensityOperator Operator (RightSpace, m, n);
// 			      OperatorMatrixElementOperation Operation(&Operator, LeftState, RightState);
// 			      Operation.ApplyOperation(Architecture.GetArchitecture());
// 			      Complex Tmp = Operation.GetScalar();
// 			      RawPrecalculatedValues[i][j].AddToMatrixElement(m, n, Tmp);
// 			      Tmp *= (Coefficients[i] * Coefficients[j]);
// 			      PrecalculatedValues.AddToMatrixElement(m, n, Tmp);
// 			      PrecalculatedValues.AddToMatrixElement(n, m, Conj(Tmp));
// 			    }
// 			}
// 		    }
// 		  delete RightSpace;
// 		}
// 	    }
// 	}
//       else
// 	{
// 	  MultiColumnASCIIFile PrecomputedElements;
// 	  if (PrecomputedElements.Parse(Manager.GetString("use-precomputed")) == false)
// 	    {
// 	      PrecomputedElements.DumpErrors(cout) << endl;
// 	      return -1;
// 	    }
// 	  int* StateIndexLeft = PrecomputedElements.GetAsIntegerArray(0);
// 	  int* StateIndexRight = PrecomputedElements.GetAsIntegerArray(1);
// 	  int* OperatorIndexLeft = PrecomputedElements.GetAsIntegerArray(2);
// 	  int* OperatorIndexRight = PrecomputedElements.GetAsIntegerArray(3);	  
// 	  Complex* OneBodyCoefficients = PrecomputedElements.GetAsComplexArray(4);	  
// 	  for (int i = 0; i < PrecomputedElements.GetNbrLines(); ++i)
// 	    {
// 	      RawPrecalculatedValues[StateIndexLeft[i]][StateIndexRight[i]].AddToMatrixElement(OperatorIndexLeft[i], OperatorIndexRight[i], OneBodyCoefficients[i]);
// 	    }
// 	  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
// 	    {
// 	      for (int m = 0; m <= ForceMaxMomentum; ++m)
// 		{
// 		  Complex Tmp = RawPrecalculatedValues[i][i][m][m];
// 		  Tmp *= (Coefficients[i] * Coefficients[i]);
// 		  PrecalculatedValues.AddToMatrixElement(m, m, Tmp);
// 		}
// 	      for (int j = i + 1; j < InputVectors.GetNbrLines(); ++j)
// 		{
// 		  for (int m = 0; m <= ForceMaxMomentum; ++m)
// 		    {
// 		      for (int n = 0; n <= ForceMaxMomentum; ++n)
// 			{
// 			  Complex Tmp = RawPrecalculatedValues[i][j][m][n];
// 			  Tmp *= (Coefficients[i] * Coefficients[j]);
// 			  PrecalculatedValues.AddToMatrixElement(m, n, Tmp);
// 			  PrecalculatedValues.AddToMatrixElement(n, m, Conj(Tmp));
// 			}
// 		    }
// 		}
// 	    }
// 	}
//       ofstream File;
//       File.precision(14);
//       if (OutputName == 0)
// 	OutputName = ReplaceExtensionToFileName(Manager.GetString("input-states"), "dat", "rhorho.dat");
//       File.open(OutputName, ios::binary | ios::out);
//       File << "# correlation coefficients for " << Manager.GetString("input-states") << endl;
//       if (CoefficientOnlyFlag == false)
// 	{
// 	  File << "#" << endl << "# m  n  c_{m,n}" << endl;
// 	  for (int i = 0; i <= ForceMaxMomentum; ++i)
// 	    for (int j = 0; j <= ForceMaxMomentum; ++j)
// 	      File << "# " << i << " " << j << " " << PrecalculatedValues[i][j] << endl;
// 	}
//       else
// 	{
// 	  File << "#" << endl << "# state_index_left state_index_right m  n  c_{m,n}" << endl;
// 	  for (int m = 0; m < InputVectors.GetNbrLines(); ++m)
// 	    for (int n = m; n < InputVectors.GetNbrLines(); ++n)
// 	      for (int i = 0; i <= ForceMaxMomentum; ++i)
// 		for (int j = 0; j <= ForceMaxMomentum; ++j)
// 		  {
// 		    File << m << " " << n << " " << i << " " << j << " " << RawPrecalculatedValues[m][n][i][j] << endl;
// 		  }
// 	}

//       if (CoefficientOnlyFlag == false)
// 	{
// 	  Complex Tmp (0.0);
// 	  double RInc = Scale / ((double) NbrSamples);
// 	  NbrSamples *= 2;      
// 	  RealVector Value(2, true);
// 	  Complex TmpValue;
// 	  Complex TmpValue2;
// 	  Value[0] = -Scale;

// 	  for (int i = 0; i <= NbrSamples; ++i)
// 	    {
// 	      Value[1] = -Scale;
// 	      for (int j = 0; j <= NbrSamples; ++j)
// 		{
// 		  Tmp = 0.0;
// 		  for (int k = 0; k <= ForceMaxMomentum; ++k)
// 		    {
// 		      Basis.GetFunctionValue(Value, TmpValue, k);
// 		      for (int l = 0; l <= ForceMaxMomentum; ++l)
// 			{
// 			  Basis.GetFunctionValue(Value, TmpValue2, l);
// 			  Tmp += PrecalculatedValues[k][l]  * Conj(TmpValue) * TmpValue2;
// 			}
// 		    }
// 		  Tmp *= exp (-0.5 * ((Value[0] * Value[0]) + (Value[1] * Value[1])));
// 		  File << Value[0] << " " << Value[1] << " " << Tmp.Re << endl;
// 		  Value[1] += RInc;
// 		}
// 	      File << endl;
// 	      Value[0] += RInc;
// 	    }
// 	}
//       File.close();
//     }
}

// get the Hilbert space and the vector state form the input file name
//
// inputState = input file name
// nbrParticles = reference on the number of particles 
// forceMaxMomentum = reference on the forced max momentum
// totalLz = reference on the total Lz value
// statistics = reference on the statistic flag
// space = reference on the pointer to the Hilbert space
// state = reference on the state vector
// referenceFile = name of reference file for the squeezed Hilbert space (0 if non squeezed basis)
// loadHilbert = name of the file where the Hilbert space is stored (0 if none) 
// hugeMode = flag for huge basis mode
// hugeMemory = amount of memory in huge mode
// return value = true if no error occured

bool FQHEDiskCorrelationGetHilbertSpace(char* inputState, int& nbrParticles, int& forceMaxMomentum, int& totalLz, bool& statistics, 
				    ParticleOnSphere*& space, RealVector& state, char* referenceFile, char* loadHilbert, bool hugeMode, long hugeMemory)
{
  if (FQHEOnDiskFindSystemInfoFromFileName(inputState, nbrParticles, forceMaxMomentum, totalLz, statistics) == false)
    {
      return false;      
    }
  if (state.ReadVector (inputState) == false)
    {
      cout << "can't open vector file " << inputState << endl;
      return false;      
    }

  int TmpMaxMomentum = 0;
  if (statistics == false)
    {
      TmpMaxMomentum = totalLz;
      if ((forceMaxMomentum >= 0) && (forceMaxMomentum < TmpMaxMomentum))
	TmpMaxMomentum = forceMaxMomentum; 
    }
  else
    {
      TmpMaxMomentum = (totalLz - (((nbrParticles - 1) * (nbrParticles - 2)) / 2));
      if ((forceMaxMomentum >= 0) && (forceMaxMomentum < TmpMaxMomentum))
	TmpMaxMomentum = forceMaxMomentum;
    }     

  if (referenceFile == 0)
    {
      if (statistics == true)
	{
#ifdef __64_BITS__
      if (TmpMaxMomentum < 63)      
#else
	if (TmpMaxMomentum < 31)
#endif
	  space = new FermionOnDisk (nbrParticles, totalLz, TmpMaxMomentum);
	else
	  space = new FermionOnDiskUnlimited (nbrParticles, totalLz, TmpMaxMomentum);      
	}
      else
	{
	  space = new BosonOnDisk(nbrParticles, totalLz, forceMaxMomentum);
	}
    }
  else
    {
      int* ReferenceState = 0;
      ConfigurationParser ReferenceStateDefinition;
      if (ReferenceStateDefinition.Parse(referenceFile) == false)
	{
	  ReferenceStateDefinition.DumpErrors(cout) << endl;
	  return false;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", nbrParticles) == false) || (nbrParticles <= 0))
	{
	  cout << "NbrParticles is not defined or as a wrong value" << endl;
	  return false;
	}
      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", forceMaxMomentum) == false) || (forceMaxMomentum <= 0))
	{
	  cout << "LzMax is not defined or as a wrong value" << endl;
	  return false;
	}
      TmpMaxMomentum = forceMaxMomentum;
      int MaxNbrLz;
      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	{
	  cout << "error while parsing ReferenceState in " << referenceFile << endl;
	  return false;     
	}
      if (MaxNbrLz != (forceMaxMomentum + 1))
	{
	  cout << "wrong LzMax value in ReferenceState" << endl;
	  return false;     
	}
      totalLz = 0;
      for (int i = 1; i <= forceMaxMomentum; ++i)
	totalLz += i * ReferenceState[i];
      if (statistics == true)
	{
	  if (loadHilbert == 0)
	    space = new FermionOnDiskHaldaneBasis (nbrParticles, totalLz, forceMaxMomentum, ReferenceState);
	  else
	    space = new FermionOnDiskHaldaneBasis (loadHilbert);
	}
      else
	{
#ifdef  __64_BITS__
	  if ((forceMaxMomentum + nbrParticles - 1) < 63)
#else
	    if ((forceMaxMomentum + nbrParticles - 1) < 31)	
#endif

	      {	  
		if (hugeMode == true)
		  {
		    if (loadHilbert == 0)
		      {
			cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
			return false;
		      }
		    space = new BosonOnDiskHaldaneHugeBasisShort (loadHilbert, hugeMemory);
		  }
		else
		  {
		    if (loadHilbert == 0)
		      space = new BosonOnDiskHaldaneBasisShort(nbrParticles, totalLz, forceMaxMomentum, ReferenceState);
		    else
		      space = new BosonOnDiskHaldaneBasisShort(loadHilbert);
		  }
	      }
	}
    }


  if (space->GetLargeHilbertSpaceDimension() != state.GetLargeVectorDimension())
    {
      cout << "dimension mismatch between the state (" << state.GetLargeVectorDimension() << ") and the Hilbert space (" << space->GetLargeHilbertSpaceDimension() << ")" << endl;
      return false;
    }
  return true;
}
