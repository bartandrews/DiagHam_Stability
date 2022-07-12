#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/WignerSmallDMatrix.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereODLRORotation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('r', "right-state", "file describing right state decomposition into L^2 eigenstates");
  (*SystemGroup) += new SingleStringOption  ('l', "left-state", "file describing left state decomposition into L^2 eigenstates");
  (*OutputGroup) += new SingleIntegerOption ('z', "lz-value", "twice the Lz value of the left state", 0);

  (*OutputGroup) += new SingleIntegerOption ('n', "nbr-points", "number of points that has to be computed", 100);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereODLRORotation -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int LzValue = Manager.GetInteger("lz-value");

  MultiColumnASCIIFile LeftStateDecomposition;
  if (LeftStateDecomposition.Parse(Manager.GetString("left-state")) == false)
    {
      LeftStateDecomposition.DumpErrors(cout);
      return -1;
    }
  MultiColumnASCIIFile RightStateDecomposition;
  if (RightStateDecomposition.Parse(Manager.GetString("right-state")) == false)
    {
      RightStateDecomposition.DumpErrors(cout);
      return -1;
    }

  double* LeftStateComponents = LeftStateDecomposition.GetAsDoubleArray(1);
  double* RightStateComponents = RightStateDecomposition.GetAsDoubleArray(1);
  int* LeftStateLValues = LeftStateDecomposition.GetAsIntegerArray(2);
  int* RightStateLValues = RightStateDecomposition.GetAsIntegerArray(2);
  int NbrLValues = 0;
  for (int i = 0; i < LeftStateDecomposition.GetNbrLines(); ++i)
    {
      int CurrentLValue = LeftStateLValues[i];
      int j = 0;
      for (; (j < RightStateDecomposition.GetNbrLines()) && (RightStateLValues[j] != CurrentLValue); ++j);
      if (j < RightStateDecomposition.GetNbrLines())
	++NbrLValues;	
    }

  int* LValues = new int [NbrLValues];
  double* Coefficients = new double [NbrLValues];
  WignerSmallDMatrix** WignerCoefficients = new WignerSmallDMatrix*[NbrLValues];
  NbrLValues = 0;
  for (int i = 0; i < LeftStateDecomposition.GetNbrLines(); ++i)
    {
      int CurrentLValue = LeftStateLValues[i];
      int j = 0;
      for (; (j < RightStateDecomposition.GetNbrLines()) && (RightStateLValues[j] != CurrentLValue); ++j);
      if (j < RightStateDecomposition.GetNbrLines())
	{
	  LValues[NbrLValues] = CurrentLValue;
	  WignerCoefficients[NbrLValues] = new WignerSmallDMatrix(CurrentLValue);
	  RealVector LeftVector;
	  if (LeftVector.ReadVector(LeftStateDecomposition(0, i)) == false)
	    {
	      cout << "can't open vector file " << LeftStateDecomposition(0, i) << endl;
	      return -1;      	      
	    }
	  RealVector RightVector;
	  if (RightVector.ReadVector(RightStateDecomposition(0, i)) == false)
	    {
	      cout << "can't open vector file " << RightStateDecomposition(0, i) << endl;
	      return -1;      	      
	    }
	  Coefficients[NbrLValues] = (LeftVector * RightVector) * (LeftStateComponents[i] * (RightStateComponents[i]));
	  ++NbrLValues;	
	}
    }

  cout.precision(14);
  int NbrPoints = Manager.GetInteger("nbr-points");
  double Theta = 0.0;
  double ThetaInc = M_PI / ((double) NbrPoints);
  for (int i = 0; i <= NbrPoints; ++i)
    {
      double Tmp = 0.0;
      for (int j = 0; j < NbrLValues; ++j)
	{
	  double Tmp2 = Coefficients[j] * (*(WignerCoefficients[j]))(LzValue, LzValue, Theta);// * (*(WignerCoefficients[j]))(LzValue, LzValue, Theta);
	  Tmp += Tmp2;
	}      
      cout << Theta << " " << Tmp << endl;
      Theta += ThetaInc;
    }

  return 0;
}


