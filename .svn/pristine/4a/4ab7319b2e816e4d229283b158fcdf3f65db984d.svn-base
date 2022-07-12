#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// compute the weight of a given orbital for a sharp real space cut
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutPosition = x position of the cut
// return value = square of the orbital weight
double FQHECylinderComputeSharpRealSpaceCutCoefficient (double orbitalIndex, double perimeter, double cutPosition);

// compute the weight of a given orbital for a sharp real space cut with patch having a finite length along the cylinder axis
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutInitialPosition = leftmost x position of the cut
// cutFinalPosition = rightmost x position of the cut
// return value = square of the orbital weight
double FQHECylinderComputeSharpRealSpaceCutCoefficient (double orbitalIndex, double perimeter, double cutInitialPosition, double cutFinalPosition);

// compute the weight for a sharp real space cut with patch having a finite length along the cylinder axis 
// and a finite width (centered around y=0)
//
// orbitalIndex1 = index of the first orbital (can be negative, zero-th orbital being centered at x=0)
// orbitalIndex2 = index of the second orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutInitialPosition = leftmost x position of the cut
// cutFinalPosition = rightmost x position of the cut
// width = patch width along the cylinder perimeter
// return value = square of the orbital weight
double FQHECylinderComputeSharpRealSpaceCutCoefficient (double orbitalIndex1, double orbitalIndex2, double perimeter, double cutInitialPosition, double cutFinalPosition, double width);

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHECylinderRealSpacePartitionCoefficients" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (if zero, assume an infinite cylinder)", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "error", "error below which a coefficient is consider as 0 (or 1)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the cylinder", 1);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "flux-insertion", "additional flux insertion along the cylinder axis (in phi0 units)", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "finite-patch", "use a patch with a finite length and optionally a finite width");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "patch-center", "position of the patch center along the cylinder axis", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "patch-length", "length of the patch along the cylinder axis", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "patch-width", "width of the patch along the cylinder perimeter (0 if it should extend along the whole cylinder perimeter)", 0.0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "optional output file name (default is realspace_cylinder_l_*_perimeter_*_2s_*.dat)");
  (*OutputGroup) += new BooleanOption ('\n', "column-output", "when having a cut preserving the translation along the cylinder perimeter, use a column formatted output instead of a single line");
  (*OutputGroup) += new BooleanOption ('\n', "show-overlap", " show the overlap for the region A  when using a finite size patch");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderRealSpacePartitionCoefficients -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");
  double Error = Manager.GetDouble("error");
  if ((NbrFluxQuanta == 0) && (Error == 0.0))
    {
      Error = MACHINE_PRECISION;
    }

  double Perimeter = Manager.GetDouble("cylinder-perimeter");
  if (Perimeter == 0.0)
    {
      if (NbrFluxQuanta == 0)
	{
	  cout << "error, nbr-flux has to be provided if cylinder-perimeter is zero" << endl;
	  return 0;
	}
      else
	{
	  Perimeter = sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) * Manager.GetDouble("aspect-ratio"));
	}
    }
  char* OutputFile = 0;

  if ((Manager.GetBoolean("finite-patch") == false) || (Manager.GetDouble("patch-width") == 0.0))
    {
      double CutPosition = 0.0;
      double CutLength = 0.0;
      double FluxInsertion = Manager.GetDouble("flux-insertion");
      if (Manager.GetString("output-file") == 0)
	{
	  if (FluxInsertion == 0.0)
	    {
	      if (Manager.GetBoolean("finite-patch") == true)
		{
		  OutputFile = new char[512];
		  sprintf (OutputFile, "realspace_cylinder_l_%.6f_perimeter_%.6f_xcenter_%.6f_length_%.6f_2s_%d.dat", CutPosition, Perimeter, 
			   Manager.GetDouble("patch-center"), Manager.GetDouble("patch-length"), NbrFluxQuanta);
		}
	      else
		{
		  OutputFile = new char[512];
		  sprintf (OutputFile, "realspace_cylinder_l_%.6f_perimeter_%.6f_2s_%d.dat", CutPosition, Perimeter, NbrFluxQuanta);
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("finite-patch") == true)
		{
		  OutputFile = new char[512];
		  sprintf (OutputFile, "realspace_cylinder_l_%.6f_perimeter_%.6f_xcenter_%.6f_length_%.6f_flux_%.6f_2s_%d.dat", CutPosition, Perimeter, 
			   Manager.GetDouble("patch-center"), Manager.GetDouble("patch-length"), FluxInsertion, NbrFluxQuanta);
		}
	      else
		{
		  OutputFile = new char[512];
		  sprintf (OutputFile, "realspace_cylinder_l_%.6f_perimeter_%.6f_flux_%.6f_2s_%d.dat", CutPosition, Perimeter, FluxInsertion, NbrFluxQuanta);
		}
	    }
	}
      else
	{
	  OutputFile = new char[strlen(Manager.GetString("output-file")) + 1];
	  strcpy (OutputFile, Manager.GetString("output-file"));
	}

      ofstream File;
      File.open(OutputFile, ios::binary | ios::out);
      File.precision(14);
      if (Manager.GetBoolean("finite-patch") == false)
	{
	  File << "# real space coefficients for a sharp cut at x=" << CutPosition << " on ";
	}
      else
	{
	  CutPosition = Manager.GetDouble("patch-center");
	  CutLength = Manager.GetDouble("patch-length");
	  File << "# real space coefficients for a sharp cut centered at x=" << CutPosition << " and of length " << CutLength << " on ";
	}
      if (NbrFluxQuanta == 0)
	{
	  File << "an infinite cylinder with perimeter L=" << Perimeter;
	}
      else
	{
	  File << "a cylinder with perimeter L=" << Perimeter << " and N_phi=" << NbrFluxQuanta;
	}
      if (FluxInsertion != 0.0)
	{
	  File << " and flux insertion phi=" << FluxInsertion;
	}
      int NbrCoefficients = 0;
      double* Coefficients = 0;
      if (NbrFluxQuanta == 0)
	{
	  int MaxNbrCoefficients = 1000; 
	  double* TmpCoefficients = new double [MaxNbrCoefficients];
	  if (Manager.GetBoolean("finite-patch") == false)
	    {
	      for (int i = 0; i < MaxNbrCoefficients; ++i)
		{
		  TmpCoefficients[i] = FQHECylinderComputeSharpRealSpaceCutCoefficient(((double) i) + FluxInsertion, Perimeter, CutPosition);
		  if (TmpCoefficients[i] < Error)
		    {
		      i = MaxNbrCoefficients;
		    }
		  else
		    {
		      ++NbrCoefficients;
		    }	    
		}
	    }
	  else
	    {
	      for (int i = 0; i < MaxNbrCoefficients; ++i)
		{
		  TmpCoefficients[i] = FQHECylinderComputeSharpRealSpaceCutCoefficient(((double) i) + FluxInsertion, Perimeter, CutPosition - (0.5 * CutLength), 
										       CutPosition + (0.5 * CutLength));
		  if (TmpCoefficients[i] < Error)
		    {
		      i = MaxNbrCoefficients;
		    }
		  else
		    {
		      ++NbrCoefficients;
		    }
		}
	    }
	  Coefficients = new double [2 * NbrCoefficients - 1];
	  Coefficients[NbrCoefficients - 1] =  TmpCoefficients[0];
	  for (int i = 1; i < NbrCoefficients; ++i)
	    {
	      Coefficients[NbrCoefficients - 1 + i] =  TmpCoefficients[i];
	      Coefficients[NbrCoefficients - 1 - i] =  1.0 - TmpCoefficients[i];
	    }
	  delete[] TmpCoefficients;
	}
      else
	{
	  Coefficients = new double [NbrFluxQuanta + 1];
	  if (Manager.GetBoolean("finite-patch") == false)
	    {
	      for (NbrCoefficients = 0; NbrCoefficients <= NbrFluxQuanta; ++NbrCoefficients)
		{
		  Coefficients[NbrCoefficients] = FQHECylinderComputeSharpRealSpaceCutCoefficient(((double) NbrCoefficients) - 0.5 * ((double) NbrFluxQuanta) + FluxInsertion, Perimeter, CutPosition);
		}
	    }
	  else
	    {
	      for (NbrCoefficients = 0; NbrCoefficients <= NbrFluxQuanta; ++NbrCoefficients)
		{
		  Coefficients[NbrCoefficients] = FQHECylinderComputeSharpRealSpaceCutCoefficient(((double) NbrCoefficients) - 0.5 * ((double) NbrFluxQuanta) + FluxInsertion, Perimeter, 
												  CutPosition - (0.5 * CutLength), CutPosition + (0.5 * CutLength));
		}
	    }
	}
      if (Manager.GetBoolean("column-output") == false)
	{
	  File << endl << "OrbitalSquareWeights =";
	  for (int i = 0; i < NbrCoefficients; ++i)
	    File << " " << Coefficients[i];
	  File << endl;
	}
      else
	{
	  File << endl;
	  for (int i = 0; i < NbrCoefficients; ++i)
	    {
	      for (int j = 0; j < NbrCoefficients; ++j)
		{
		  if (i == j)
		    {
		      File << i << " " << j << " " << sqrt(Coefficients[i]) << " " << sqrt(1.0 - Coefficients[i]) << endl;
		    }
		  else
		    {
		      File << i << " " << j << " " << 0.0 << " " << 0.0 << endl;
		    }
		}
	    }
	}
      File.close();
      delete[] Coefficients;
    }
  else
    {
      double CutPosition = Manager.GetDouble("patch-center");
      double CutLength = Manager.GetDouble("patch-length");
      double CutWidth = Manager.GetDouble("patch-width");
      if (Manager.GetString("output-file") == 0)
	{
	  OutputFile = new char[512];
	  sprintf (OutputFile, "realspace_cylinder_l_%.6f_perimeter_%.6f_xcenter_%.6f_length_%.6f_width_%.6f_2s_%d.dat", CutPosition, Perimeter, 
		   CutPosition, CutLength, CutWidth, NbrFluxQuanta);
	}
      else
	{
	  OutputFile = new char[strlen(Manager.GetString("output-file")) + 1];
	  strcpy (OutputFile, Manager.GetString("output-file"));
	}

      ofstream File;
      File.open(OutputFile, ios::binary | ios::out);
      File.precision(14);
      File << "# real space coefficients for a sharp cut centered at x=" << CutPosition << " and of length=" << CutLength << " and width=" << CutWidth << " on ";
      if (NbrFluxQuanta == 0)
	{
	  File << "an infinite cylinder with perimeter L=" << Perimeter;
	}
      else
	{
	  File << "a cylinder with perimeter L=" << Perimeter << " and N_phi=" << NbrFluxQuanta;
	}
      File << endl;
      if (NbrFluxQuanta == 0)
	{
	}
      else
	{
	  RealSymmetricMatrix TmpOverlapMatrix (NbrFluxQuanta + 1, true);
	  for (int i = 0; i <= NbrFluxQuanta; ++i)
	    {
	      TmpOverlapMatrix.SetMatrixElement(i, i, (CutWidth / Perimeter) *  FQHECylinderComputeSharpRealSpaceCutCoefficient(((double) i) - 0.5 * ((double) NbrFluxQuanta), Perimeter, 
																CutPosition - (0.5 * CutLength), CutPosition + (0.5 * CutLength)));
	      for (int j = i + 1; j <= NbrFluxQuanta; ++j)
		{
		  TmpOverlapMatrix.SetMatrixElement(i, j, FQHECylinderComputeSharpRealSpaceCutCoefficient(((double) i) - 0.5 * ((double) NbrFluxQuanta), 
													  ((double) j) - 0.5 * ((double) NbrFluxQuanta), Perimeter, 
													  CutPosition - (0.5 * CutLength), CutPosition + (0.5 * CutLength), CutWidth));
		}
	    }
	  if (Manager.GetBoolean("show-overlap") == true)
	    {
	      cout.precision(14);
	      cout << TmpOverlapMatrix << endl;
	    }
	  RealMatrix TmpTransformationMatrix1 (NbrFluxQuanta + 1, NbrFluxQuanta + 1);
	  TmpTransformationMatrix1.SetToIdentity();
	  RealDiagonalMatrix TmpDiag1 = RealDiagonalMatrix(NbrFluxQuanta + 1, true);  
	  RealDiagonalMatrix TmpDiag2 = RealDiagonalMatrix(NbrFluxQuanta + 1, true);  
#ifdef __LAPACK__
	  TmpOverlapMatrix.LapackDiagonalize(TmpDiag1, TmpTransformationMatrix1);
#else
	  TmpOverlapMatrix.Diagonalize(TmpDiag1, TmpTransformationMatrix1);
#endif 	  
	  for (int i = 0; i <= NbrFluxQuanta; ++i)
	    {
	      TmpDiag2[i] = sqrt(1.0 - TmpDiag1[i]);	      
	      TmpDiag1[i] = sqrt(TmpDiag1[i]);
	    }
// 	  cout << TmpTransformationMatrix1 << endl;
// 	  cout << TmpDiag1 <<endl;
	  RealMatrix TmpTransformationMatrix2 = TmpTransformationMatrix1.DuplicateAndTranspose();
 	  RealMatrix TmpTransformationMatrix3 = TmpTransformationMatrix1 * TmpDiag2;
	  for (int i = 0; i < TmpTransformationMatrix1.GetNbrColumn(); ++i)
	    {
	      TmpTransformationMatrix1[i] *= TmpDiag1[i];
	    }
	  RealMatrix TmpTransformationMatrix4 = TmpTransformationMatrix1 * TmpTransformationMatrix2;
	  RealMatrix TmpTransformationMatrix5 = TmpTransformationMatrix3 * TmpTransformationMatrix2;
//	  cout << TmpTransformationMatrix4 << endl;
	  for (int i = 0; i <= NbrFluxQuanta; ++i)
	    {
	      for (int j = 0; j <= NbrFluxQuanta; ++j)
		{
		  double Tmp1;
		  double Tmp2;
		  TmpTransformationMatrix4.GetMatrixElement(i, j, Tmp1);
		  TmpTransformationMatrix5.GetMatrixElement(i, j, Tmp2);
		  File << i << " " << j << " " << Tmp1 << " " << Tmp2 << endl;
		}
	    }
	}
      File.close();
    }
  return 0;
}


// compute the weight of a given orbital for a sharp real space cut
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutPosition = x position of the cut
// return value = square of the orbital weight

double FQHECylinderComputeSharpRealSpaceCutCoefficient (double orbitalIndex, double perimeter, double cutPosition)
{  
  return (0.5 * (1.0 + erf(cutPosition - (orbitalIndex * 2.0 * M_PI / perimeter))));
}

// compute the weight of a given orbital for a sharp real space cut with patch having a finite length along the cylinder axis
//
// orbitalIndex = index of the orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutInitialPosition = leftmost x position of the cut
// cutFinalPosition = rightmost x position of the cut
// return value = square of the orbital weight

double FQHECylinderComputeSharpRealSpaceCutCoefficient (double orbitalIndex, double perimeter, double cutInitialPosition, double cutFinalPosition)
{  
  return (0.5 * (erf(cutFinalPosition - (orbitalIndex * 2.0 * M_PI / perimeter))
		 - erf(cutInitialPosition - (orbitalIndex * 2.0 * M_PI / perimeter))));
}

// compute the weight for a sharp real space cut with patch having a finite length along the cylinder axis 
// and a finite width (centered around y=0)
//
// orbitalIndex1 = index of the first orbital (can be negative, zero-th orbital being centered at x=0)
// orbitalIndex2 = index of the second orbital (can be negative, zero-th orbital being centered at x=0)
// perimeter = cylinder perimeter
// cutInitialPosition = leftmost x position of the cut
// cutFinalPosition = rightmost x position of the cut
// width = patch width along the cylinder perimeter
// return value = square of the orbital weight

double FQHECylinderComputeSharpRealSpaceCutCoefficient (double orbitalIndex1, double orbitalIndex2, double perimeter, double cutInitialPosition, double cutFinalPosition, double width)
{  
  return (0.5 * (erf(cutFinalPosition - (orbitalIndex1 + orbitalIndex2) * M_PI / perimeter)
		 - erf(cutInitialPosition - (orbitalIndex1 + orbitalIndex2) * M_PI / perimeter))
	  * exp (-(M_PI * M_PI) / (perimeter * perimeter) * (orbitalIndex1 - orbitalIndex2) * (orbitalIndex1 - orbitalIndex2))
	  * sin (M_PI * width / perimeter * (orbitalIndex1 - orbitalIndex2)) / (M_PI * (orbitalIndex1 - orbitalIndex2)));
}
