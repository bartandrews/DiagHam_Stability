#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>



using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;



int main(int argc, char** argv)
{

  OptionManager Manager ("FQHESphereFermionsWithSpinEntanglementEntropyParticlePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "name of the file describing the system ground state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-na", "minimum size of the particles whose entropy has to be evaluated", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-na", "maximum size of the particles whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new BooleanOption  ('\n', "compute-lvalue", "compute the L value of each reduced density matrix eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n', "compute-svalue", "compute the S value of each reduced density matrix eigenstate");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the reduced density matrices in the a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "lza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Lz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "sza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Sz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
  (*ToolsGroup) += new SingleDoubleOption  ('\n', "diag-precision", "convergence precision in non LAPACK mode", 1e-7);
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpinEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterLza = Manager.GetInteger("lza-eigenstate");
  int FilterSza = Manager.GetInteger("sza-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");

  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  bool ComputeLValueFlag = Manager.GetBoolean("compute-lvalue");
  bool ComputeSValueFlag = Manager.GetBoolean("compute-svalue");
  char* FileName = Manager.GetString("input-file");
  bool SVDFlag = Manager.GetBoolean("use-svd");
  
  if (FileName == 0)
    {
      cout << " an input file has to be provided" << endl;
      return -1;
    }

  int NbrParticles=0;
  int LzMax=0;
  int TotalLz=0;
  int TotalSz=0;
  bool Statistics = true;

  if ((ComputeLValueFlag == true) && (DensityMatrixFileName == 0))
    {
      cout << "compute-lvalue only valid when density-matrix is activated" << endl;
      return - 1;
    }
  if ((ComputeSValueFlag == true) && (DensityMatrixFileName == 0))
    {
      cout << "compute-svalue only valid when density-matrix is activated" << endl;
      return - 1;
    }


  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(FileName, NbrParticles, LzMax, TotalLz, TotalSz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << FileName << endl;
      return -1;
    }
  int NbrParticlesUp = (NbrParticles + TotalSz) >> 1;
  int NbrParticlesDown = (NbrParticles - TotalSz) >> 1;
  ParticleOnSphereWithSpin* Space = 0;
  
#ifdef __64_BITS__
  if (LzMax <= 31)
#else
    if (LzMax <= 15)
#endif
      {
	Space = new FermionOnSphereWithSpin  (NbrParticles, TotalLz, LzMax, TotalSz);
      }
    else
      {
	cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	return 0;
      }	
  
  
  RealVector GroundState;
  if (GroundState.ReadVector (FileName) == false)
    {
      cout << "can't open vector file " << FileName << endl;
      return -1;      
    }   
  
  if (Space->GetHilbertSpaceDimension() != GroundState.GetVectorDimension())
    {
      cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
      return -1;
    }



  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "# N    Lz    Sz    lambda";
      if (ComputeLValueFlag == true)
	DensityMatrixFile << " L^2 L";
      if (ComputeSValueFlag == true)
	DensityMatrixFile << " S^2 S";
      DensityMatrixFile.close();
    }
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName;
      TmpFileName = ReplaceExtensionToFileName(FileName, "vec", "partent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << FileName << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);

  int MaxSubsystemNbrParticles = (NbrParticles >> 1) + (NbrParticles & 1);
  cout <<"MaxSubsystemNbrParticles = "<< MaxSubsystemNbrParticles << endl;
  if (Manager.GetInteger("max-na") > 0)
    MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
  int SubsystemNbrParticles = Manager.GetInteger("min-na");
  cout << "SubsystemNbrParticles = "<< SubsystemNbrParticles <<endl;
  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
       double EntanglementEntropy = 0.0;
       double DensitySum = 0.0;
       
       int SubsystemTotalSz = 0;
       int SubsystemMaxTotalSz = SubsystemNbrParticles;
       SubsystemTotalSz = -SubsystemNbrParticles; 
       cout <<"Loop On NA"<<endl;
       for (; SubsystemTotalSz <= SubsystemMaxTotalSz; SubsystemTotalSz += 2)
	 {
	   cout <<"Loop On SA"<<endl;
	   int SubsystemTotalLz = 0;
	   int SubsystemNbrParticlesUp = (SubsystemNbrParticles + SubsystemTotalSz) >> 1;
	   int SubsystemNbrParticlesDown = (SubsystemNbrParticles - SubsystemTotalSz) >> 1;
	   int ComplementarySubsystemNbrParticlesUp = NbrParticlesUp - SubsystemNbrParticlesUp;
	   int ComplementarySubsystemNbrParticlesDown = NbrParticlesDown - SubsystemNbrParticlesDown;
	   if ((SubsystemNbrParticlesUp >= 0) && (SubsystemNbrParticlesDown >= 0) &&
	       (ComplementarySubsystemNbrParticlesUp >= 0) && (ComplementarySubsystemNbrParticlesDown >= 0))
	     {
	       int SubsystemMaxTotalLz = (((SubsystemNbrParticlesUp * ((LzMax << 1) - SubsystemNbrParticlesUp + 1))
					   + (SubsystemNbrParticlesDown * ((LzMax << 1) - SubsystemNbrParticlesDown + 1)))) >> 1;
	       int ComplementarySubsystemMinTotalLz = ((ComplementarySubsystemNbrParticlesUp * (ComplementarySubsystemNbrParticlesUp - 1))
						       + (ComplementarySubsystemNbrParticlesDown * (ComplementarySubsystemNbrParticlesDown - 1))) >> 1;
	       int ComplementarySubsystemMaxTotalLz = (LzMax * (ComplementarySubsystemNbrParticlesUp + ComplementarySubsystemNbrParticlesDown)) - ComplementarySubsystemMinTotalLz;
	       ComplementarySubsystemMinTotalLz += (ComplementarySubsystemNbrParticlesUp + ComplementarySubsystemNbrParticlesDown) * LzMax;
	       ComplementarySubsystemMaxTotalLz += (ComplementarySubsystemNbrParticlesUp + ComplementarySubsystemNbrParticlesDown) * LzMax;
	       SubsystemTotalLz = ((SubsystemNbrParticlesUp * (SubsystemNbrParticlesUp - 1))
				   + (SubsystemNbrParticlesDown * (SubsystemNbrParticlesDown - 1))) >> 1; 
	       int ShiftedTotalLz = (TotalLz + (NbrParticles * LzMax)) >> 1;
	       for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz++)
		 {
		   cout <<"Loop On LzA"<<endl;
		   int SubsystemTrueTotalLz = ((SubsystemTotalLz << 1) - (SubsystemNbrParticles * LzMax));

 		      if (((ShiftedTotalLz - SubsystemTotalLz) <= ComplementarySubsystemMaxTotalLz) &&
 			  ((ShiftedTotalLz - SubsystemTotalLz) >= ComplementarySubsystemMinTotalLz) && 
 			  ((EigenstateFlag == false) || ((FilterLza == SubsystemTrueTotalLz) && (FilterSza == SubsystemTotalSz))))
 			{
 			  cout << "processing subsystem size=" << SubsystemSize << "  subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTrueTotalLz << " subsystem total Sz=" << SubsystemTotalSz << endl;
 			  RealSymmetricMatrix PartialDensityMatrix = Space->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTrueTotalLz, SubsystemTotalSz, GroundState);
 			  if (PartialDensityMatrix.GetNbrRow() > 1)
 			    {
 			      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
 #ifdef __LAPACK__
 			      if (LapackFlag == true)
 				{
 				  if ((EigenstateFlag == true) && (FilterLza == SubsystemTrueTotalLz) && (FilterSza == SubsystemTotalSz))
 				    {
 				      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
 								PartialDensityMatrix.GetNbrRow(), true);
 				      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
 					TmpEigenstates[i][i] = 1.0;
 				      PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
 				      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
 				      char* TmpEigenstateName = new char[512];
 				      int MaxNbrEigenstates = NbrEigenstates;
 				      if (NbrEigenstates == 0)
 					MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
 				      for (int i = 0; i < MaxNbrEigenstates; ++i)
 					{
 					  if (TmpDiag[i] > 1e-14)
 					    {
 					      sprintf (TmpEigenstateName,
 						       "fermions_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d_sza_%d_.%d.vec",
 						       NbrParticles, LzMax, TotalLz, SubsystemSize,
 						       SubsystemNbrParticles, SubsystemTrueTotalLz, SubsystemTotalSz, i);
 					      TmpEigenstates[i].WriteVector(TmpEigenstateName);
 					    }
 					}
 				      delete[] TmpEigenstateName;
				    }
 				  else
				    {
				      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
 				      TmpDiag.SortMatrixDownOrder();
 				    }
 				}
 			      else
 				{
 				  if ((EigenstateFlag == true) && (FilterLza == SubsystemTrueTotalLz ) && (FilterSza == SubsystemTotalSz))
 				    {
 				      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
 								PartialDensityMatrix.GetNbrRow(), true);
 				      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
 					TmpEigenstates[i][i] = 1.0;
 				      PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
 				      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
 				      char* TmpEigenstateName = new char[512];
 				      int MaxNbrEigenstates = NbrEigenstates;
 				      if (NbrEigenstates == 0)
 					MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
 				      for (int i = 0; i < MaxNbrEigenstates; ++i)
 					{
 					  if (TmpDiag[i] > 1e-14)
 					    {
					      sprintf (TmpEigenstateName,
 						       "fermions_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d_sza_%d.%d.vec",
 						       NbrParticles, LzMax, TotalLz, SubsystemSize,
 						       SubsystemNbrParticles, SubsystemTrueTotalLz, SubsystemTotalSz, i);
 					      TmpEigenstates[i].WriteVector(TmpEigenstateName);
 					    }
 					}
 				      delete[] TmpEigenstateName;
 				    }
 				  else
 				    {
 				      PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
 				      TmpDiag.SortMatrixDownOrder();
 				    }
 				}
#else
 			      if ((EigenstateFlag == true) && (FilterLza == SubsystemTrueTotalLz) && (FilterSza == SubsystemTotalSz))
 				{
 				  if (PartialDensityMatrix.GetNbrRow() == 1)
 				    {
 				      PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
 				      TmpDiag.SortMatrixDownOrder();
 				    }
 				  else
 				    {
 				      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
 								PartialDensityMatrix.GetNbrRow(), true);
 				      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
 					TmpEigenstates[i][i] = 1.0;
 				      PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
 				      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
 				      char* TmpEigenstateName = new char[512];
 				      int MaxNbrEigenstates = NbrEigenstates;
 				      if (NbrEigenstates == 0)
 					MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
 				      for (int i = 0; i < MaxNbrEigenstates; ++i)
 					{
 					  if (TmpDiag[i] > 1e-14)
 					    {
 					      sprintf (TmpEigenstateName,
 						       "fermions_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d_sza_%d.%d.vec",
 						       NbrParticles, LzMax, TotalLz, SubsystemSize,
 						       SubsystemNbrParticles, SubsystemTrueTotalLz, SubsystemTotalSz, i);
 					      TmpEigenstates[i].WriteVector(TmpEigenstateName);
 					    }
 					}
 				      delete[] TmpEigenstateName;
 				    }
 				}
 			      else
 				{
 				  PartialDensityMatrix.Diagonalize(TmpDiag, Manager.GetDouble("diag-precision"));
 				  TmpDiag.SortMatrixDownOrder();
 				}
#endif		  
 			      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
 				{
 				  if (TmpDiag[i] > 1e-14)
				    {
 				      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
 				      DensitySum += TmpDiag[i];
 				    }
 				}
 			      if (DensityMatrixFileName != 0)
				{
 				  ofstream DensityMatrixFile;
 				  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
 				  DensityMatrixFile.precision(14);
 				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
 				    DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTrueTotalLz << " " << SubsystemTotalSz << " " << TmpDiag[i] << endl;
 				  DensityMatrixFile.close();
 				}
 			    }
 			  else
 			    if (PartialDensityMatrix.GetNbrRow() == 1)
 			      {
 				double TmpValue = PartialDensityMatrix(0,0);
 				if (TmpValue > 1e-14)
 				  {
 				    EntanglementEntropy += TmpValue * log(TmpValue);
 				    DensitySum += TmpValue;
 				  }
 				if (DensityMatrixFileName != 0)
 				  {
 				    ofstream DensityMatrixFile;
 				    DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
 				    DensityMatrixFile.precision(14);
 				    DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTrueTotalLz << " " << SubsystemTotalSz << " " << TmpValue << endl;
 				    DensityMatrixFile.close();
 				  }		  
 			      }
			}
		 }
	     }
	 }
       File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << endl;
    }
  return 0;
}
