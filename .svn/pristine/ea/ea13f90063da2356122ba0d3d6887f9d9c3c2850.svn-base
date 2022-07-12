#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealBandDiagonalSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"

#include "GeneralTools/Endian.h"


#include <iostream>
#include <fstream>
#include <cstring>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("ReplayFastLanczos" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += LanczosGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*LanczosGroup) += new SingleDoubleOption ('p', "lanczos-precision", "define Lanczos precision for eigenvalues (0 if automatically defined by the program)", 0);  
  (*LanczosGroup) += new SingleDoubleOption ('s', "lanczos-shift", "define the shift applied to the hamiltonian during the original Lanczos process", 0);
  (*LanczosGroup) += new BooleanOption  ('e', "eigenstate", "compute the ground state", false);  
  (*LanczosGroup) += new  SingleStringOption ('o', "ground-filename", "name of the file where the ground state has to be stored (in default ground.vec.XX, with XX=nbr iteration)", 0);
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-iter", "set a new number of lanczos iteration (0 if the one of the lanczos.dat has to be kept)", 0);
  (*LanczosGroup) += new BooleanOption ('c', "complex-lanczos", "indicate whether a complex Lanczos algorithm was used");
  (*LanczosGroup) += new BooleanOption  ('\n', "projector-lanczos", "projectors were used in Lanczos algorithm", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type ReplayFastLanczos -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  double Shift = ((SingleDoubleOption*) Manager["lanczos-shift"])->GetDouble();
  bool EigenstateFlag = Manager.GetBoolean("eigenstate");
  int NbrIter = Manager.GetInteger("nbr-iter");
  bool ProjectorFlag = Manager.GetBoolean("projector-lanczos");  
  bool BlockLanczosFlag = Manager.GetBoolean("block-lanczos");  
  bool LapackFlag = Manager.GetBoolean("use-lapack");  
  bool AutomaticProjectorConstructionFlag = false;
  int NbrProjectors;
  int LanczosIndex;
  int InitialNbrProjectors;
  double *ProjectorEigenvalues = NULL;
  double PreviousLastWantedEigenvalue;
  int TmpDimension;
  RealTriDiagonalSymmetricMatrix TridiagonalizedMatrix(4000, true);

  char *OutputName;

  if (Manager.GetString("ground-filename") != NULL)
    {
      OutputName = new char[strlen(Manager.GetString("ground-filename"))+16];
      sprintf(OutputName,"%s",Manager.GetString("ground-filename"));
    }
  else
    {
      ifstream File;
      int LanczosIndex;
      File.open("lanczos.dat", ios::binary | ios::in);
      if (!File.is_open())
	{
	  cout << "No disk storage found: Could not open file lanczos.dat"<<endl;
	  exit(-1);
	}
      ReadLittleEndian(File, LanczosIndex);
      File.close();
      OutputName = new char[30];
      sprintf(OutputName, "ground.vec.%d", LanczosIndex);
    }
  
  if (!Manager.GetBoolean("complex-lanczos"))
    {
      if (BlockLanczosFlag == false)
	{
	  ifstream File;
	  File.open("lanczos.dat", ios::binary | ios::in);
	  if (!File.is_open())
	    {
	      cout << "No disk storage found: Could not open file lanczos.dat"<<endl;
	      exit(-1);
	    }
	  ReadLittleEndian(File, LanczosIndex);
	  ReadLittleEndian(File, PreviousLastWantedEigenvalue);
	  if (ProjectorFlag == true)
	    {
	      ReadLittleEndian(File, AutomaticProjectorConstructionFlag);
	      ReadLittleEndian(File, NbrProjectors);
	      ReadLittleEndian(File, InitialNbrProjectors);
	      if (ProjectorEigenvalues != 0)
		{
		  for (int i = 0; i < NbrProjectors; ++i)
		    ReadLittleEndian(File, ProjectorEigenvalues[i]);
		}
	    }

	  ReadLittleEndian(File, TmpDimension);
	  TridiagonalizedMatrix.Resize(TmpDimension, TmpDimension);
	  --TmpDimension;
	  for (int i = 0; i <= TmpDimension; ++i)
	    {
	      ReadLittleEndian(File, TridiagonalizedMatrix.DiagonalElement(i));
	    }
	  for (int i = 0; i < TmpDimension; ++i)
	    {
	      ReadLittleEndian(File, TridiagonalizedMatrix.UpperDiagonalElement(i));
	    }
	  cout << "Previous last wanted eigenvalue is: " << (PreviousLastWantedEigenvalue - Shift) << endl;
	  cout << "Lanczos index is: " << LanczosIndex << endl;
	  cout << "current dimension of the tridiagonal matrix: " << TmpDimension << endl;
      
	  double PreviousGroundStateEnergy = 0.0;
	  for (int i = 4; i <= TmpDimension; ++i)
	    {
	      RealTriDiagonalSymmetricMatrix DiagonalizedMatrix(i, true);      
	      TridiagonalizedMatrix.Resize(i, i);
	      int Dimension = TridiagonalizedMatrix.GetNbrRow();
	      DiagonalizedMatrix.Copy(TridiagonalizedMatrix);
	      DiagonalizedMatrix.Diagonalize(50);
	      double GroundStateEnergy = DiagonalizedMatrix.DiagonalElement(0);      
	      for (int DiagPos = 1; DiagPos < Dimension; DiagPos++)
		if (DiagonalizedMatrix.DiagonalElement(DiagPos) < GroundStateEnergy)
		  GroundStateEnergy = DiagonalizedMatrix.DiagonalElement(DiagPos);  
	      GroundStateEnergy -= Shift;
	      cout << GroundStateEnergy << " " << GroundStateEnergy << " " << ((PreviousGroundStateEnergy - GroundStateEnergy) / GroundStateEnergy) << " (step " << (i + 1) << ")" << endl;
	      PreviousGroundStateEnergy = GroundStateEnergy;
	    }
      
	  if ((NbrIter > 0) && (NbrIter <= TmpDimension))
	    {
	      TridiagonalizedMatrix.Resize(NbrIter - 1, NbrIter - 1);      
	    }
      
	  if (EigenstateFlag == true)
	    {
	      RealMatrix TmpEigenvector (TridiagonalizedMatrix.GetNbrRow(), TridiagonalizedMatrix.GetNbrRow(), true);
	      for (int i = 0; i < TridiagonalizedMatrix.GetNbrRow(); ++i)
		TmpEigenvector(i, i) = 1.0;
	      RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (TridiagonalizedMatrix.GetNbrRow());
	      SortedDiagonalizedMatrix.Copy(TridiagonalizedMatrix);
	      SortedDiagonalizedMatrix.Diagonalize(TmpEigenvector);
	      SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);
	      double* TmpComponents = new double [TridiagonalizedMatrix.GetNbrRow()];
	      for (int j = 0; j < TridiagonalizedMatrix.GetNbrRow(); ++j)
		{
		  TmpComponents[j] = TmpEigenvector(j, 0);
		}
	  
	      RealVector V1;
	      RealVector GroundState;
	      V1.ReadVector("vector.0");	      
	      GroundState.Copy(V1, TmpComponents[0]);
	      char* TmpVectorName = new char [256];
	      for (int i = 1; i < TridiagonalizedMatrix.GetNbrRow(); ++i)
		{
		  sprintf(TmpVectorName, "vector.%d", i);
		  V1.ReadVector(TmpVectorName);	      
		  GroundState.AddLinearCombination(TmpComponents[i], V1);	      
		  cout << ".";
		  cout.flush();
		}	  
	      delete[] TmpVectorName;      
	      cout << endl;
	      GroundState /= GroundState.Norm();
	      GroundState.WriteVector(OutputName);
	    }
	}
      else
	{
	  ifstream File;
	  File.open("lanczos.dat", ios::binary | ios::in);
	  if (!File.is_open())
	    {
	      cout << "No disk storage found: Could not open file lanczos.dat"<<endl;
	      exit(-1);
	    }	  
	  int NbrEigenvalue = 0;
	  int BlockSize = 0;
	  int MaximumNumberIteration = 0;
	  ReadLittleEndian(File, LanczosIndex);
	  ReadLittleEndian(File, NbrEigenvalue);      
	  ReadLittleEndian(File, BlockSize);            
	  ReadLittleEndian(File, PreviousLastWantedEigenvalue);
	  ReadLittleEndian(File, MaximumNumberIteration);
	  ReadLittleEndian(File, TmpDimension);

	  cout << "number of wanted eigenvalues : " << NbrEigenvalue << endl;
	  cout << "block size : " << BlockSize << endl;
	  cout << "number of iterations : " << LanczosIndex << endl;
	  cout << "matrix size : " << TmpDimension << endl;
	  int TwiceBlockSize = 2 * BlockSize;
	  int TmpMax = TmpDimension - TwiceBlockSize;

	  RealBandDiagonalSymmetricMatrix ReducedMatrix(TmpDimension, BlockSize, true);
	  for (int i = 0; i < TmpMax; ++i)    
	    {    
	      for (int j = 0; j < TwiceBlockSize; ++j)
		ReadLittleEndian(File, ReducedMatrix(i, i + j));
	    }  
	  for (int i = TmpMax; i < TmpDimension; ++i)    
	    {    
	      for (int j = TmpMax + 1; j < TmpDimension; ++j)
		ReadLittleEndian(File, ReducedMatrix(i, j));
	    }  
	  File.close();  

	  int Dimension = ReducedMatrix.GetNbrRow();
	  if ((NbrIter > 0) && ((NbrIter * BlockSize) <Dimension ))
	    {
	      Dimension = NbrIter * BlockSize;
	      ReducedMatrix.Resize(Dimension, Dimension);
	    }
	  RealBandDiagonalSymmetricMatrix TemporaryReducedMatrix (TmpDimension, BlockSize, true);
	  RealTriDiagonalSymmetricMatrix TridiagonalizedMatrix (TmpDimension, true);
	  RealTriDiagonalSymmetricMatrix DiagonalizedMatrix (TmpDimension, true);
	  TemporaryReducedMatrix.Copy(ReducedMatrix);
	  if (EigenstateFlag == true)
	    {
	      cout << "diagonalizing matrix"  << endl;
	      RealVector* LanczosVectors = new RealVector [3 * BlockSize];
	      RealVector* Eigenstates = new RealVector [BlockSize];
	      RealMatrix TmpEigenvector (ReducedMatrix.GetNbrRow(), ReducedMatrix.GetNbrRow(), true);
	      for (int i = 0; i < ReducedMatrix.GetNbrRow(); ++i)
		TmpEigenvector(i, i) = 1.0;
	      
	      RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (ReducedMatrix.GetNbrRow());
	      TemporaryReducedMatrix.Copy(ReducedMatrix);
#ifdef __LAPACK__
	      if (LapackFlag == true)
		{
		  RealDiagonalMatrix TmpDiag (SortedDiagonalizedMatrix.GetNbrColumn());
		  TemporaryReducedMatrix.LapackDiagonalize(TmpDiag, TmpEigenvector);
		  for (int i = 0; i < SortedDiagonalizedMatrix.GetNbrColumn(); ++i)
		    SortedDiagonalizedMatrix.DiagonalElement(i) = TmpDiag[i];
		}
	      else
		{
#endif
		  TemporaryReducedMatrix.Tridiagonalize(SortedDiagonalizedMatrix, 1e-7, TmpEigenvector);
		  SortedDiagonalizedMatrix.Diagonalize(TmpEigenvector);
#ifdef __LAPACK__
		}
#endif
	      SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);

	      cout << "computing ground states"  << endl;
	      double* TmpCoefficents = new double [BlockSize];
	      char* TmpVectorName = new char [256];
	      for (int i = 0; i < BlockSize; ++i)
		{
		  for (int j = 0; j < BlockSize; ++j)
		    {
		      sprintf(TmpVectorName, "vector.%d", j);
		      LanczosVectors[j].ReadVector(TmpVectorName);
		    }
		  Eigenstates[i].Copy(LanczosVectors[0], TmpEigenvector(0, i));
		  for (int j = 1; j < BlockSize; ++j)
		    TmpCoefficents[j - 1] = TmpEigenvector(j, i);	  
		  AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), &(LanczosVectors[1]), BlockSize - 1,  TmpCoefficents);
		  Operation.ApplyOperation(Architecture.GetArchitecture());
		}       
	      for (int i = 1; i < NbrIter; ++i)
		{
		  for (int j = 0; j < BlockSize; ++j)
		    {
		      sprintf(TmpVectorName, "vector.%d", ((i * BlockSize) + j));
		      LanczosVectors[j].ReadVector(TmpVectorName);
		    }
		  for (int k = 0; k < BlockSize; ++k)
		    {
		      for (int j = 0; j < BlockSize; ++j)
			TmpCoefficents[j] = TmpEigenvector((i * BlockSize) + j, k);	  
		      AddRealLinearCombinationOperation Operation (&(Eigenstates[k]), LanczosVectors, BlockSize,  TmpCoefficents);
		      Operation.ApplyOperation(Architecture.GetArchitecture());
		    }
		  cout << i << "/" << LanczosIndex << "           \r";
		  cout.flush();
		}
	      cout << endl;
	    
	      for (int i = 0; i < BlockSize; ++i)
		{
		  Eigenstates[i] /= Eigenstates[i].Norm();
		  sprintf(TmpVectorName, "groundvector.%d.vec", i);
		  Eigenstates[i].WriteVector(TmpVectorName);
		  cout << SortedDiagonalizedMatrix.DiagonalElement(i) << endl;
		}
	      delete[] TmpVectorName;
	      delete[] TmpCoefficents;
	    }
	  else
	    {
#ifdef __LAPACK__
	      if (LapackFlag == true)
		{
		  RealDiagonalMatrix TmpDiag (TemporaryReducedMatrix.GetNbrColumn());
		  TemporaryReducedMatrix.LapackDiagonalize(TmpDiag);
		  DiagonalizedMatrix.Resize(TemporaryReducedMatrix.GetNbrColumn(), TemporaryReducedMatrix.GetNbrColumn());
		  for (int i = 0; i < TemporaryReducedMatrix.GetNbrColumn(); ++i)
		    DiagonalizedMatrix.DiagonalElement(i) = TmpDiag[i];
		}
	      else
		{
#endif
		  TemporaryReducedMatrix.Tridiagonalize(DiagonalizedMatrix, 1e-7);
		  DiagonalizedMatrix.Diagonalize();
#ifdef __LAPACK__
		}
#endif
	      DiagonalizedMatrix.SortMatrixUpOrder();
	      for (int i = 0; i < NbrEigenvalue; ++i)
		cout << DiagonalizedMatrix.DiagonalElement(i) << endl;
	    }
	}
    }
  else // have complex Lanczos data
    {
      ifstream File;
      File.open("lanczos.dat", ios::binary | ios::in);
      if (!File.is_open())
	{
	  cout << "No disk storage found: Could not open file lanczos.dat"<<endl;
	  exit(-1);
	}
      ReadLittleEndian(File, LanczosIndex);
      ReadLittleEndian(File, PreviousLastWantedEigenvalue);
      if (ProjectorFlag == true)
	{
	  ReadLittleEndian(File, AutomaticProjectorConstructionFlag);
	  ReadLittleEndian(File, NbrProjectors);
	  ReadLittleEndian(File, InitialNbrProjectors);
	  if (ProjectorEigenvalues != 0)
	    {
	      for (int i = 0; i < NbrProjectors; ++i)
		ReadLittleEndian(File, ProjectorEigenvalues[i]);
	    }
	}

      ReadLittleEndian(File, TmpDimension);
      TridiagonalizedMatrix.Resize(TmpDimension, TmpDimension);
      --TmpDimension;
      for (int i = 0; i <= TmpDimension; ++i)
	{
	  ReadLittleEndian(File, TridiagonalizedMatrix.DiagonalElement(i));
	}
      for (int i = 0; i < TmpDimension; ++i)
	{
	  ReadLittleEndian(File, TridiagonalizedMatrix.UpperDiagonalElement(i));
	}
      cout << "Previous last wanted eigenvalue is: " << (PreviousLastWantedEigenvalue - Shift) << endl;
      cout << "Lanczos index is: " << LanczosIndex << endl;
      cout << "current dimension of the tridiagonal matrix: " << TmpDimension << endl;
      
      double PreviousGroundStateEnergy = 0.0;
      for (int i = 4; i <= TmpDimension; ++i)
	{
	  RealTriDiagonalSymmetricMatrix DiagonalizedMatrix(i, true);      
	  TridiagonalizedMatrix.Resize(i, i);
	  int Dimension = TridiagonalizedMatrix.GetNbrRow();
	  DiagonalizedMatrix.Copy(TridiagonalizedMatrix);
	  DiagonalizedMatrix.Diagonalize(50);
	  double GroundStateEnergy = DiagonalizedMatrix.DiagonalElement(0);      
	  for (int DiagPos = 1; DiagPos < Dimension; DiagPos++)
	    if (DiagonalizedMatrix.DiagonalElement(DiagPos) < GroundStateEnergy)
	      GroundStateEnergy = DiagonalizedMatrix.DiagonalElement(DiagPos);  
	  GroundStateEnergy -= Shift;
	  cout << GroundStateEnergy << " " << GroundStateEnergy << " " << ((PreviousGroundStateEnergy - GroundStateEnergy) / GroundStateEnergy) << " (step " << (i + 1) << ")" << endl;
	  PreviousGroundStateEnergy = GroundStateEnergy;
	}
      
      if ((NbrIter > 0) && (NbrIter <= TmpDimension))
	{
	  TridiagonalizedMatrix.Resize(NbrIter - 1, NbrIter - 1);      
	}
      
      if (EigenstateFlag == true)
	{
	  RealMatrix TmpEigenvector (TridiagonalizedMatrix.GetNbrRow(), TridiagonalizedMatrix.GetNbrRow(), true);
	  for (int i = 0; i < TridiagonalizedMatrix.GetNbrRow(); ++i)
	    TmpEigenvector(i, i) = 1.0;
	  
	  RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (TridiagonalizedMatrix.GetNbrRow());
	  SortedDiagonalizedMatrix.Copy(TridiagonalizedMatrix);
	  SortedDiagonalizedMatrix.Diagonalize(TmpEigenvector);
	  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);
	  double* TmpComponents = new double [TridiagonalizedMatrix.GetNbrRow()];
	  for (int j = 0; j < TridiagonalizedMatrix.GetNbrRow(); ++j)
	    {
	      TmpComponents[j] = TmpEigenvector(j, 0);
	    }
	  
	  ComplexVector V1;
	  ComplexVector GroundState;
	  V1.ReadVector("vector.0");	      
	  GroundState.Copy(V1, TmpComponents[0]);
	  char* TmpVectorName = new char [256];
	  for (int i = 1; i < TridiagonalizedMatrix.GetNbrRow(); ++i)
	    {
	      sprintf(TmpVectorName, "vector.%d", i);
	      V1.ReadVector(TmpVectorName);	      
	      GroundState.AddLinearCombination(TmpComponents[i], V1);	      
	      cout << ".";
	      cout.flush();
	    }	  
	  delete[] TmpVectorName;      
	  cout << endl;
	  GroundState /= GroundState.Norm();
	  GroundState.WriteVector(OutputName);
	}
    }
  delete [] OutputName;
  return 0;  
}
    
