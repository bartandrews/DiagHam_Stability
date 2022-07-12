#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/Endian.h"


#include <iostream>
#include <fstream>

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
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += LanczosGroup;
  Manager += MiscGroup;

  (*LanczosGroup) += new SingleDoubleOption ('p', "lanczos-precision", "define Lanczos precision for eigenvalues (0 if automatically defined by the program)", 0);
  (*LanczosGroup) += new SingleDoubleOption ('s', "lanczos-shift", "define the shift applied to the hamiltonian during the original Lanczos process", 0);
  (*LanczosGroup) += new BooleanOption  ('e', "eigenstate", "compute the ground state", false);  
  (*LanczosGroup) += new  SingleStringOption ('o', "ground-filename", "name of the file where the ground state has to be stored", "ground.vec");
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-iter", "set a new number of lanczos iteration (0 if the one of the lanczos.dat has to be kept)", 0);
  (*LanczosGroup) += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type ReplayFastLanczos -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  double Shift = ((SingleDoubleOption*) Manager["lanczos-shift"])->GetDouble();
  bool EigenstateFlag = ((BooleanOption*) Manager["eigenstate"])->GetBoolean();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  bool BlockLanczosFlag = ((BooleanOption*) Manager["block-lanczos"])->GetBoolean();

  int LanczosIndex;
  double PreviousLastWantedEigenvalue;
  int TmpDimension;
  RealTriDiagonalSymmetricMatrix TridiagonalizedMatrix(4000, true);

  if (BlockLanczosFlag == false)
    {
      ifstream File;
      File.open("lanczos.dat", ios::binary | ios::in);
      ReadLittleEndian(File, LanczosIndex);
      ReadLittleEndian(File, PreviousLastWantedEigenvalue);
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
	  GroundState.WriteVector(((SingleStringOption*) Manager["ground-filename"])->GetString());
	}
    }
  else
    {
      ifstream File;
      File.open("lanczos.dat", ios::binary | ios::in);
      int NbrEigenvalue = 0;
      int BlockSize = 0;
      int MaximumNumberIteration = 0;
      ReadLittleEndian(File, LanczosIndex);
      ReadLittleEndian(File, NbrEigenvalue);      
      ReadLittleEndian(File, BlockSize);            
      ReadLittleEndian(File, PreviousLastWantedEigenvalue);
      ReadLittleEndian(File, MaximumNumberIteration);
      ReadLittleEndian(File, TmpDimension);
//       TridiagonalizedMatrix.Resize(TmpDimension, TmpDimension);
//       --TmpDimension;
//       int TwiceBlockSize = 2 * this->BlockSize;
//       int TmpMax = TmpDimension - TwiceBlockSize;
//       for (int i = 0; i < TmpMax; ++i)    
// 	{    
// 	  for (int j = 0; j < TwiceBlockSize; ++j)
// 	    ReadLittleEndian(File, this->TridiagonalizedMatrix(i, i + j));
// 	}  
//       for (int i = TmpMax; i < TmpDimension; ++i)    
// 	{    
// 	  for (int j = TmpMax + 1; j < TmpDimension; ++j)
// 	    ReadLittleEndian(File, this->TridiagonalizedMatrix(i, j));
// 	}  
      File.close();  
    }
  return 0;  
}
