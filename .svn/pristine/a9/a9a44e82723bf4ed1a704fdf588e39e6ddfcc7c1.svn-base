#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/FermionOnSphereMPSWrapper.h"
#include "HilbertSpace/FermionOnCylinderMPSWrapper.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSLaughlinMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Matrix/SparseRealMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "GeneralTools/ArrayTools.h"

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


int main(int argc, char** argv)
{
  cout.precision(14); 
  
  OptionManager Manager ("FQHESphereMPSEntanglementSpectrum" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager;

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new BooleanOption  ('\n', "use-padding", "root partitions use the extra zero padding");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "la", "number of orbitals in subsystem A", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "na", "number of particles in subsystem A", 0);
  (*SystemGroup) += new BooleanOption ('\n', "all-na", "print all charge sectors");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMPSEntanglementSpectrum -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0; 
  int NbrFluxQuanta = 0;
  int TotalLz = 0;
  int* ReferenceState = 0;
  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
    return -1;

  int EntCut = Manager.GetInteger("la");
  int Na = Manager.GetInteger("na");

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");
  double AspectRatio = Manager.GetDouble("aspect-ratio");
  double kappa = 0.0;
  if (CylinderFlag)
    {
       kappa = (2.0 * M_PI)/sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) * AspectRatio);
       cout << "Cylinder geometry, kappa= "<<kappa<<endl;
    }

  int LandauLevel = 0;

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int NbrBMatrices = 2;
  SparseRealMatrix* BMatrices = MPSMatrix->GetMatrices();
  SparseRealMatrix* ConjugateBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    ConjugateBMatrices[i] = BMatrices[i].Transpose();

  cout << "B matrix size = " << BMatrices[0].GetNbrRow() << "x" << BMatrices[0].GetNbrColumn() << endl;
  
  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;
  if (Manager.GetBoolean("use-padding") == true)
    {
      if (Manager.GetBoolean("k-2") == true)
	{
	  if ((Manager.GetInteger("r-index") & 1) == 0)
	    MPSRowIndex = Manager.GetInteger("p-truncation") + (Manager.GetInteger("r-index") / 2);
	  else
	    MPSRowIndex = 2 * Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index") - 1;
	}
      else
	{
	  if (Manager.GetBoolean("rr-3") == true)
	    {
	      MPSRowIndex = 3 * (Manager.GetInteger("p-truncation") + 1);
	    }
	  else
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + ((Manager.GetInteger("laughlin-index") - 1) / 2);
	    }
	}
      MPSColumnIndex = MPSRowIndex;
    }
  else
    {
      if (Manager.GetBoolean("k-2") == true)
	{
	  if ((Manager.GetInteger("r-index") & 1) == 0)
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index");
	      MPSColumnIndex = Manager.GetInteger("p-truncation");
	    }
	  else
	    {
	      MPSRowIndex = 2 * (Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index"));
	      MPSColumnIndex = 2 * Manager.GetInteger("p-truncation");
	    }
	}
      else
	{
	  if (Manager.GetBoolean("rr-3") == true)
	    {
	      MPSRowIndex = 3 * (Manager.GetInteger("p-truncation") + 2);
	      MPSColumnIndex = 3 * Manager.GetInteger("p-truncation");
	    }
	  else
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + (Manager.GetInteger("laughlin-index") - 1);
	      MPSColumnIndex = Manager.GetInteger("p-truncation");
	    }
	}
    }

  ofstream File;
  File.precision(14);
  
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = new char [512];
      char* StateName = new char [256];
      if (Manager.GetBoolean("k-2") == true)
	{
	  sprintf (StateName, "clustered_k_2_r_%ld", Manager.GetInteger("r-index"));
	}
      else
	{
	  if (Manager.GetBoolean("rr-3") == true)
	    {
	      sprintf (StateName, "readrezayi3");
	    }
	  else
	    {
	      sprintf (StateName, "laughlin%ld", Manager.GetInteger("laughlin-index"));
	    }
	}      
      if (CylinderFlag == true)
	{
	  sprintf(TmpFileName, "fermions_cylinder_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.ent", StateName,
		  Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	}
      else
	{
	  sprintf(TmpFileName, "fermions_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.ent", StateName,
		  Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	}
      File.open(TmpFileName, ios::binary | ios::out);     
   }

  //cout<<"B0"<<endl;
  //cout<<BMatrices[0]<<endl;
  //cout<<"B1"<<endl;
  //cout<<BMatrices[1]<<endl;

  int MatDim = BMatrices[0].GetNbrRow();
  int LambdaMax = Manager.GetInteger("p-truncation");
  int LaughlinIndex = Manager.GetInteger("laughlin-index");

  cout << "B matrix size = " << MatDim << "x" << MatDim << endl;

  double CutOff = 1e-14;
  double TmpTrace;

  double* NormalizationCoefficients = new double[NbrFluxQuanta + 1];
  BinomialCoefficients Binomial(NbrFluxQuanta);
  for (int i = 0; i <= NbrFluxQuanta; ++i)
    {      
      NormalizationCoefficients[i] = ((double) (NbrFluxQuanta + 1)) / Binomial.GetNumericalCoefficient(NbrFluxQuanta, i);
    }
  if (CylinderFlag == false)
   {
     cout<<"Sphere geometry, normalization factors: "<<endl;
     for (int i = 0; i <= NbrFluxQuanta; ++i)
       cout << NormalizationCoefficients[i] << " ";
     cout << endl;
   }    
  


  RealMatrix* TmpBUD;
  RealMatrix* UMatrix;
  RealMatrix* VMatrix;
  RealMatrix* NewDMatrix;
  RealMatrix* NewDVMatrix;
  RealMatrix* NewUDMatrix;
  double* SingularValues; 
  int MaxDim;

  RealMatrix DenseB0 (BMatrices[0]);
  RealMatrix DenseB1 (BMatrices[1]);


  cout<<"----------Start SVD from right end-------------"<<endl;

  TmpBUD = new RealMatrix (MatDim, 2, true);
  double Tmp;
  for (int i = 0; i < MatDim; i++)
    {
      BMatrices[0].GetMatrixElement(i, MPSColumnIndex, Tmp);
      TmpBUD->SetMatrixElement(i, 0, Tmp);
      BMatrices[1].GetMatrixElement(i, MPSColumnIndex, Tmp);
      TmpBUD->SetMatrixElement(i, 1, Tmp);
    }
  //cout<<"TmpBUD"<<endl;
  //cout<<*TmpBUD<<endl;

/*
  RealMatrix TestSVD (10, 2, true);
  TestSVD.SetMatrixElement(0, 0, 0.76966541249324);
  TestSVD.SetMatrixElement(3, 1, 0.87730576909835);

  cout<<"TestSVD "<<TestSVD<<endl;
  RealMatrix* TestU = new RealMatrix(10,10,true);

  RealMatrix* TestV = new RealMatrix(2,2,true);
  double* SingTest = TestSVD.SingularValueDecomposition(*TestU, *TestV, false);

  RealMatrix* TestD = new RealMatrix (10, 2, true);
  for (int i = 0; i < 2; i++)
    TestD->SetMatrixElement(i,i,SingTest[i]);  

  cout<<"D "<<*TestD<<endl;

  RealMatrix TestUDV (10, 2, true);
  TestUDV.Copy(*TestU);
  
  TestUDV.Multiply(*TestD);
  TestUDV.Multiply(*TestV);
  cout<<"Test"<<endl;
  cout<<TestUDV<<endl;
  exit(1);
*/

  //cout<<"U"<<endl;
  //cout<<(*TestU)<<endl;
  //cout<<"V"<<endl;
  //cout<<(*TestV)<<endl;
  

  UMatrix = new RealMatrix (MatDim, MatDim, true);
  VMatrix = new RealMatrix (2, 2, true);

  SingularValues = TmpBUD->SingularValueDecomposition(*UMatrix, *VMatrix, false);

  int SingularDimension = 2;
  if (MatDim < 2)
    SingularDimension = MatDim;
 
  //cout<<"Singular values: ";
  //for (int i = 0; i < SingularDimension; i++)
  //  cout<<SingularValues[i]<<" ";
  //cout<<endl;

  int NbrNonZeroSingularValues = 0;
  TmpTrace = 0.0;
  for (int i = 0; i < SingularDimension; i++)
   if (SingularValues[i] > CutOff)
    {
     NbrNonZeroSingularValues++;
     TmpTrace += SingularValues[i];
    }

  //cout<<"non zero = "<<NbrNonZeroSingularValues<<endl;

  //cout<<"UMatrix "<<endl;
  //cout<<*UMatrix<<endl;

  //cout<<"VMatrix "<<endl;
  //cout<<*VMatrix<<endl;

  UMatrix->Resize(MatDim, NbrNonZeroSingularValues);
  VMatrix->Resize(NbrNonZeroSingularValues, 2);
  
  NewDMatrix = new RealMatrix (NbrNonZeroSingularValues, NbrNonZeroSingularValues, true);
  for (int i = 0; i < NbrNonZeroSingularValues; ++i)
    NewDMatrix->SetMatrixElement(i, i, SingularValues[i]);

  NewUDMatrix = new RealMatrix (MatDim, NbrNonZeroSingularValues, true);
  NewUDMatrix->Copy(UMatrix->Multiply(*NewDMatrix));

  //cout<<"NewUD"<<*NewUDMatrix<<endl;

  delete UMatrix;
  delete VMatrix;
  delete NewDMatrix;
  delete TmpBUD;
  delete[] SingularValues;


  int OldNbrNonZeroSingularValues = NbrNonZeroSingularValues;

  for (int TmpCut = NbrFluxQuanta; TmpCut > EntCut; TmpCut--)
    {
       cout << "Bond " << TmpCut << endl;

       RealMatrix TmpB0 (MatDim, MatDim);
       TmpB0.Copy(DenseB0);
       TmpB0.Multiply(*NewUDMatrix);
       TmpB0.Resize(MatDim, OldNbrNonZeroSingularValues);
       RealMatrix TmpB1 (MatDim, MatDim);
       TmpB1.Copy(DenseB1);
       TmpB1.Multiply(*NewUDMatrix);
       TmpB1.Resize(MatDim, OldNbrNonZeroSingularValues);

       TmpBUD = new RealMatrix (MatDim, 2 * OldNbrNonZeroSingularValues, true);

       for (int i = 0; i < MatDim; i++)
         for (int j = 0; j < OldNbrNonZeroSingularValues; ++j)
           {
             TmpB0.GetMatrixElement(i, j, Tmp);
             TmpBUD->SetMatrixElement(i, j, Tmp);
             TmpB1.GetMatrixElement(i, j, Tmp);
             TmpBUD->SetMatrixElement(i, j + OldNbrNonZeroSingularValues, Tmp);
           }

       if (fabs(TmpTrace) > CutOff)  
          (*TmpBUD) /= TmpTrace;
       else
         { 
           cout << "Warning: trying to normalize with 0! " << endl;
           exit(1);
         } 

       delete NewUDMatrix;

       //cout<<"TmpBUD"<<endl;
       //cout<<*TmpBUD<<endl;


       UMatrix = new RealMatrix (MatDim, MatDim, true);
       VMatrix = new RealMatrix (2 * OldNbrNonZeroSingularValues, 2 * OldNbrNonZeroSingularValues, true);

       SingularValues = TmpBUD->SingularValueDecomposition(*UMatrix, *VMatrix, false);
 
       SingularDimension = 2 * OldNbrNonZeroSingularValues;
       if (MatDim < 2 * OldNbrNonZeroSingularValues)
          SingularDimension = MatDim;

       //cout<<"Singular values: ";
       //for (int i = 0; i < SingularDimension; i++)
       //  cout<<SingularValues[i]<<" ";
       //cout<<endl;

       int NbrNonZeroSingularValues = 0;
       TmpTrace = 0.0;
       for (int i = 0; i < SingularDimension; i++)
         if (SingularValues[i] > CutOff)
          {
            NbrNonZeroSingularValues++; 
            TmpTrace += SingularValues[i];
          }

       //cout<<"non zero = "<<NbrNonZeroSingularValues<<endl;

/*

       cout<<"UMatrix "<<endl;
       cout<<*UMatrix<<endl;

       cout<<"VMatrix "<<endl;
       cout<<*VMatrix<<endl;

        RealMatrix DTest (MatDim, 2 * OldNbrNonZeroSingularValues, true);
        for (int i = 0; i < SingularDimension; ++i)
          DTest.SetMatrixElement(i, i, SingularValues[i]);

       cout<<"DTest "<<endl;
       cout<<DTest<<endl;


        RealMatrix SingTest (MatDim, 2 * OldNbrNonZeroSingularValues, true);
        SingTest.Copy(*UMatrix);
        SingTest.Multiply(DTest);
        SingTest.Multiply(*VMatrix);
        cout<<"Test SVD"<<endl;
        cout<<SingTest<<endl;

        exit(1);

*/
        UMatrix->Resize(MatDim, NbrNonZeroSingularValues);
        VMatrix->Resize(NbrNonZeroSingularValues, 2 * OldNbrNonZeroSingularValues);
       
        NewDMatrix = new RealMatrix (NbrNonZeroSingularValues, NbrNonZeroSingularValues, true);
        for (int i = 0; i < NbrNonZeroSingularValues; ++i)
          NewDMatrix->SetMatrixElement(i, i, SingularValues[i]);
    
        NewUDMatrix = new RealMatrix (MatDim, 2 * OldNbrNonZeroSingularValues, true);
        NewUDMatrix->Copy(UMatrix->Multiply(*NewDMatrix));
        NewUDMatrix->Resize(MatDim, NbrNonZeroSingularValues);

        //cout<<"NewUDMatrix "<<*NewUDMatrix<<endl;
 
        delete UMatrix;
        delete VMatrix;
        delete TmpBUD;
        delete NewDMatrix;
        delete[] SingularValues;

        OldNbrNonZeroSingularValues = NbrNonZeroSingularValues;
       
    }

  
  RealMatrix NewUDOld (*NewUDMatrix);
  //cout<<"Carry away NewUD "<<endl;
  //cout<<NewUDOld<<endl;

  
  cout<<"----------Start SVD from left end-------------"<<endl;



  TmpBUD = new RealMatrix (2, MatDim, true);
  for (int i = 0; i < MatDim; i++)
    {
      BMatrices[0].GetMatrixElement(MPSRowIndex, i, Tmp);
      TmpBUD->SetMatrixElement(0, i, Tmp);
      BMatrices[1].GetMatrixElement(MPSRowIndex, i, Tmp);
      TmpBUD->SetMatrixElement(1, i, Tmp);
    }
  //cout<<"TmpBUD"<<endl;
  //cout<<*TmpBUD<<endl;


  UMatrix = new RealMatrix (2, 2, true);
  VMatrix = new RealMatrix (MatDim, MatDim, true);

  SingularValues = TmpBUD->SingularValueDecomposition(*UMatrix, *VMatrix, false);

  SingularDimension = 2;
  if (MatDim < 2)
    SingularDimension = MatDim;
 
  //cout<<"Singular values: ";
  //for (int i = 0; i < SingularDimension; i++)
  //  cout<<SingularValues[i]<<" ";
  //cout<<endl;

  NbrNonZeroSingularValues = 0;
  TmpTrace = 0.0;
  for (int i = 0; i < SingularDimension; i++)
   if (SingularValues[i] > CutOff)
    { 
      NbrNonZeroSingularValues++;
      TmpTrace += SingularValues[i];
    }  

  //cout<<"non zero = "<<NbrNonZeroSingularValues<<endl;

  //cout<<"UMatrix "<<endl;
  //cout<<*UMatrix<<endl;

  //cout<<"VMatrix "<<endl;
  //cout<<*VMatrix<<endl;


  UMatrix->Resize(2, NbrNonZeroSingularValues);
  VMatrix->Resize(NbrNonZeroSingularValues, MatDim);
  
  NewDMatrix = new RealMatrix (NbrNonZeroSingularValues, NbrNonZeroSingularValues, true);
  for (int i = 0; i < NbrNonZeroSingularValues; ++i)
    NewDMatrix->SetMatrixElement(i, i, SingularValues[i]);


  NewDVMatrix = new RealMatrix (NbrNonZeroSingularValues, MatDim, true);
  NewDVMatrix->Copy(NewDMatrix->Multiply(*VMatrix));

  //cout<<"NewDV"<<*NewDVMatrix<<endl;


  delete UMatrix;
  delete VMatrix;
  delete NewDMatrix;
  delete TmpBUD;
  delete[] SingularValues;

  OldNbrNonZeroSingularValues = NbrNonZeroSingularValues;

  for (int TmpCut = 2; TmpCut <= EntCut; TmpCut++)
    {
       cout << "Bond " << TmpCut << endl;

       RealMatrix TmpB0 (MatDim, MatDim);
       TmpB0.Copy(*NewDVMatrix);
       TmpB0.Multiply(DenseB0);
       TmpB0.Resize(OldNbrNonZeroSingularValues,MatDim);

       RealMatrix TmpB1 (MatDim, MatDim);
       TmpB1.Copy(*NewDVMatrix);
       TmpB1.Multiply(DenseB1);
       TmpB1.Resize(OldNbrNonZeroSingularValues, MatDim);

       TmpBUD = new RealMatrix (2 * OldNbrNonZeroSingularValues, MatDim, true);

       for (int i = 0; i < OldNbrNonZeroSingularValues; i++)
         for (int j = 0; j < MatDim; ++j)
           {
             TmpB0.GetMatrixElement(i, j, Tmp);
             TmpBUD->SetMatrixElement(i, j, Tmp);
             TmpB1.GetMatrixElement(i, j, Tmp);
             TmpBUD->SetMatrixElement(i + OldNbrNonZeroSingularValues, j, Tmp);
           }

       if (fabs(TmpTrace) > CutOff)  
          (*TmpBUD) /= TmpTrace;
       else
         { 
          cout << "Warning: trying to normalize with 0! " << endl;
          exit(1);
         } 


       delete NewDVMatrix;

       //cout<<"TmpB0"<<endl;
       //cout<<TmpB0<<endl;
       //cout<<"TmpB1"<<endl;
       //cout<<TmpB1<<endl;

       //cout<<"TmpBUD"<<endl;
       //cout<<*TmpBUD<<endl;


       UMatrix = new RealMatrix (2 * OldNbrNonZeroSingularValues, 2 * OldNbrNonZeroSingularValues, true);
       VMatrix = new RealMatrix (MatDim, MatDim, true);


       //cout<<"TmpBUD "<<TmpBUD->GetNbrRow()<<" "<<TmpBUD->GetNbrColumn()<<" ; "<<(2 * OldNbrNonZeroSingularValues)<<" "<<MatDim<<endl;

       SingularValues = TmpBUD->SingularValueDecomposition(*UMatrix, *VMatrix, false);
 
       SingularDimension = 2 * OldNbrNonZeroSingularValues;
       if (MatDim < 2 * OldNbrNonZeroSingularValues)
          SingularDimension = MatDim;

       //cout<<"Singular values: ";
       //for (int i = 0; i < SingularDimension; i++)
       //  cout<<SingularValues[i]<<" ";
       //cout<<endl;

       int NbrNonZeroSingularValues = 0;
       TmpTrace = 0.0;
       for (int i = 0; i < SingularDimension; i++)
         if (SingularValues[i] > CutOff)
          { 
            NbrNonZeroSingularValues++; 
            TmpTrace += SingularValues[i];
          }

       //cout<<"non zero = "<<NbrNonZeroSingularValues<<endl;

       UMatrix->Resize(2 * OldNbrNonZeroSingularValues, NbrNonZeroSingularValues);
       VMatrix->Resize(NbrNonZeroSingularValues, MatDim);

       //cout<<"UMatrix "<<endl;
       //cout<<*UMatrix<<endl;

       //cout<<"VMatrix "<<endl;
       //cout<<*VMatrix<<endl;
       

       NewDMatrix = new RealMatrix (NbrNonZeroSingularValues, NbrNonZeroSingularValues, true);
       for (int i = 0; i < NbrNonZeroSingularValues; ++i)
          NewDMatrix->SetMatrixElement(i, i, SingularValues[i]);

       //cout<<"DMatrix "<<endl;
       //cout<<*NewDMatrix<<endl;

       NewDVMatrix = new RealMatrix (NbrNonZeroSingularValues, MatDim, true);
       NewDVMatrix->Copy(NewDMatrix->Multiply(*VMatrix));

       //cout<<"NewDV"<<*NewDVMatrix<<endl;

        delete UMatrix;
        delete VMatrix;
        delete TmpBUD;
        delete NewDMatrix;
        delete[] SingularValues;

        OldNbrNonZeroSingularValues = NbrNonZeroSingularValues;    
 
    }

  cout<<"Completed..........................."<<endl;

  RealMatrix EntMatrixOnSite (MatDim, MatDim, true);
  EntMatrixOnSite.Copy(*NewDVMatrix);
  EntMatrixOnSite.Multiply(NewUDOld);
  EntMatrixOnSite.Resize(NewDVMatrix->GetNbrRow(), NewUDOld.GetNbrColumn());

  UMatrix = new RealMatrix (NewDVMatrix->GetNbrRow(), NewDVMatrix->GetNbrRow(), true);
  VMatrix = new RealMatrix (NewUDOld.GetNbrColumn(), NewUDOld.GetNbrColumn(), true);

  SingularValues = EntMatrixOnSite.SingularValueDecomposition(*UMatrix, *VMatrix, false);

  SingularDimension = NewDVMatrix->GetNbrRow();
  if (NewUDOld.GetNbrColumn() < NewDVMatrix->GetNbrRow())
    SingularDimension = NewUDOld.GetNbrColumn();

  double Trace = 0.0;
  for (int i = 0; i < SingularDimension; i++)
    {
      SingularValues[i] *= SingularValues[i];
      Trace += SingularValues[i];
    }  

  cout<<"Entanglement levels: "<<endl;
  for (int i = 0; i < SingularDimension; i++)
    {
      if (fabs(SingularValues[i]) > CutOff)
        cout<<SingularValues[i]/Trace<<endl;
    } 
  cout<<endl; 
  cout<<"Tr= "<<Trace<<endl;
 
  double CutOff2 = 1e-10;

  int LeftRow = NewDVMatrix->GetNbrRow();
  int RightColumn = NewUDOld.GetNbrColumn(); 
  
  cout<<"Left rows " << LeftRow << " right columns " << RightColumn << endl;

  double* IndexMat = new double [LeftRow];
  for (int i = 0; i < LeftRow; i++) 
    IndexMat[i] = 0.0; 

  for (int i = 0; i < LeftRow; i++)
   for (int j = 0; j < RightColumn; j++)
     {
       bool found = false;
       for (int k = 0; (k < MatDim) && (found == false); k++)
         {
            double Tmp1, Tmp2;
            NewDVMatrix->GetMatrixElement(i, k, Tmp1);
            NewUDOld.GetMatrixElement(k, j, Tmp2);
            if ((fabs(Tmp1) > CutOff2) && (fabs(Tmp2) > CutOff2))
              {
                IndexMat[i] = k;
                found = true;
              }
         }
     }

  cout<<"IndexMat = ";
  for (int i = 0; i < LeftRow; i++)
    cout << IndexMat[i] << " ";
  cout<<endl;

  //cout<<"Indices"<<endl;
  //for (int i = 0; i < MatDim; i++)
  // {
  //   int TmpP, TmpN;
  //   MPSMatrix->GetChargeAndPLevelFromMatrixIndex(i, TmpP, TmpN);
  //   cout<<"i= "<<i<<" P= "<<TmpP<<" N= "<<TmpN<<endl;
  // }

  int p = 0;
  int q = 0;
  if (Manager.GetBoolean("k-2") == true) 
    {
       p = 2;
       q = 2 + Manager.GetInteger("r-index");
    }
  else
   {
    if (Manager.GetBoolean("rr-3") == true)
      {
        p = 3;
        q = 5;
      }
     else
      {
        p = 1;
        q = Manager.GetInteger("laughlin-index");
      }  
   } 
  cout << "Filling factor nu = p/q = "<<p<<"/"<<q<<endl;
  int gcd = FindGCD(p, q);
  if (gcd > 1)
    {
      p /= gcd;
      q /= gcd;
      cout << "Filling factor nu* = p/q = "<<p<<"/"<<q<<endl;
    }

  int MinQValue = 0;
  int MaxQValue = 0;
  cout << "warning, code not compatible with multi Q range" << endl;
  MPSMatrix->GetChargeIndexRange(0, MinQValue, MaxQValue);
  cout<<"MinQvalue = "<< MinQValue << " MaxQValue= "<< MaxQValue << endl;

  bool* RowIndicator = new bool[UMatrix->GetNbrRow()];
  for (int i = 0; i < UMatrix->GetNbrRow(); ++i)
    RowIndicator[i] = false;
  int TmpParticlesA, TmpP, TmpN;

  cout<<"Entanglement levels with quantum numbers La,Na,P: "<<endl;

  for (int i = 0; i < SingularDimension; i++)
    {
      if (fabs(SingularValues[i]) > CutOff)
       {
          UMatrix->GetMatrixElement (i, i, Tmp);
          if (fabs(Tmp) > CutOff2)
            {
               RowIndicator[i] = true; 
               //MPSMatrix->GetPNFromMatrixIndex (IndexMat[i], TmpN, TmpP);
               MPSMatrix->GetChargeAndPLevelFromMatrixIndex(IndexMat[i], TmpP, TmpN);

               if ((p == 2) && (q == 4)) //Moore-Read
                 TmpN -= (MaxQValue - 1)/2;
               else
                 TmpN -= MaxQValue/2;

               if (Manager.GetBoolean("use-padding") == false)
                  TmpN -= p;

               TmpN =  (p * EntCut - TmpN)/q; 

               if (Manager.GetBoolean("all-na"))
                 {
                    cout << "Na= " << TmpN << " P= " << TmpP << " " << SingularValues[i]/Trace << endl;
                    File << TmpN << " " << TmpP << " " << SingularValues[i]/Trace << endl;
                 }
               else if (TmpN == Na)
                 {
                  cout << "Na= " << TmpN << " P= " << TmpP << " " << SingularValues[i]/Trace << endl;
                  File << TmpN << " " << TmpP << " " << SingularValues[i]/Trace << endl;
                 }
               //cout << "N= " << TmpN << " P= " << TmpP << " " << SingularValues[i]/Trace << endl;
            }
          else
           {
             bool found = false;
             for (int j = 0; ((j < UMatrix->GetNbrRow()) && (found == false)); ++j)
               {
                 if (RowIndicator[j] == false)
                   {
                     UMatrix->GetMatrixElement(j, i, Tmp);
                     if ( (fabs(Tmp) > CutOff2) && (IndexMat[j] != 0) )
                       {
                          found = true;
                          RowIndicator[j] = true;
                          //MPSMatrix->GetPNFromMatrixIndex (IndexMat[j], TmpN, TmpP);
                          MPSMatrix->GetChargeAndPLevelFromMatrixIndex(IndexMat[j], TmpP, TmpN);
                          if ((p == 2) && (q == 4)) //Moore-Read
                             TmpN -= (MaxQValue - 1)/2;
                          else
                             TmpN -= MaxQValue/2;

                          if (Manager.GetBoolean("use-padding") == false)
                             TmpN -= p;

                          TmpN =  (p * EntCut - TmpN)/q; 

                          if (Manager.GetBoolean("all-na"))
                           {
                             cout << "Na= " << TmpN << " P= " << TmpP << " " << SingularValues[i]/Trace << endl;
                             File << TmpN << " " << TmpP << " " << SingularValues[i]/Trace << endl;
                           }
                          else if (TmpN == Na)
                           { 
                              cout << "Na= " << TmpN << " P= " << TmpP << " " << SingularValues[i]/Trace << endl;
                              File << TmpN << " " << TmpP << " " << SingularValues[i]/Trace << endl;
                           }
                          //cout << "N= " << TmpN << " P= " << TmpP << " " << SingularValues[i]/Trace << endl;
                       }
                   }
               } 
           }  
       }
    } 
  cout<<endl; 

  delete[] RowIndicator;
  delete UMatrix;
  delete VMatrix;
  delete[] IndexMat;
  delete NewDVMatrix;
  delete[] SingularValues;  

  File.close();
 
  return 0;
}
