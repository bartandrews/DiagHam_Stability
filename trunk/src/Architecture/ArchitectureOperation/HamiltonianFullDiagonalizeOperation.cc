////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of hamiltonian full diagonalization operation           //
//                                                                            //
//                        last modification : 06/01/2012                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "Architecture/ArchitectureOperation/HamiltonianFullDiagonalizeOperation.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "Matrix/RealMatrix.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>


using std::cout;
using std::endl;


#ifdef __SCALAPACK__

// binding to the BLACS SL_init function
//
extern "C" void FORTRAN_NAME(sl_init)(const int* context, const int* nprow, const int* npcol);

// binding to the BLACS BLACS_GRIDINFO function
//
extern "C" void FORTRAN_NAME(blacs_gridinfo)(const int* context, const int* nprow, const int* npcol, const int* myrow, const int* mycol);

// binding to the BLACS BLACS_GRIDEXIT function
//
extern "C" void FORTRAN_NAME(blacs_gridexit)(const int* context); 

// binding to the BLACS ZGEBS2D function
//
extern "C" void FORTRAN_NAME(zgebs2d)(const int* context, const char* scope, const char* top, const int* m, const int* n, const doublecomplex* matrix, const int* lda);

// binding to the BLACS ZGEBR2D function
//
extern "C" void FORTRAN_NAME(zgebr2d)(const int* context, const char* scope, const char* top, const int* m, const int* n, const doublecomplex* matrix, const int* lda, const int* myrow, const int* mycol);

// binding to the SCALAPACK DESCINIT function
//
extern "C" void FORTRAN_NAME(descinit)(const int* desc, const int* m, const int* n, const int* mblock, const int* nblock, 
				       const int* rsrc, const int* csrc, const int* context, const int* mxllda, const int* info); 

// binding to the SCALAPACK NUMROC function
//
extern "C" int FORTRAN_NAME(numroc)(const int* n, const int* nblock, const int* iproc, const int* isrcproc, const int* nprocs); 

// binding to the SCALAPACK PDLAMCH function
//
extern "C" double FORTRAN_NAME(pdlamch)(const int* context, const char* cmach);

// binding to the SCALAPACK PZELSET function
//
extern "C" double FORTRAN_NAME(pzelset) (const doublecomplex* localMatrix, const int* rowIndex, const int* columnIndex, const int* desc, const doublecomplex* element);

// binding to the SCALAPACK PDELSET function
//
extern "C" double FORTRAN_NAME(pdelset) (const double* localMatrix, const int* rowIndex, const int* columnIndex, const int* desc, const double* element);


// binding to the SCALAPACK PZHEEVX function
//
extern "C" void FORTRAN_NAME(pzheevx) (const char* jobz, const char* range, const char* uplo, 
				       const int* dimension, const doublecomplex* matrix, 
				       const int* localStartingRowIndex, const int* localStartingColumnIndex, const int* desc, 
				       const double* spectrumLowerBound, const double* spectrumUpperBound, 
				       const int* spectrumLowerIndex, const int* spectrumUpperIndex, 
				       const double* absoluteErrorTolerance, const int* nbrFoundEigenvalues,
				       const int* nbrFoundEigenstates, const double* eigenvalues, 
				       const double* orthogonalizationFactor, const doublecomplex* eigenstates, 
				       const int* localRowEigenstateIndex, const int* localColumnEigenstateIndex, const int* descEigenstateMatrix, 
				       const doublecomplex* work, const int* lwork, 
				       const double* rwork, const int* lrwork, 
				       const int* iwork, const int* liwork, 
				       const int* iFail, const int* iCluster, const double* gap, const int* info);

// binding to the SCALAPACK PZHEEV function
//
extern "C" void FORTRAN_NAME(pzheev) (const char* jobz, const char* uplo, 
				      const int* dimension, const doublecomplex* matrix, 
				      const int* localStartingRowIndex, const int* localStartingColumnIndex, const int* desc, 
				      const double* eigenvalues, const doublecomplex* eigenstates, 
				      const int* localRowEigenstateIndex, const int* localColumnEigenstateIndex, const int* descEigenstateMatrix, 
				      const doublecomplex* work, const int* lwork, 
				      const doublecomplex* rwork, const int* lrwork, 
				      const int* info);


// binding to the SCALAPACK PDSYEVX function
//
extern "C" void FORTRAN_NAME(pdsyevx) (const char* jobz, const char* range, const char* uplo, 
				       const int* dimension, const double* matrix, 
				       const int* localStartingRowIndex, const int* localStartingColumnIndex, const int* desc, 
				       const double* spectrumLowerBound, const double* spectrumUpperBound, 
				       const int* spectrumLowerIndex, const int* spectrumUpperIndex, 
				       const double* absoluteErrorTolerance, const int* nbrFoundEigenvalues,
				       const int* nbrFoundEigenstates, const double* eigenvalues, 
				       const double* orthogonalizationFactor, const double* eigenstates, 
				       const int* localRowEigenstateIndex, const int* localColumnEigenstateIndex, const int* descEigenstateMatrix, 
				       const double* work, const int* lwork, 
				       const int* iwork, const int* liwork, 
				       const int* iFail, const int* iCluster, const double* gap, const int* info);


// binding to the SCALAPACK PDSYEV function
//
extern "C" void FORTRAN_NAME(pdsyev) (const char* jobz, const char* uplo, 
				      const int* dimension, const double* matrix, 
				      const int* localStartingRowIndex, const int* localStartingColumnIndex, const int* desc, 
				      const double* eigenvalues, const double* eigenstates, 
				      const int* localRowEigenstateIndex, const int* localColumnEigenstateIndex, const int* descEigenstateMatrix, 
				      const double* work, const int* lwork, 
				      const int* info);


#endif


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// complexFlag = true if the hamiltonian is complex
// eigenstateFlag = true if the eigenstates have to be computed
// nbrEigenstates = number of eigenstates that have to be computed (<=0 if all eigenstates have to be computed)

HamiltonianFullDiagonalizeOperation::HamiltonianFullDiagonalizeOperation (AbstractHamiltonian* hamiltonian, bool complexFlag, bool eigenstateFlag, int nbrEigenstates)
{
  this->Hamiltonian = hamiltonian;
  if (this->Hamiltonian->IsHermitian() == true)
    {
      this->UseHermitianFlag = true;
    }
  this->OperationType = AbstractArchitectureOperation::HamiltonianFullDiagonalize;
  this->ComplexFlag = complexFlag;
  this->EigenstateFlag = eigenstateFlag;
  this->NbrEigenstates = nbrEigenstates;
  if (this->NbrEigenstates <= 0)
    this->NbrEigenstates = this->Hamiltonian->GetHilbertSpaceDimension();
}

// copy constructor 
//
// operation = reference on operation to copy

HamiltonianFullDiagonalizeOperation::HamiltonianFullDiagonalizeOperation(const HamiltonianFullDiagonalizeOperation& operation)
{
  this->Hamiltonian = operation.Hamiltonian;
  if (this->Hamiltonian->IsHermitian() == true)
    {
      this->UseHermitianFlag = true;
    }
  this->OperationType = AbstractArchitectureOperation::HamiltonianFullDiagonalize;
  this->ComplexFlag = operation.ComplexFlag;
  this->EigenstateFlag = operation.EigenstateFlag;
  this->NbrEigenstates = operation.NbrEigenstates;
}
  
// constructor from a master node information
//
// hamiltonian = pointer to the hamiltonian to use
// architecture = pointer to the distributed architecture to use for communications

HamiltonianFullDiagonalizeOperation::HamiltonianFullDiagonalizeOperation(AbstractHamiltonian* hamiltonian, SimpleMPIArchitecture* architecture)
{
  this->Hamiltonian = hamiltonian;
  if (this->Hamiltonian->IsHermitian() == true)
    {
      this->UseHermitianFlag = true;
    }
  this->OperationType = AbstractArchitectureOperation::HamiltonianFullDiagonalize;
  int TmpFlag = 0;
  architecture->BroadcastToSlaves(TmpFlag);
  if (TmpFlag == 1)
    this->ComplexFlag = true;
  else
    this->ComplexFlag = false;
  TmpFlag = 0;
  architecture->BroadcastToSlaves(TmpFlag);
  if (TmpFlag == 1)
    this->EigenstateFlag = true;
  else
    this->EigenstateFlag = false;
  TmpFlag = 0;
  architecture->BroadcastToSlaves(this->NbrEigenstates);
}
  
// destructor
//

HamiltonianFullDiagonalizeOperation::~HamiltonianFullDiagonalizeOperation()
{
}
  
// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* HamiltonianFullDiagonalizeOperation::Clone()
{
  return new HamiltonianFullDiagonalizeOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool HamiltonianFullDiagonalizeOperation::RawApplyOperation()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  gettimeofday (&(TotalStartingTime2), 0);

  gettimeofday (&(TotalEndingTime2), 0);
  this->ExecutionTime = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  return true;
}

// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool HamiltonianFullDiagonalizeOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __SCALAPACK__
  if (architecture->IsMasterNode())
    {
      if (architecture->RequestOperation(this->OperationType) == false)
	{
	  return false;
	}
      
      int TmpFlag = 0;
      if (this->ComplexFlag == true)
	TmpFlag = 1;
      architecture->BroadcastToSlaves(TmpFlag);
      TmpFlag = 0;
      if (this->EigenstateFlag == true)
	TmpFlag = 1;
      architecture->BroadcastToSlaves(TmpFlag);
      architecture->BroadcastToSlaves(this->NbrEigenstates);
    }
  
  
  long TmpMinimumIndex = 0;
  long TmpMaximumIndex = 0;
  architecture->GetTypicalRange(TmpMinimumIndex, TmpMaximumIndex);
  
  int Context;
  int NbrNodePerColumn = architecture->GetNbrNodes();
  int NbrNodePerRow = (int) sqrt((double) NbrNodePerColumn);
  while ((NbrNodePerRow > 1) && ((NbrNodePerColumn % NbrNodePerRow) != 0))
    --NbrNodePerRow;
  NbrNodePerColumn /= NbrNodePerRow;
  if (architecture->IsMasterNode())
    {
      cout << "nbr of nodes per row = " << NbrNodePerRow << endl;
      cout << "nbr of nodes per column = " << NbrNodePerColumn << endl;
    }
  FORTRAN_NAME(sl_init) (&Context, &NbrNodePerRow, &NbrNodePerColumn);
  
  int LocalNodeRow;
  int LocalNodeColumn;
  FORTRAN_NAME(blacs_gridinfo) (&Context, &NbrNodePerRow, &NbrNodePerColumn, &LocalNodeRow, &LocalNodeColumn);
  
  int* Desc = new int[9];
  int TmpGlobalNbrRow = this->Hamiltonian->GetHilbertSpaceDimension();
  int TmpGlobalNbrColumn = this->Hamiltonian->GetHilbertSpaceDimension();
  int NbrRowPerBlock = 64;
  int NbrColumnPerBlock = 64;
  int Information = 0;
  int TmpZero = 0;
  int LocalLeadingDimensionRow = FORTRAN_NAME(numroc) (&TmpGlobalNbrRow, &NbrRowPerBlock, &LocalNodeRow, &TmpZero, &NbrNodePerRow);
  int LocalLeadingDimensionColumn = FORTRAN_NAME(numroc) (&TmpGlobalNbrColumn, &NbrColumnPerBlock, &LocalNodeColumn, &TmpZero, &NbrNodePerColumn);
  FORTRAN_NAME(descinit) (Desc, &TmpGlobalNbrRow, &TmpGlobalNbrColumn, &NbrRowPerBlock, &NbrColumnPerBlock,
			  &TmpZero, &TmpZero, &Context, &LocalLeadingDimensionRow, &Information);
  
  
  int* DescEingenstateMatrix = new int[9];
  FORTRAN_NAME(descinit) (DescEingenstateMatrix, &TmpGlobalNbrRow, &TmpGlobalNbrColumn, &NbrRowPerBlock, &NbrColumnPerBlock,
			  &TmpZero, &TmpZero, &Context, &LocalLeadingDimensionRow, &Information);
  
  
  if (this->ComplexFlag == true)
    {
      doublecomplex* LocalScalapackMatrix = new doublecomplex[((long) LocalLeadingDimensionRow) * ((long) LocalLeadingDimensionColumn)];
      
      timeval TotalStartingTime;
      if (architecture->IsMasterNode())
	gettimeofday (&TotalStartingTime, 0);
      
      Complex Tmp;
      doublecomplex TmpElement;
      ComplexVector InputVector (this->Hamiltonian->GetHilbertSpaceDimension());
      ComplexVector OutputVector (this->Hamiltonian->GetHilbertSpaceDimension());
      for (int j = 1; j <= TmpGlobalNbrRow; ++j)
	{
	  int TmpNode = architecture->GetNodeIDFromIndex(j - 1);
	  if (TmpNode ==  architecture->GetNodeNbr())
	    {
	      if (this->Hamiltonian->IsHermitian() == true)
		{
		  InputVector.ClearVector();
		  InputVector[j - 1] = 1.0;
		  this->Hamiltonian->HermitianMultiply(InputVector, OutputVector, j - 1, 1);
		}
	      else
		{
		  InputVector[j - 1] = 1.0;
		  this->Hamiltonian->Multiply(InputVector, OutputVector, j - 1, 1);
		}
	    }
	  architecture->BroadcastVector(TmpNode, OutputVector);
	  for (int i = 1; i <= TmpGlobalNbrRow; ++i)
	    {
	      TmpElement.r = OutputVector[i - 1].Re;
	      TmpElement.i = OutputVector[i - 1].Im;
	      FORTRAN_NAME(pzelset) (LocalScalapackMatrix, &i, &j, Desc, &TmpElement);
	    }
	}

      if (architecture->IsMasterNode())
	{
	  timeval TotalEndingTime;
	  gettimeofday (&TotalEndingTime, 0);
	  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
			(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));	
	  if (architecture->VerboseMode())
	    {
	      char TmpString[256];
	      sprintf (TmpString, "HamiltonianFullDiagonalizeOperation fill matrix operation done in %.3f seconds", Dt);
	      architecture->AddToLog(TmpString, true);
	    }
	  cout << "HamiltonianFullDiagonalizeOperation fill matrix operation done in " << Dt << " seconds" << endl;
	  gettimeofday (&TotalStartingTime, 0);
	}
      
      Information = 0; 
      const char* JobZ = "N";
      if (this->EigenstateFlag == true)
	JobZ = "V";
      const char* UpperLower = "U";
      int LocalStartingRowIndex = 1;
      int LocalStartingColumnIndex = 1;
      double* Eigenvalues = new double [TmpGlobalNbrRow];
      doublecomplex* Eigenstates = 0;
      if (this->EigenstateFlag == true)
	{
	  Eigenstates = new doublecomplex [((long) LocalLeadingDimensionRow) * ((long) LocalLeadingDimensionColumn)];
	}
      int LocalRowEigenstateIndex = 1;
      int LocalColumnEigenstateIndex = 1;
      doublecomplex* ScalapackWorkingArea = new doublecomplex[1];
      int ScalapackWorkingAreaSize = -1;
      doublecomplex* ScalapackRWorkingArea = new  doublecomplex[1];
      int ScalapackRWorkingAreaSize = -1; 
      
      for (int i = 0; i < TmpGlobalNbrRow; ++i)
	Eigenvalues[i] = 0.0;
      
      if (architecture->IsMasterNode())
	{
	  cout << "starting diagonalization" << endl;
	}
      
      FORTRAN_NAME(pzheev)(JobZ, UpperLower, 
			   &TmpGlobalNbrRow, LocalScalapackMatrix, 
			   &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
			   Eigenvalues,
			   Eigenstates, 
			   &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
			   ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
			   ScalapackRWorkingArea, &ScalapackRWorkingAreaSize, 
			   &Information);  
      
      ScalapackWorkingAreaSize = (int) ScalapackWorkingArea[0].r;
      ScalapackRWorkingAreaSize = (int) ScalapackRWorkingArea[0].r;
      delete[] ScalapackWorkingArea;
      delete[] ScalapackRWorkingArea;
      ScalapackWorkingArea = new doublecomplex[ScalapackWorkingAreaSize];
      ScalapackRWorkingArea = new doublecomplex [ScalapackRWorkingAreaSize];
      
      FORTRAN_NAME(pzheev)(JobZ, UpperLower, 
			   &TmpGlobalNbrRow, LocalScalapackMatrix, 
			   &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
			   Eigenvalues,
			   Eigenstates, 
			   &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
			   ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
			   ScalapackRWorkingArea, &ScalapackRWorkingAreaSize, 
			   &Information);  

      if (architecture->IsMasterNode())
	{
	  timeval TotalEndingTime;
	  gettimeofday (&TotalEndingTime, 0);
	  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 			
			(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
	  cout << "HamiltonianFullDiagonalizeOperation diagonalization operation done in " << Dt << "  seconds" << endl;
	  if (architecture->VerboseMode())
	    {
	      char TmpString[256];
	      sprintf (TmpString, "HamiltonianFullDiagonalizeOperation diagonalization operation done in %.3f seconds", Dt);
	      architecture->AddToLog(TmpString, true);
	    }
	  gettimeofday (&TotalStartingTime, 0);
	}
      
      
      if (architecture->IsMasterNode())
 	{
 	  this->DiagonalizedMatrix = RealDiagonalMatrix (Eigenvalues, this->Hamiltonian->GetHilbertSpaceDimension());
	  if (this->EigenstateFlag == true)
	    {
	      if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
		gettimeofday (&TotalStartingTime, 0);
	      this->ComplexEigenstates = ComplexMatrix (this->Hamiltonian->GetHilbertSpaceDimension(), this->NbrEigenstates, true);
	      for (int SlaveID = 0; SlaveID < architecture->GetNbrSlaveNodes(); ++SlaveID)
		{
		  int TmpLocalNodeRow = 0;
		  int TmpLocalNodeColumn = 0;
		  int TmpLocalLeadingDimensionRow = 0;
		  int TmpLocalLeadingDimensionColumn = 0;
		  int TmpNbr = 1;
		  architecture->ReceiveFromSlave(SlaveID, &TmpLocalNodeRow, TmpNbr);
		  architecture->ReceiveFromSlave(SlaveID, &TmpLocalNodeColumn, TmpNbr);
		  architecture->ReceiveFromSlave(SlaveID, &TmpLocalLeadingDimensionRow, TmpNbr);
		  architecture->ReceiveFromSlave(SlaveID, &TmpLocalLeadingDimensionColumn, TmpNbr);
		  long TmpSize = ((long) TmpLocalLeadingDimensionRow) * ((long) TmpLocalLeadingDimensionColumn);
		  doublecomplex* TmpEigenstates = new doublecomplex [TmpSize];
		  architecture->ReceiveFromSlave(SlaveID, TmpEigenstates, TmpSize);
		  int TmpNbrBlockPerRow = TmpLocalLeadingDimensionRow / NbrRowPerBlock;
		  if ((TmpNbrBlockPerRow * NbrRowPerBlock) !=  TmpLocalLeadingDimensionRow)
		    ++TmpNbrBlockPerRow;
		  int TmpNbrBlockPerColumn = TmpLocalLeadingDimensionColumn / NbrColumnPerBlock;
		  if ((TmpNbrBlockPerColumn * NbrColumnPerBlock) !=  TmpLocalLeadingDimensionColumn)
		    ++TmpNbrBlockPerColumn;
		  Complex Tmp;
		  cout << "SlaveID = " << SlaveID << " TmpLocalNodeRow = " << TmpLocalNodeRow
		       << " TmpLocalNodeColumn = " << TmpLocalNodeColumn << " TmpLocalLeadingDimensionRow = " << TmpLocalLeadingDimensionRow
		       << " TmpLocalLeadingDimensionColumn = " << TmpLocalLeadingDimensionColumn << endl;
		  cout << TmpNbrBlockPerColumn << " " <<  NbrColumnPerBlock  << " " << TmpLocalLeadingDimensionColumn << endl;
		  for (int j = 0; j < TmpNbrBlockPerColumn; ++j)
		    {
		      int TmpGlobalColumnIndex = (TmpLocalNodeColumn + j * NbrNodePerColumn) * NbrColumnPerBlock;
		      if (TmpGlobalColumnIndex < this->NbrEigenstates)
			{
			  int TmpLocalNbrColumnInBlock = NbrColumnPerBlock;
			  int TmpGlobalLastColumnIndex = TmpGlobalColumnIndex + TmpLocalNbrColumnInBlock;
			  if (TmpGlobalLastColumnIndex > this->NbrEigenstates)
			    {
			      TmpLocalNbrColumnInBlock -= (TmpGlobalLastColumnIndex - this->NbrEigenstates);
			      TmpGlobalLastColumnIndex = this->NbrEigenstates;
			    }

			  for (int i = 0; i < TmpNbrBlockPerRow; ++i)
			    {
			      int TmpGlobalRowIndex = (TmpLocalNodeRow + i * NbrNodePerRow) * NbrRowPerBlock;
			      int TmpLocalNbrRowInBlock = NbrRowPerBlock;
			      int TmpGlobalLastRowIndex = TmpGlobalRowIndex + TmpLocalNbrRowInBlock;
			      if (TmpGlobalLastRowIndex > this->Hamiltonian->GetHilbertSpaceDimension())
				{
				  TmpLocalNbrRowInBlock -= TmpGlobalLastRowIndex - this->Hamiltonian->GetHilbertSpaceDimension();
				  TmpGlobalLastRowIndex = this->Hamiltonian->GetHilbertSpaceDimension();
				}
			      for (int k = 0; k < TmpLocalNbrColumnInBlock; ++k)
				{
				  long TmpLocalStartingIndex = (((j * NbrColumnPerBlock)  + k) * ((long) TmpLocalLeadingDimensionRow)) + (i * NbrRowPerBlock);				  
				  for (int l = 0; l < TmpLocalNbrRowInBlock; ++l)
				    {
				      Tmp.Re = TmpEigenstates[TmpLocalStartingIndex + l].r;
				      Tmp.Im = TmpEigenstates[TmpLocalStartingIndex + l].i;
				      this->ComplexEigenstates.SetMatrixElement(TmpGlobalRowIndex + l, TmpGlobalColumnIndex + k, Tmp);
				    }
				}
			    }
			}
		    }
		  delete[] TmpEigenstates;	      
		}
	      cout << "done with slaves" << endl;

	      int TmpNbrBlockPerRow = LocalLeadingDimensionRow / NbrRowPerBlock;
	      if ((TmpNbrBlockPerRow * NbrRowPerBlock) !=  LocalLeadingDimensionRow)
		++TmpNbrBlockPerRow;
	      int TmpNbrBlockPerColumn = LocalLeadingDimensionColumn / NbrColumnPerBlock;
	      if ((TmpNbrBlockPerColumn * NbrColumnPerBlock) !=  LocalLeadingDimensionColumn)
		++TmpNbrBlockPerColumn;
	      Complex Tmp;
	      for (int j = 0; j < TmpNbrBlockPerColumn; ++j)
		{
		  int TmpGlobalColumnIndex = (LocalNodeColumn + j * NbrNodePerColumn) * NbrColumnPerBlock;
		  if (TmpGlobalColumnIndex < this->NbrEigenstates)
		    {
		      int TmpLocalNbrColumnInBlock = NbrColumnPerBlock;
		      int TmpGlobalLastColumnIndex = TmpGlobalColumnIndex + TmpLocalNbrColumnInBlock;
		      if (TmpGlobalLastColumnIndex > this->NbrEigenstates)
			{
			  TmpLocalNbrColumnInBlock -= (TmpGlobalLastColumnIndex - this->NbrEigenstates);
			  TmpGlobalLastColumnIndex = this->NbrEigenstates;
			}
		      for (int i = 0; i < TmpNbrBlockPerRow; ++i)
			{
			  int TmpGlobalRowIndex = (LocalNodeRow + i * NbrNodePerRow) * NbrRowPerBlock;
			  int TmpLocalNbrRowInBlock = NbrRowPerBlock;
			  int TmpGlobalLastRowIndex = TmpGlobalRowIndex + TmpLocalNbrRowInBlock;
			  if (TmpGlobalLastRowIndex > this->Hamiltonian->GetHilbertSpaceDimension())
			    {
			      TmpLocalNbrRowInBlock -= TmpGlobalLastRowIndex - this->Hamiltonian->GetHilbertSpaceDimension();
			      TmpGlobalLastRowIndex = this->Hamiltonian->GetHilbertSpaceDimension();
			    }
			  for (int k = 0; k < TmpLocalNbrColumnInBlock; ++k)
			    {
			      long TmpLocalStartingIndex = (((j * NbrColumnPerBlock)  + k) * ((long) LocalLeadingDimensionRow)) + (i * NbrRowPerBlock);
			      for (int l = 0; l < TmpLocalNbrRowInBlock; ++l)
				{
				  Tmp.Re = Eigenstates[TmpLocalStartingIndex + l].r;
				  Tmp.Im = Eigenstates[TmpLocalStartingIndex + l].i;
				  this->ComplexEigenstates.SetMatrixElement(TmpGlobalRowIndex + l, TmpGlobalColumnIndex + k, Tmp);
				}
			    }
			}
		    }
		}
	      if (architecture->IsMasterNode())
		{
		  timeval TotalEndingTime;
		  gettimeofday (&TotalEndingTime, 0);
		  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
				(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
		  if (architecture->VerboseMode())
		    {
		      char TmpString[256];
		      sprintf (TmpString, "HamiltonianFullDiagonalizeOperation reassembling eigenstates done in %.3f seconds", Dt);
		      architecture->AddToLog(TmpString, true);
		    }
		  cout << "HamiltonianFullDiagonalizeOperation reassembling eigenstates done in " << Dt << " seconds" << endl;
		}      
	    }
 	}
      else
	{
	  if (this->EigenstateFlag == true)
	    {
	      long TmpSize = ((long) LocalLeadingDimensionRow) * ((long) LocalLeadingDimensionColumn);
	      int TmpNbrElement = 1;
	      architecture->SendToMaster(&LocalNodeRow, TmpNbrElement);
	      architecture->SendToMaster(&LocalNodeColumn, TmpNbrElement);
	      architecture->SendToMaster(&LocalLeadingDimensionRow, TmpNbrElement);
	      architecture->SendToMaster(&LocalLeadingDimensionColumn, TmpNbrElement);
	      architecture->SendToMaster(Eigenstates, TmpSize);
	    }
	  delete[] Eigenvalues;
	}
      if (this->EigenstateFlag == true)
	delete[] Eigenstates;
      delete[] ScalapackWorkingArea;
      delete[] ScalapackRWorkingArea;
      delete[] LocalScalapackMatrix;
    }
  else
    {
      timeval TotalStartingTime;
      if (architecture->IsMasterNode())
	gettimeofday (&TotalStartingTime, 0);
      
      double* LocalScalapackMatrix = new double[((long) LocalLeadingDimensionRow) * ((long) LocalLeadingDimensionColumn)];
      
      double Tmp;
      RealVector InputVector (this->Hamiltonian->GetHilbertSpaceDimension());
      RealVector OutputVector (this->Hamiltonian->GetHilbertSpaceDimension());
      for (int j = 1; j <= TmpGlobalNbrRow; ++j)
	{
	  int TmpNode = architecture->GetNodeIDFromIndex(j - 1);
	  if (TmpNode ==  architecture->GetNodeNbr())
	    {
	      if (this->Hamiltonian->IsHermitian() == true)
		{
		  InputVector.ClearVector();
		  InputVector[j - 1] = 1.0;
		  this->Hamiltonian->HermitianMultiply(InputVector, OutputVector, j - 1, 1);
		}
	      else
		{
		  InputVector[j - 1] = 1.0;
		  this->Hamiltonian->Multiply(InputVector, OutputVector, j - 1, 1);
		}
	    }
	  architecture->BroadcastVector(TmpNode, OutputVector);
  	  for (int i = 1; i <= TmpGlobalNbrRow; ++i)
	    {
	      FORTRAN_NAME(pdelset) (LocalScalapackMatrix, &i, &j, Desc, &(OutputVector[i - 1]));
	    }
	}
      
      if (architecture->IsMasterNode())
	{
	  timeval TotalEndingTime;
	  gettimeofday (&TotalEndingTime, 0);
	  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
			(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
	  if (architecture->VerboseMode())
	    {
	      char TmpString[256];
	      sprintf (TmpString, "HamiltonianFullDiagonalizeOperation fill matrix operation done in %.3f seconds", Dt);
	      architecture->AddToLog(TmpString, true);
	    }
	  cout << "HamiltonianFullDiagonalizeOperation fill matrix operation done in " << Dt << " seconds" << endl;
	  gettimeofday (&TotalStartingTime, 0);
	}

      Information = 0; 
      const char* JobZ = "N";
      if (this->EigenstateFlag == true)
	JobZ = "V";
      
      const char* UpperLower = "U";
      int LocalStartingRowIndex = 1;
      int LocalStartingColumnIndex = 1;
      int NbrFoundEigenvalues = 0;
      int NbrFoundEigenstates = 0;
      double* Eigenvalues = new double [TmpGlobalNbrRow];
      double* Eigenstates = 0;
      if (this->EigenstateFlag == true)
	Eigenstates = new double [((long) LocalLeadingDimensionRow) * ((long) LocalLeadingDimensionColumn)];
      int LocalRowEigenstateIndex = 1;
      int LocalColumnEigenstateIndex = 1;
      double* ScalapackWorkingArea = new double[1];
      int ScalapackWorkingAreaSize = -1;
  
      for (int i = 0; i < TmpGlobalNbrRow; ++i)
	Eigenvalues[i] = 0.0;
      
      if (architecture->IsMasterNode())
	{
	  cout << "starting diagonalization" << endl;
	}

      FORTRAN_NAME(pdsyev)(JobZ, UpperLower, 
			   &TmpGlobalNbrRow, LocalScalapackMatrix, 
			   &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
			   Eigenvalues,
			   Eigenstates, 
			   &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
			   ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
			   &Information); 

      long ScalapackWorkingAreaSizeLong = (long) ScalapackWorkingArea[0];
      ScalapackWorkingAreaSize = (int) ScalapackWorkingArea[0];
      delete[] ScalapackWorkingArea;
      ScalapackWorkingArea = new double[ScalapackWorkingAreaSizeLong];
      cout << "Scalapack working area size = " << ScalapackWorkingAreaSizeLong << "(long), " << ScalapackWorkingAreaSize << "(int)" << endl;  
      FORTRAN_NAME(pdsyev)(JobZ, UpperLower, 
			   &TmpGlobalNbrRow, LocalScalapackMatrix, 
			   &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
			   Eigenvalues,
			   Eigenstates, 
			   &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
			   ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
			   &Information); 

      if (architecture->IsMasterNode())
	{
	  timeval TotalEndingTime;
	  gettimeofday (&TotalEndingTime, 0);
	  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
			(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
	  if (architecture->VerboseMode())
	    {
	      char TmpString[256];
	      sprintf (TmpString, "HamiltonianFullDiagonalizeOperation diagonalization operation done in %.3f seconds", Dt);
	      architecture->AddToLog(TmpString, true);
	    }
	  cout << "HamiltonianFullDiagonalizeOperation diagonalization operation done in " << Dt << " seconds" << endl;
	  gettimeofday (&TotalStartingTime, 0);
	}

      if (architecture->IsMasterNode())
 	{
 	  this->DiagonalizedMatrix = RealDiagonalMatrix (Eigenvalues, this->Hamiltonian->GetHilbertSpaceDimension());
	  if (this->EigenstateFlag == true)
	    {
	      if ((architecture->IsMasterNode()) && (architecture->VerboseMode()))
		gettimeofday (&TotalStartingTime, 0);
	      this->RealEigenstates = RealMatrix (this->Hamiltonian->GetHilbertSpaceDimension(), this->NbrEigenstates, true);
	      for (int SlaveID = 0; SlaveID < architecture->GetNbrSlaveNodes(); ++SlaveID)
		{
		  int TmpLocalNodeRow = 0;
		  int TmpLocalNodeColumn = 0;
		  int TmpLocalLeadingDimensionRow = 0;
		  int TmpLocalLeadingDimensionColumn = 0;
		  int TmpNbr = 1;
		  architecture->ReceiveFromSlave(SlaveID, &TmpLocalNodeRow, TmpNbr);
		  architecture->ReceiveFromSlave(SlaveID, &TmpLocalNodeColumn, TmpNbr);
		  architecture->ReceiveFromSlave(SlaveID, &TmpLocalLeadingDimensionRow, TmpNbr);
		  architecture->ReceiveFromSlave(SlaveID, &TmpLocalLeadingDimensionColumn, TmpNbr);
		  long TmpSize = ((long) TmpLocalLeadingDimensionRow) * ((long) TmpLocalLeadingDimensionColumn);
		  double* TmpEigenstates = new double [TmpSize];
		  architecture->ReceiveFromSlave(SlaveID, TmpEigenstates, TmpSize);
		  int TmpNbrBlockPerRow = TmpLocalLeadingDimensionRow / NbrRowPerBlock;
		  if ((TmpNbrBlockPerRow * NbrRowPerBlock) !=  TmpLocalLeadingDimensionRow)
		    ++TmpNbrBlockPerRow;
		  int TmpNbrBlockPerColumn = TmpLocalLeadingDimensionColumn / NbrColumnPerBlock;
		  if ((TmpNbrBlockPerColumn * NbrColumnPerBlock) !=  TmpLocalLeadingDimensionColumn)
		    ++TmpNbrBlockPerColumn;
		  for (int j = 0; j < TmpNbrBlockPerColumn; ++j)
		    {
		      int TmpGlobalColumnIndex = (TmpLocalNodeColumn + j * NbrNodePerColumn) * NbrColumnPerBlock;
		      if (TmpGlobalColumnIndex < this->NbrEigenstates)
			{
			  int TmpLocalNbrColumnInBlock = NbrColumnPerBlock;
			  int TmpGlobalLastColumnIndex = TmpGlobalColumnIndex + TmpLocalNbrColumnInBlock;
			  if (TmpGlobalLastColumnIndex > this->NbrEigenstates)
			    {
			      TmpLocalNbrColumnInBlock -= (TmpGlobalLastColumnIndex - this->NbrEigenstates);
			      TmpGlobalLastColumnIndex = this->NbrEigenstates;
			    }

			  for (int i = 0; i < TmpNbrBlockPerRow; ++i)
			    {
			      int TmpGlobalRowIndex = (TmpLocalNodeRow + i * NbrNodePerRow) * NbrRowPerBlock;
			      int TmpLocalNbrRowInBlock = NbrRowPerBlock;
			      int TmpGlobalLastRowIndex = TmpGlobalRowIndex + TmpLocalNbrRowInBlock;
			      if (TmpGlobalLastRowIndex > this->Hamiltonian->GetHilbertSpaceDimension())
				{
				  TmpLocalNbrRowInBlock -= TmpGlobalLastRowIndex - this->Hamiltonian->GetHilbertSpaceDimension();
				  TmpGlobalLastRowIndex = this->Hamiltonian->GetHilbertSpaceDimension();
				}
			      for (int k = 0; k < TmpLocalNbrColumnInBlock; ++k)
				{
				  long TmpLocalStartingIndex = (((j * NbrColumnPerBlock)  + k) * ((long) TmpLocalLeadingDimensionRow)) + (i * NbrRowPerBlock);				  
				  for (int l = 0; l < TmpLocalNbrRowInBlock; ++l)
				    {
				      this->RealEigenstates.SetMatrixElement(TmpGlobalRowIndex + l, TmpGlobalColumnIndex + k, TmpEigenstates[TmpLocalStartingIndex + l]);
				    }
				}
			    }
			}
		    }
		  delete[] TmpEigenstates;	      
		}

	      int TmpNbrBlockPerRow = LocalLeadingDimensionRow / NbrRowPerBlock;
	      if ((TmpNbrBlockPerRow * NbrRowPerBlock) !=  LocalLeadingDimensionRow)
		++TmpNbrBlockPerRow;
	      int TmpNbrBlockPerColumn = LocalLeadingDimensionColumn / NbrColumnPerBlock;
	      if ((TmpNbrBlockPerColumn * NbrColumnPerBlock) !=  LocalLeadingDimensionColumn)
		++TmpNbrBlockPerColumn;
	      for (int j = 0; j < TmpNbrBlockPerColumn; ++j)
		{
		  int TmpGlobalColumnIndex = (LocalNodeColumn + j * NbrNodePerColumn) * NbrColumnPerBlock;
		  if (TmpGlobalColumnIndex < this->NbrEigenstates)
		    {
		      int TmpLocalNbrColumnInBlock = NbrColumnPerBlock;
		      int TmpGlobalLastColumnIndex = TmpGlobalColumnIndex + TmpLocalNbrColumnInBlock;
		      if (TmpGlobalLastColumnIndex > this->NbrEigenstates)
			{
			  TmpLocalNbrColumnInBlock -= (TmpGlobalLastColumnIndex - this->NbrEigenstates);
			  TmpGlobalLastColumnIndex = this->NbrEigenstates;
			}
		      for (int i = 0; i < TmpNbrBlockPerRow; ++i)
			{
			  int TmpGlobalRowIndex = (LocalNodeRow + i * NbrNodePerRow) * NbrRowPerBlock;
			  int TmpLocalNbrRowInBlock = NbrRowPerBlock;
			  int TmpGlobalLastRowIndex = TmpGlobalRowIndex + TmpLocalNbrRowInBlock;
			  if (TmpGlobalLastRowIndex > this->Hamiltonian->GetHilbertSpaceDimension())
			    {
			      TmpLocalNbrRowInBlock -= TmpGlobalLastRowIndex - this->Hamiltonian->GetHilbertSpaceDimension();
			      TmpGlobalLastRowIndex = this->Hamiltonian->GetHilbertSpaceDimension();
			    }
			  for (int k = 0; k < TmpLocalNbrColumnInBlock; ++k)
			    {
			      long TmpLocalStartingIndex = (((j * NbrColumnPerBlock)  + k) * ((long) LocalLeadingDimensionRow)) + (i * NbrRowPerBlock);
			      for (int l = 0; l < TmpLocalNbrRowInBlock; ++l)
				{
				  this->RealEigenstates.SetMatrixElement(TmpGlobalRowIndex + l, TmpGlobalColumnIndex + k, Eigenstates[TmpLocalStartingIndex + l]);
				}
			    }
			}
		    }
		}
	      if (architecture->IsMasterNode())
		{
		  timeval TotalEndingTime;
		  gettimeofday (&TotalEndingTime, 0);
		  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
				(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
		  if (architecture->VerboseMode())
		    {
		      char TmpString[256];
		      sprintf (TmpString, "HamiltonianFullDiagonalizeOperation reassembling eigenstates done in %.3f seconds", Dt);
		      architecture->AddToLog(TmpString, true);
		    }
		  cout <<  "HamiltonianFullDiagonalizeOperation reassembling eigenstates done in " << Dt << " seconds" << endl;
		}      
	    }
 	}
      else
	{
	  if (this->EigenstateFlag == true)
	    {
	      long TmpSize = ((long) LocalLeadingDimensionRow) * ((long) LocalLeadingDimensionColumn);
	      int TmpNbrElement = 1;
	      architecture->SendToMaster(&LocalNodeRow, TmpNbrElement);
	      architecture->SendToMaster(&LocalNodeColumn, TmpNbrElement);
	      architecture->SendToMaster(&LocalLeadingDimensionRow, TmpNbrElement);
	      architecture->SendToMaster(&LocalLeadingDimensionColumn, TmpNbrElement);
	      architecture->SendToMaster(Eigenstates, TmpSize);
	    }
	  delete[] Eigenvalues;
	}
      if (this->EigenstateFlag == true)
	delete[] Eigenstates;
      delete[] ScalapackWorkingArea;
      delete[] LocalScalapackMatrix;
    }

  return true;
#else
  return this->RawApplyOperation();
#endif
}

  
