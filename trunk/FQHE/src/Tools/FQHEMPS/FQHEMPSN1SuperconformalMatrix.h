////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the N=1 superconformal states           //
//                                                                            //
//                        last modification : 11/03/2013                      //
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


#ifndef FQHEMPSN1SUPERCONFORMALMATRIX_H
#define FQHEMPSN1SUPERCONFORMALMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Vector/LongRationalVector.h"


class LongRationalMatrix;
class BosonOnDiskShort;


class FQHEMPSN1SuperconformalMatrix : public FQHEMPSClustered2RMatrix
{

  friend class FQHEMPSEvaluateCFTOperation;

 protected:


  // effective P level truncation that has to be applied to the descendant of the primary field
  int EffectivePLevel;

 public:
  
  // default constructor 
  //
  FQHEMPSN1SuperconformalMatrix();

  // constructor from a file describing the state
  //
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // fileName = name of the file that contains the state description
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSN1SuperconformalMatrix(int pLevel, int nbrBMatrices, char* fileName,  bool trimChargeIndices, bool cylinderFlag = false, double kappa = 1.0, 
				AbstractArchitecture* architecture = 0);

  // constructor from stored B matrices
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSN1SuperconformalMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool trimChargeIndices, bool cylinderFlag = false, double kappa = 1.0);

  // destructor
  //
  ~FQHEMPSN1SuperconformalMatrix();
  
  // create the B matrices for the laughlin state
  //
  // cftDirectory = an optional path to the directory where all the CFT matrices are stored
  // architecture = architecture to use for precalculation
  virtual void CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture);

 protected:

  //// load the specific informations from the file header
  //// 
  //// file = reference on the input file stream
  //// return value = true if no error occurred  
  //virtual bool LoadHeader (ifstream& file);

  //// save the specific informations to the file header 
  //// 
  //// file = reference on the output file stream
  //// return value = true if no error occurred  
  //virtual bool SaveHeader (ofstream& file);

  // compute the scalar product matrices of the Virasoro descendant
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // invCentralCharge3 = reference on the value of three divided by the central charge
  // weight = weight of the primary field that is considered
  // return value = scalar product
  virtual LongRational ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
							       LongRational& centralCharge12, LongRational& invCentralCharge3, LongRational& weight);
  
  // compute the scalar product matrices of the Virasoro descendant, using information from previous levels
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // invCentralCharge3 = reference on the value of three divided by the central charge
  // weight = weight of the primary field that is considered
  // precomputedScalarProduct = matrices where scalar product matrix elements computed for previous levels are stored
  // precomputedScalarProductMaxPLevel = maxixum P level that can be accessed through precomputedScalarProduct
  // basis = basis that related the partitions to their index
  // temporaryOccupationNumber = local temporary to store the occupation numbers 
  // return value = scalar product  
  virtual LongRational ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
							       LongRational& centralCharge12, LongRational& invCentralCharge3, LongRational& weight,
							       LongRationalMatrix* precomputedScalarProduct, int precomputedScalarProductMaxPLevel, 
							       BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber);

  // compute the scalar product matrices of the Virasoro descendant
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // invCentralCharge3 = reference on the value of three divided by the central charge
  // weight = weight of the primary field that is considered
  // return value = scalar product  
  double ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
						 double& centralCharge12, double& invCentralCharge3, double& weight);

  // compute the scalar product matrices of the Virasoro descendant, using information from previous levels
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // invCentralCharge3 = reference on the value of three divided by the central charge
  // weight = weight of the primary field that is considered
  // precomputedScalarProduct = matrices where scalar product matrix elements computed for previous levels are stored
  // precomputedScalarProductMaxPLevel = maxixum P level that can be accessed through precomputedScalarProduct
  // basis = basis that related the partitions to their index
  // temporaryOccupationNumber = local temporary to store the occupation numbers 
  // return value = scalar product
  double ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
						 double& centralCharge12, double& invCentralCharge3, double& weight,
						 RealSymmetricMatrix* precomputedScalarProduct, int precomputedScalarProductMaxPLevel, 
						 BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber);

  // compute the matrix elements of any primary field in the Virasoro descendant basis
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // descendantPosition = location of the primary field
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // weight1 = weight of the primary field that is considered for the left state
  // weight2 = weight of the primary field that is considered for the right state
  // weight = weight of the primary field whose matrix elements are computed
  // return value = matrix element
  virtual LongRational ComputeDescendantMatrixElement (long* partition, int partitionLength, int descendantPosition, int position, 
						       LongRational& centralCharge12, LongRational& weight1, LongRational& weight2, 
						       LongRational& weight);
  
  // compute the matrix elements of any primary field in the Virasoro descendant basis
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // descendantPosition = location of the primary field
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // weight1 = weight of the primary field that is considered for the left state
  // weight2 = weight of the primary field that is considered for the right state
  // weight = weight of the primary field whose matrix elements are computed
  // precomputedDescendantMatrixElement = matrices where matrix elements computed for previous levels are stored
  // precomputedDescendantMatrixElementMaxLeftPLevel = maxixum P level that can be accessed through precomputedDescendantMatrixElement for the left entry
  // precomputedDescendantMatrixElementMaxRightPLevel = maxixum P level that can be accessed through precomputedDescendantMatrixElement for the right entry
  // basis = basis that related the partitions to their index
  // temporaryOccupationNumber = local temporary to store the occupation numbers 
  // return value = matrix element
  virtual LongRational ComputeDescendantMatrixElement (long* partition, int partitionLength, int descendantPosition, int position, 
						       LongRational& centralCharge12, LongRational& weight1, LongRational& weight2, 
						       LongRational& weight, LongRationalMatrix** precomputedDescendantMatrixElement,
						       int precomputedDescendantMatrixElementMaxLeftPLevel, 
						       int precomputedDescendantMatrixElementMaxRightPLevel, 
						       BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber);


};

  

#endif
