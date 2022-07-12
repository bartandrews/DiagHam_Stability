////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                        to particles on a sphere with                       //
//                                                                            //
//                        last modification : 23/03/2003                      //
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


#ifndef ABSTRACTQHEONSPHEREHAMILTONIAN_H
#define ABSTRACTQHEONSPHEREHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;
class ParticleOnSphereL2Hamiltonian;


class AbstractQHEOnSphereHamiltonian : public AbstractQHEHamiltonian
{

  friend class QHEParticlePrecalculationOperation;
  friend class LinearlySuperposedQHEOnSphereHamiltonian;

 protected:
  
  // Hilbert space associated to the system
  ParticleOnSphere* Particles;

  // number of particles
  int NbrParticles;

  // maximum Lz value reached by a particle in the state
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;

  // array containing all interaction factors 
  double* InteractionFactors;
  // number of interaction factors
  int NbrInteractionFactors;
  // arrays for indices attached to each interaction factor
  int* M1Value;
  int* M2Value;
  int* M3Value;
  // alternative method used to store indices attached to each interaction factor
  int NbrM12Indices;
  int* NbrM3Values;
  int** M3Values;

  // flag to indicate if there is any one body terms in the Hamiltonian
  bool OneBodyTermFlag;
  // array containing all factors fo the one body terms
  double* OneBodyInteractionFactors;
  // number of one body terms
  int NbrOneBodyInteractionFactors;
  // arrays for indices attached to each one body term
  int* OneBodyMValues;
  int* OneBodyNValues;

  // pointer to an optional L^2 operator in the Hamiltonian 
  ParticleOnSphereL2Hamiltonian* L2Operator;

  // shift to apply to go from precalculation index to the corresponding index in the HilbertSpace
  int PrecalculationShift;

  // amount of memory (in bytes) that can be used to store precalculated matrix elements
  long Memory;

  // flag for implementation of hermitian symmetry
  bool HermitianSymmetryFlag;
  
  // flag for fast multiplication algorithm
  bool FastMultiplicationFlag;
  // step between each precalculated index (main part: start at 0, FastMultiplicationStep, 2*FastMultiplicationStep, ...)
  int FastMultiplicationStep;
  // number of non-null term in the hamiltonian for each state
  int* NbrInteractionPerComponent;
  // index of the state obtained for each term of the hamiltonian when applying on a given state
  int** InteractionPerComponentIndex;
  // multiplicative coefficient obtained for each term of the hamiltonian when applying on a given state and with a given destination state
  double** InteractionPerComponentCoefficient;

  // number of tasks for load balancing
  int NbrBalancedTasks;
  // load balancing array for parallelisation, indicating starting indices
  long* LoadBalancingArray;
  // cumulative count of non-zero matrix elements
  
  
  // flag to indicate if a hamiltonian is temporary stored on disk
  bool DiskStorageFlag;
  // name of the file that contains hamiltonian matrix elements
  char* DiskStorageFileName;
  // index of the first row that appears in the on-disk hamiltonian
  int DiskStorageStart;
  // maximum number of non-null terms in the hamiltonian for each state
  int MaxNbrInteractionPerComponent;
  // size of the in-memory temporary buffer
  long BufferSize;

  // shift to apply to the Hamiltonian diagonal elements
  double HamiltonianShift;



 public:

  // default constructor
  //
  AbstractQHEOnSphereHamiltonian();

  // destructor
  //
  virtual ~AbstractQHEOnSphereHamiltonian() = 0;

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  virtual AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  virtual int GetHilbertSpaceDimension ();
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  virtual void ShiftHamiltonian (double shift);
  
  // ask if Hamiltonian implements methods using hermitian symmetry 
  //
  virtual bool IsHermitian();

  // ask if Hamiltonian implements methods applying the conjugate of the Hamiltonian
  //
  virtual bool IsConjugate();

  // symmetrize interaction factors to enable hermitian matrix multiplication
  // return = true upon success
  virtual bool HermitianSymmetrizeInteractionFactors();

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  virtual Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
                                          int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
                                                  int firstComponent, int nbrComponent);
  
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);


  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						   int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* ConjugateLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							   int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						   int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							   int firstComponent, int nbrComponent);

  
  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  virtual List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  virtual List<Matrix*> RightInteractionOperators();  

  // save precalculations in a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be stored
  // return value = true if no error occurs
  virtual bool SavePrecalculation (char* fileName);

 protected:
 
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
                                                             int firstComponent, int nbrComponent);
  
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
                                                             int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
                                                     int firstComponent, int nbrComponent);
  
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
                                                     int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
                                                                     int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
                                                             int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& ConjugateLowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
                                                             int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& ConjugateLowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
							      int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations,
									      int nbrVectors, int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* ConjugateLowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
                                                             int firstComponent, int nbrComponent);


  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
                                                             int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& HermitianLowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
                                                     int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using partial fast multiply option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations,
									      int nbrVectors, int firstComponent, int nbrComponent);

  // multiply a et of vectors by the current hamiltonian for a given range of indices 
  // and add result to another et of vectors, low level function (no architecture optimization)
  // using disk storage option
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* HermitianLowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
                                                             int firstComponent, int nbrComponent);

  
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
  // nbrThreads = number of threads requested
  // segmentIndices = array returning the reference to an array of the first index of each of the segments
  //
  virtual bool GetLoadBalancing(int nbrTasks, long* &segmentIndices);

  // test the amount of memory needed for fast multiplication algorithm
  //
  // allowedMemory = amount of memory that cam be allocated for fast multiplication
  // return value = amount of memory needed
  virtual long FastMultiplicationMemory(long allowedMemory);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = number of components that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int nbrComponent);

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication();

  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = number of components that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int nbrComponent);
  
  // enable fast multiplication algorithm using on disk cache 
  //
  // fileName = prefix of the name of the file where temporary matrix elements will be stored
  virtual void EnableFastMultiplicationWithDiskStorage(char* fileName);

  // load precalculations from a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be read
  // return value = true if no error occurs
  virtual bool LoadPrecalculation (char* fileName);

  // core part of the FastMultiplication method involving 2-body term
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
						    int* indexArray, double* coefficientArray, long& position);

};

// core part of the FastMultiplication method involving 2-body term
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void AbstractQHEOnSphereHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphere* particles, int index, 
											 int* indexArray, double* coefficientArray, long& position)
{
  if (this->NbrM12Indices == 0)
    {
      //      indexArray = this->InteractionPerComponentIndex[position];
      //      coefficientArray = this->InteractionPerComponentCoefficient[position];
      int Pos = 0;
      int m1;
      int m2;
      int m3;
      int m4;
      int Index = 0;
      double Coefficient = 0.0;
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  m4 = m1 + m2 - m3;
	  Index = particles->AdAdAA(index, m1, m2, m3, m4, Coefficient);
	  if (Index < particles->GetHilbertSpaceDimension())
	    {
	      indexArray[Pos] = Index;
	      coefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
	      ++Pos;
	    }
	}
      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      m1 = this->OneBodyMValues[j];
	      m2 = this->OneBodyNValues[j];
	      Index = particles->AdA(index, m1, m2, Coefficient);
	      if (Index < particles->GetHilbertSpaceDimension())
		{
		  indexArray[Pos] = Index;
		  coefficientArray[Pos] = Coefficient * this->OneBodyInteractionFactors[j];
		  ++Pos;
		}
	    }
	}
      ++position;
    }
  else
    {
      double Coefficient2;
      int SumIndices;
      int TmpNbrM3Values;
      int* TmpM3Values;
      int ReducedNbrInteractionFactors;
      //      indexArray = this->InteractionPerComponentIndex[position];
      //      coefficientArray = this->InteractionPerComponentCoefficient[position];
      int Pos = 0;
      int Index = 0;
      double Coefficient = 0.0;
      ReducedNbrInteractionFactors = 0;
      int AbsoluteIndex = index + this->PrecalculationShift;
      for (int m1 = 0; m1 < this->NbrM12Indices; ++m1)
	{
	  Coefficient = particles->AA(AbsoluteIndex, this->M1Value[m1], this->M2Value[m1]);	  
	  if (Coefficient != 0.0)
	    {
	      SumIndices = this->M1Value[m1] + this->M2Value[m1];
	      TmpM3Values = this->M3Values[m1];
	      TmpNbrM3Values = this->NbrM3Values[m1];
	      for (int m3 = 0; m3 < TmpNbrM3Values; ++m3)
		{
		  Index = particles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
		  if (Index < particles->GetHilbertSpaceDimension())
		    {
		      indexArray[Pos] = Index;
		      coefficientArray[Pos] = Coefficient * Coefficient2 * this->InteractionFactors[ReducedNbrInteractionFactors];
// 		      cout << (this->M1Value[m1]*2 - this->LzMax) << " " << (this->M2Value[m1]*2 - this->LzMax) << " " << (m3*2 - this->LzMax) << " " << (this->M1Value[m1] + this->M2Value[m1] - m3) << " " << coefficientArray[Pos] << endl;
	      //cout << index << " -> " << (this->M1Value[m1]) << " " << (this->M2Value[m1]) << " " << (TmpM3Values[m3]) << " " << (SumIndices - TmpM3Values[m3]) << " : " << (Coefficient) << " " << (Coefficient2) << endl;
		      ++Pos;
		    }		      
		  ++ReducedNbrInteractionFactors;
		}    
	    }
	  else
	    ReducedNbrInteractionFactors += this->NbrM3Values[m1];
	}
      if (this->OneBodyTermFlag == true)
	{
	  int m1;
	  int m2;
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      m1 = this->OneBodyMValues[j];
	      m2 = this->OneBodyNValues[j];
	      Index = particles->AdA(AbsoluteIndex, m1, m2, Coefficient);
	      if (Index < particles->GetHilbertSpaceDimension())
		{
		  indexArray[Pos] = Index;
		  coefficientArray[Pos] = Coefficient * this->OneBodyInteractionFactors[j];
		  ++Pos;
		}
	    }
	}
      ++position;
    }      
}


#endif
