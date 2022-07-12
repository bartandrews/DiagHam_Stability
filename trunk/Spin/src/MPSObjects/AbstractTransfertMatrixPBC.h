#ifndef _ABSTRACTTRANSFERTMATRIXPBC_H
#define _ABSTRACTTRANSFERTMATRIXPBC_H


#include "Hamiltonian/AbstractHamiltonian.h"
#include "MPSObjects/AbstractMPSSite.h"
#include "Tensor/Tensor3.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

class AbstractMPSSite;

class AbstractTransfertMatrixPBC : public AbstractHamiltonian
{
 protected:

  unsigned int NbrNonZeroElements;
  unsigned int PhysicalDimension;
  unsigned int MPOBondDimension;
  AbstractSpinChain * HilbertSpace;

  int * TemporaryArray;
  int ChainLength;
  int * PowerD;
  
  AbstractArchitecture * Architecture;
  
  double *** ValuesNonZeroTensorElementTopLeft;
  int ** NbrNonZeroTensorElementTopLeft;
  int *** IndiceRightNonZeroTensorElementTopLeft;
  int *** IndiceBottomNonZeroTensorElementTopLeft;
  
  // global shift to apply to the diagonal matrix elements
  double HamiltonianShift;
  
  
 public:
  
  AbstractTransfertMatrixPBC ();
  AbstractTransfertMatrixPBC(MultiColumnASCIIFile & tensorElementsFile,AbstractArchitecture * architecture);
  ~AbstractTransfertMatrixPBC ();  
  
  virtual void PrintTensorElements();
  
  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  virtual void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);
  
  // set site to be acted on
  //
  // site = pointer to the siteto use 
  //void SetSite (AbstractMPSSite* site);
  
  // set site to be acted on
  //
  // site = pointer to the siteto use 
  //  void SetSiteLeftAndRight (AbstractMPSSite* siteLeft,AbstractMPSSite* siteRight);

  // set site to be acted on
  //
  // site = pointer to the siteto use 
  //  void SetDMRGFlag (bool newFlag){this->DMRGFlag=newFlag;};
  
  // set site to be acted on
  //
  // site = pointer to the siteto use 
  //  void SetIDMRGFlag (bool newFlag){this->IDMRGFlag=newFlag;};
  

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  virtual AbstractHilbertSpace* GetHilbertSpace ();
  
  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  virtual int GetHilbertSpaceDimension ();
  int GetTwoSitesHilbertSpaceDimension ();  
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);
  
  inline int GetMPODimension() const {return  this->MPOBondDimension;};
  inline int GetBondDimension() const {return  this->MPOBondDimension;};
  inline int GetPhysicalDimension() const {return  this->PhysicalDimension;};
  
  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination);
  //  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  //  virtual RealVector& LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination);
  //  virtual ComplexVector& LowLevelMultiplyTwoSites(ComplexVector& vSource, ComplexVector& vDestination);
  

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);
  //  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);


  //  virtual RealVector& LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);
  //  virtual ComplexVector& LowLevelMultiplyTwoSites(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);
  
  //  virtual RealSymmetricMatrix& GetTwoSitesHamiltonian (RealSymmetricMatrix & M);
  //  virtual HermitianMatrix & GetTwoSitesHamiltonian (HermitianMatrix & M);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  //virtual RealVector& LowLevelMultiplyOneSite(RealVector& vSource, RealVector& vDestination,  int firstComponent, int nbrComponent);
  //virtual ComplexVector& LowLevelMultiplyOneSite(ComplexVector& vSource, ComplexVector& vDestination,  int firstComponent, int nbrComponent);

  
  virtual void LowLevelMultiplyCoreFirst(Tensor3<Complex> * result, Tensor3<Complex> * source , ComplexVector & vSource, int firstComponent, int nbrComponent)
  {cout <<"using undefined function"<<endl; };
  virtual void LowLevelMultiplyCoreFirst(Tensor3<double> * result, Tensor3<double> * source , RealVector & vSource, int firstComponent, int nbrComponent)
  {cout <<"using undefined function"<<endl; };
  virtual void LowLevelMultiplyCoreTwoSitesFirst(Tensor3<double> * result, Tensor3<double> * source , RealVector & vSource, int firstComponent, int nbrComponent)
  {cout <<"using undefined function"<<endl; };
  virtual void LowLevelMultiplyCoreTwoSitesFirst(Tensor3<Complex> * result, Tensor3<Complex> * source , ComplexVector & vSource, int firstComponent, int nbrComponent)
  {cout <<"using undefined function"<<endl; };
  
  virtual void LowLevelMultiplyCoreSecond(Tensor3<Complex> * result, Tensor3<Complex> * source , ComplexVector & vSource, int firstComponent, int nbrComponent)
  {cout <<"using undefined function"<<endl; };
  virtual void LowLevelMultiplyCoreSecond(Tensor3<double> * result, Tensor3<double> * source , RealVector & vSource, int firstComponent, int nbrComponent)
  {cout <<"using undefined function"<<endl; };
  virtual void LowLevelMultiplyCoreTwoSitesSecond(Tensor3<Complex> * result, Tensor3<Complex> * source , ComplexVector & vSource, int firstComponent, int nbrComponent)
  {cout <<"using undefined function"<<endl; };
  virtual void LowLevelMultiplyCoreTwoSitesSecond(Tensor3<double> * result, Tensor3<double> * source , RealVector & vSource, int firstComponent, int nbrComponent)
  {cout <<"using undefined function"<<endl; };
  
  
  virtual int GenerateResultingStateAndCoefficient(int indiceTop, int chainSize, int lastIndice, double * coefArray, unsigned long * stateArray, unsigned long pos);
  
  virtual int EvaluateNbrResultingState(int indiceTop, int chainSize, int lastIndice);
  
 protected:
  
  inline unsigned int GetIndiceDownFromTensorIndex(unsigned int tensorIndex);

  inline unsigned int GetIndiceUpFromTensorIndex(unsigned int tensorIndex);
  
  inline unsigned int GetIndiceLeftFromTensorIndex(unsigned int tensorIndex);
  
  inline unsigned int GetIndiceRightFromTensorIndex(unsigned int tensorIndex);
  
  inline unsigned int GetTensorIndexFromAllIndices(unsigned int indexDown, unsigned int indexUp, unsigned int indexLeft, unsigned int indexRight);

  inline void GetAllIndicesFromTensorIndex(unsigned int tensorIndex, unsigned int & indexDown, unsigned int & indexUp, unsigned int & indexLeft, unsigned int & indexRight);
  
  virtual void InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile);
  
};

inline unsigned int AbstractTransfertMatrixPBC::GetIndiceDownFromTensorIndex(unsigned int tensorIndex)
{
  return  tensorIndex%this->PhysicalDimension;    
}


inline unsigned int AbstractTransfertMatrixPBC::GetIndiceUpFromTensorIndex(unsigned int tensorIndex)
{
  return  (tensorIndex/this->PhysicalDimension)%this->PhysicalDimension;    
}


inline unsigned int AbstractTransfertMatrixPBC::GetIndiceRightFromTensorIndex(unsigned int tensorIndex)
{
  return  (tensorIndex/(this->PhysicalDimension*this->PhysicalDimension)% this-> MPOBondDimension);    
}


inline unsigned int AbstractTransfertMatrixPBC::GetIndiceLeftFromTensorIndex(unsigned int tensorIndex)
{
  return  (tensorIndex/(this->PhysicalDimension*this->PhysicalDimension *this->MPOBondDimension));    
}


inline unsigned int AbstractTransfertMatrixPBC::GetTensorIndexFromAllIndices(unsigned int indexDown, unsigned int indexUp, unsigned int indexLeft, unsigned int indexRight)
{
  return  ((indexLeft * this-> MPOBondDimension  + indexRight)*this->PhysicalDimension +indexUp)*this->PhysicalDimension + indexDown;
}


inline void AbstractTransfertMatrixPBC::GetAllIndicesFromTensorIndex(unsigned int tensorIndex, unsigned int & indexDown, unsigned int & indexUp, unsigned int & indexLeft, unsigned int & indexRight)
{
   indexDown = tensorIndex%this->PhysicalDimension;     
   tensorIndex/=this->PhysicalDimension;
   indexUp =  tensorIndex%this->PhysicalDimension;
   tensorIndex/=this->PhysicalDimension;
   indexRight = tensorIndex% this-> MPOBondDimension;
   indexLeft = tensorIndex/this->MPOBondDimension;
}


#endif

