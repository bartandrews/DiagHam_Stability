#ifndef _ABSTRACTMPOPERATOROBC_H
#define _ABSTRACTMPOPERATOROBC_H

#include "Hamiltonian/AbstractHamiltonian.h"
#include "AbstractMPSSite.h"
#include "Tensor/Tensor3.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Architecture/AbstractArchitecture.h"

class AbstractMPSSite;

class AbstractMPOperatorOBC : public AbstractHamiltonian
{

 protected:
  unsigned int NbrNonZeroElements;
  unsigned int * IndexValues;
  unsigned int PhysicalDimension;
  unsigned int MPOBondDimension;
  AbstractHilbertSpace* HilbertSpace;
  bool IDMRGFlag;
  AbstractMPSSite * Site;
  AbstractMPSSite * SiteLeft;
  AbstractMPSSite * SiteRight;
  AbstractArchitecture * Architecture;
  AbstractMPOperatorOBC ();
  ~AbstractMPOperatorOBC ();
  
  virtual void InitializeTensorsElements(){ cout <<"using  not defined AbstractMPOperatorOBC::InitializeTensorsElements()"<<endl;}

  // global shift to apply to the diagonal matrix elements
  double HamiltonianShift;

 
 public:
  
  virtual void ComputeL(Tensor3<double> & L);
  virtual void ComputeLBis(Tensor3<double> & L){};
  virtual void ComputeLBis(Tensor3<Complex> & L){};
  virtual void ComputeL(Tensor3<Complex> & L);
  virtual void ComputeR(Tensor3<double> & R);
  virtual void ComputeR(Tensor3<Complex> & R);
  virtual void PrintTensorElements() = 0;

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // set site to be acted on
  //
  // site = pointer to the siteto use 
  void SetSite (AbstractMPSSite* site);

  // set site to be acted on
  //
  // site = pointer to the siteto use 
  void SetSiteLeftAndRight (AbstractMPSSite* siteLeft,AbstractMPSSite* siteRight);

  // set site to be acted on
  //
  // site = pointer to the siteto use 
  void SetDMRGFlag (bool newFlag){this->IDMRGFlag=newFlag;};
  

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();
  
  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();
  int GetTwoSitesHilbertSpaceDimension ();  

  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);
  
  inline int GetMPODimension() const {return  MPOBondDimension;};
  inline int GetPhysicalDimension() const {return  PhysicalDimension;};

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination);
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  virtual RealVector& LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination);
  virtual ComplexVector& LowLevelMultiplyTwoSites(ComplexVector& vSource, ComplexVector& vDestination);


  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);
  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);


  virtual RealVector& LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);
  virtual ComplexVector& LowLevelMultiplyTwoSites(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);

  virtual RealSymmetricMatrix& GetTwoSitesHamiltonian (RealSymmetricMatrix & M);
  virtual HermitianMatrix & GetTwoSitesHamiltonian (HermitianMatrix & M);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelMultiplyOneSite(RealVector& vSource, RealVector& vDestination,  int firstComponent, int nbrComponent);
  virtual ComplexVector& LowLevelMultiplyOneSite(ComplexVector& vSource, ComplexVector& vDestination,  int firstComponent, int nbrComponent);

  virtual void MPOApplyOnTensorOnTheLeftCore(Tensor3<double> * result, Tensor3<double> * source, int firstComponent, int nbrComponent){cout <<"using undefined function"<<endl; };
  virtual void MPOApplyOnTensorOnTheRightCore(Tensor3<double> * result, Tensor3<double> * source, int firstComponent, int nbrComponent){cout <<"using undefined function"<<endl; };
  virtual void MPOApplyOnTensorOnTheLeftCore(Tensor3<Complex> * result, Tensor3<Complex> * source, int firstComponent, int nbrComponent){cout <<"using undefined function"<<endl; };
  virtual void MPOApplyOnTensorOnTheRightCore(Tensor3<Complex> * result, Tensor3<Complex> * source, int firstComponent, int nbrComponent){cout <<"using undefined function"<<endl; };


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

 protected:
  
  inline unsigned int GetIndiceDownFromTensorIndex(unsigned int tensorIndex);

  inline unsigned int GetIndiceUpFromTensorIndex(unsigned int tensorIndex);
  
  inline unsigned int GetIndiceLeftFromTensorIndex(unsigned int tensorIndex);
  
  inline unsigned int GetIndiceRightFromTensorIndex(unsigned int tensorIndex);
  
  inline unsigned int GetTensorIndexFromAllIndices(unsigned int indexDown, unsigned int indexUp, unsigned int indexLeft, unsigned int indexRight);

  inline void GetAllIndicesFromTensorIndex(unsigned int tensorIndex, unsigned int & indexDown, unsigned int & indexUp, unsigned int & indexLeft, unsigned int & indexRight);
 
};

inline unsigned int AbstractMPOperatorOBC::GetIndiceDownFromTensorIndex(unsigned int tensorIndex)
{
  return  tensorIndex%this->PhysicalDimension;    
}


inline unsigned int AbstractMPOperatorOBC::GetIndiceUpFromTensorIndex(unsigned int tensorIndex)
{
  return  (tensorIndex/this->PhysicalDimension)%this->PhysicalDimension;    
}


inline unsigned int AbstractMPOperatorOBC::GetIndiceRightFromTensorIndex(unsigned int tensorIndex)
{
  return  (tensorIndex/(this->PhysicalDimension*this->PhysicalDimension)% this-> MPOBondDimension);    
}




inline unsigned int AbstractMPOperatorOBC::GetIndiceLeftFromTensorIndex(unsigned int tensorIndex)
{
  return  (tensorIndex/(this->PhysicalDimension*this->PhysicalDimension *this->MPOBondDimension));    
}


inline unsigned int AbstractMPOperatorOBC::GetTensorIndexFromAllIndices(unsigned int indexDown, unsigned int indexUp, unsigned int indexLeft, unsigned int indexRight)
{
  return  ((indexLeft * this-> MPOBondDimension  + indexRight)*this->PhysicalDimension +indexUp)*this->PhysicalDimension + indexDown;
}


inline void AbstractMPOperatorOBC::GetAllIndicesFromTensorIndex(unsigned int tensorIndex, unsigned int & indexDown, unsigned int & indexUp, unsigned int & indexLeft, unsigned int & indexRight)
{
   indexDown = tensorIndex%this->PhysicalDimension;     
   tensorIndex/=this->PhysicalDimension;
   indexUp =  tensorIndex%this->PhysicalDimension;
   tensorIndex/=this->PhysicalDimension;
   indexRight = tensorIndex% this-> MPOBondDimension;
   indexLeft = tensorIndex/this->MPOBondDimension;
}


#endif
