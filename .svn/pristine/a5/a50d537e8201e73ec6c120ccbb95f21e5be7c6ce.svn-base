#ifndef _REALMPOPERATOROBC_H
#define _REALMPOPERATOROBC_H

#include "AbstractMPOperatorOBC.h"
#include "RealMPSSite.h"
#include "Tensor/Tensor3.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Matrix/RealSymmetricMatrix.h"

class RealMPOperatorOBC : public AbstractMPOperatorOBC
{
 protected:

  double * LeftVector;
  double * RightVector;
  double * ElementsValues;

  
  RealMPOperatorOBC ();
  ~RealMPOperatorOBC ();
  
  virtual void InitializeTensorsElements(){ cout <<"using  not defined RealMPOperatorOBC::InitializeTensorsElements()"<<endl;}

 
 public:
  
  virtual void ComputeL(Tensor3<double> & L);
  virtual void ComputeLCore(Tensor3<double> & L);
  virtual void ComputeLCoreBis(Tensor3<double> & L);

  virtual void ComputeR(Tensor3<double> & R);
  virtual void ComputeRCore(Tensor3<double> & R);
  virtual void ComputeRCoreBis(Tensor3<double> & R);
  virtual  void PrintTensorElements();


  virtual RealVector& LowLevelMultiplyCore(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent);
 

  virtual RealVector& LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);

  virtual RealVector& LowLevelMultiplyTwoSitesCore(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);
  virtual RealVector& LowLevelMultiplyTwoSitesCoreBis(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);

 virtual RealSymmetricMatrix& GetTwoSitesHamiltonian (RealSymmetricMatrix & M);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelMultiplyOneSite(RealVector& vSource, RealVector& vDestination,   int firstComponent, int nbrComponent);
  

  virtual void MPOApplyOnTensorOnTheLeftCore(Tensor3<double> * result, Tensor3<double> * source, int firstComponent, int nbrComponent);
  virtual void MPOApplyOnTensorOnTheRightCore(Tensor3<double> * result, Tensor3 <double> * source, int firstComponent, int nbrComponent);
  virtual void LowLevelMultiplyCoreFirst(Tensor3<double> * result, Tensor3<double> * source , RealVector & vSource, int firstComponent, int nbrComponent);
  virtual void LowLevelMultiplyCoreTwoSitesFirst(Tensor3<double> * result, Tensor3<double> * source , RealVector & vSource, int firstComponent, int nbrComponent);
  virtual void LowLevelMultiplyCoreSecond(Tensor3<double> * result, Tensor3<double> * source , RealVector & vSource, int firstComponent, int nbrComponent);
  virtual void LowLevelMultiplyCoreTwoSitesSecond(Tensor3<double> * rightTensor, Tensor3<double> * source , RealVector & vDestination, int firstComponent, int nbrComponent);

};

#endif
