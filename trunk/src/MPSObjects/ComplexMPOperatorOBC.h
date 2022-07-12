#ifndef _COMPLEXMPOPERATOROBC_H
#define _COMPLEXMPOPERATOROBC_H

#include "AbstractMPOperatorOBC.h"
#include "ComplexMPSSite.h"
#include "Tensor/Tensor3.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Matrix/HermitianMatrix.h"

class ComplexMPOperatorOBC : public AbstractMPOperatorOBC
{
 protected:

  Complex * LeftVector;
  Complex * RightVector;
  Complex * ElementsValues;

  
  ComplexMPOperatorOBC ();
  ~ComplexMPOperatorOBC ();
  
  virtual void InitializeTensorsElements(){ cout <<"using  not defined ComplexMPOperatorOBC::InitializeTensorsElements()"<<endl;}

 
 public:

  virtual void PrintTensorElements();  

  virtual void ComputeL(Tensor3<Complex> & L);
  virtual void ComputeLCore(Tensor3<Complex> & L);

  virtual void ComputeR(Tensor3<Complex> & R);
  virtual void ComputeRCore(Tensor3<Complex> & R);
  


  virtual ComplexVector& LowLevelMultiplyCore(ComplexVector & vSource, ComplexVector & vDestination, 
				       int firstComponent, int nbrComponent);

  virtual ComplexVector & LowLevelMultiplyTwoSites(ComplexVector & vSource, ComplexVector & vDestination, int firstComponent, int nbrComponent);

  virtual ComplexVector & LowLevelMultiplyTwoSitesCore(ComplexVector & vSource, ComplexVector & vDestination, int firstComponent, int nbrComponent);

  virtual HermitianMatrix& GetTwoSitesHamiltonian (HermitianMatrix & M);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector & LowLevelMultiplyOneSite(ComplexVector & vSource, ComplexVector & vDestination,   int firstComponent, int nbrComponent);
  virtual void LowLevelMultiplyCoreFirst(Tensor3<Complex> * result, Tensor3<Complex> * source , ComplexVector & vSource, int firstComponent, int nbrComponent);
  virtual void LowLevelMultiplyCoreSecond(Tensor3<Complex> * leftTensor, Tensor3<Complex> * source , ComplexVector & vDestination, int firstComponent, int nbrComponent); 

  virtual void LowLevelMultiplyCoreTwoSitesFirst(Tensor3<Complex> * result, Tensor3<Complex> * source , ComplexVector & vSource, int firstComponent, int nbrComponent);
  virtual void LowLevelMultiplyCoreTwoSitesSecond(Tensor3<Complex> * rightTensor, Tensor3<Complex> * source , ComplexVector & vDestination, int firstComponent, int nbrComponent);

  virtual void MPOApplyOnTensorOnTheLeftCore(Tensor3<Complex> * result, Tensor3<Complex> * source, int firstComponent, int nbrComponent);
  virtual void MPOApplyOnTensorOnTheRightCore(Tensor3<Complex> * result, Tensor3 <Complex> * source, int firstComponent, int nbrComponent);

};

#endif
