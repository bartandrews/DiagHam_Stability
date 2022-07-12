
#ifndef _COMPLEXMPSSITE_H
#define _COMPLEXMPSSITE_H

#include "AbstractMPOperatorOBC.h"
#include "AbstractMPSSite.h"
#include "Tensor/Tensor3.h"
#include "Matrix/ComplexMatrix.h"
#include "GeneralTools/GarbageFlag.h"

class ComplexMPSSite : public AbstractMPSSite
{
 protected:

   ComplexMatrix * M;  

   Tensor3<Complex> * L;
   Tensor3<Complex> * R;
      
 public:

   ComplexMPSSite();
   ComplexMPSSite(unsigned int physicalDimension,  ComplexMPSSite * siteOnLeft, ComplexMPSSite * siteOnRight , unsigned int bondDimension, AbstractMPOperatorOBC * mPOperator);

   ~ComplexMPSSite();

   // assignement (without duplicating datas)
   //
   // M = matrix to copy
   // return value = reference on modified matrix
   virtual ComplexMPSSite & operator = (const ComplexMPSSite & site);

   bool CheckLeftNormalization();   
   bool CheckRightNormalization();
   
   void InitializeWithRandomMatrices();
   
   void InitializeLeft(ComplexMatrix * newA);
   void InitializeRight(ComplexMatrix * newB);
   void UpdateFromVector(ComplexVector * psi);
   void GetMatrixInVectorForm(ComplexVector *& resultInvector );          
   void BringMInLeftCanonicalForm();
   void BringMInRightCanonicalForm();
   void ComputeDensityMatrixRight();
   void ComputeDensityMatrixLeft();
   virtual ComplexVector *  StatePrediction(ComplexMPSSite * rightSite, RealDiagonalMatrix & SingularValues, RealDiagonalMatrix & OldSingularValues);
   void SymmetricUpdateOfTwoSites(ComplexMPSSite * rightSite, ComplexVector * psi, RealDiagonalMatrix & singularValues );

   inline Tensor3<Complex> & GetPreviousL ()const
     { return (* ((ComplexMPSSite*)this->SiteOnLeft)->L);}
   
   inline Tensor3<Complex> & GetNextR ()const
     { return (* ((ComplexMPSSite*)this->SiteOnRight)->R);}
   inline ComplexMatrix * GetM() const
     { return this->M;}
   
};

#endif
