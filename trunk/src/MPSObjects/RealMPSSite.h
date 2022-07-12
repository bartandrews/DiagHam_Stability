
#ifndef _REALMPSSITE_H
#define _REALMPSSITE_H

#include "AbstractMPOperatorOBC.h"
#include "AbstractMPSSite.h"
#include "Tensor/Tensor3.h"
#include "Matrix/RealMatrix.h"
#include "GeneralTools/GarbageFlag.h"

class RealMPSSite : public AbstractMPSSite
{
 protected:

   RealMatrix * M;  

   Tensor3<double> * L;
   Tensor3<double> * R;
      
 public:

   RealMPSSite();
   RealMPSSite(unsigned int physicalDimension,  RealMPSSite * siteOnLeft,  RealMPSSite * siteOnRight , unsigned int bondDimension, AbstractMPOperatorOBC * mPOperator);

   ~RealMPSSite();

   // assignement (without duplicating datas)
   //
   // M = matrix to copy
   // return value = reference on modified matrix
   virtual RealMPSSite & operator = (const RealMPSSite & site);

   bool CheckLeftNormalization();   
   bool CheckRightNormalization();
   
   void InitializeWithRandomMatrices();
   
   void InitializeLeft(RealMatrix * newA);
   void InitializeRight(RealMatrix * newB);
   void UpdateFromVector(RealVector * psi);
   void GetMatrixInVectorForm(RealVector *& resultInvector );          
   void BringMInLeftCanonicalForm();
   void BringMInRightCanonicalForm();
   void ComputeDensityMatrixRight();
   void ComputeDensityMatrixLeft();
   void SymmetricUpdateOfTwoSites(RealMPSSite * rightSite, RealVector * psi, RealDiagonalMatrix & singularValues );
   RealVector *  StatePrediction(RealMPSSite * rightSite, RealDiagonalMatrix & SingularValues, RealDiagonalMatrix & OldSingularValues);
  inline Tensor3<double> & GetPreviousL ()const
  { return (* ((RealMPSSite*)this->SiteOnLeft)->L);}

   inline Tensor3<double> & GetNextR ()const
     { return (* ((RealMPSSite*)this->SiteOnRight)->R);}
   inline RealMatrix * GetM() const
     { return this->M;}
};

#endif
