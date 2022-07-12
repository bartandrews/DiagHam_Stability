#ifndef _ABSTRACTMPSSITE_H
#define _ABSTRACTMPSSITE_H

#include "AbstractMPOperatorOBC.h"
#include "Tensor/Tensor3.h"
#include "GeneralTools/GarbageFlag.h"

class AbstractMPOperatorOBC;

class AbstractMPSSite 
{
 protected:

   // garbage flag to avoid data duplication
   GarbageFlag Flag;
   AbstractMPSSite * SiteOnLeft;
   AbstractMPSSite * SiteOnRight;

   unsigned int PhysicalDimension;
   unsigned int SquarePhysicalDimension;
   unsigned int BondDimensionLeft;
   unsigned int BondDimensionRight;
   unsigned int MaxBondDimension;

   
   AbstractMPOperatorOBC * OperatorToBeMinimized;
   
 public:

   AbstractMPSSite();
   AbstractMPSSite(unsigned int physicalDimension, unsigned int bondDimension, AbstractMPOperatorOBC * mPOperator);

   ~AbstractMPSSite();

   // assignement (without duplicating datas)
   //
   // M = matrix to copy
   // return value = reference on modified matrix
   virtual AbstractMPSSite & operator = (const AbstractMPSSite & site);

   
   virtual bool CheckLeftNormalization();
   virtual bool CheckRightNormalization();
   
   virtual void InitializeWithRandomMatrices() = 0;
   
  //   void InitializeLeft(RealMatrix * newA);
   //   void InitializeRight(RealMatrix * newB);
   //   void UpdateFromVector(RealVector * psi);
//   void GetMatrixInVectorForm(RealVector *& resultInvector );          

   virtual void BringMInLeftCanonicalForm();
   virtual void BringMInRightCanonicalForm();
   virtual void BringMInLeftCanonicalFormCareful();
   virtual void BringMInRightCanonicalFormCareful();
   void SetBondDimension(int bondDimensionLeft, int bondDimensionRight);
   inline void SetRightDimension(int bondDimensionRight);
   inline void SetLeftDimension(int bondDimensionLeft);

   //void ComputeDensityMatrixRight();
//   void ComputeDensityMatrixLeft();

   //void SymmetricUpdateOfTwoSites(AbstractMPSSite * leftSite , AbstractMPSSite * rightSite, RealVector * psi, RealDiagonalMatrix & singularValues );

   inline unsigned int GetBondDimensionLeft() const
     {return this->BondDimensionLeft;}
   inline unsigned int GetBondDimensionRight() const
     {return this->BondDimensionRight;}

   inline void SetLeftSite (AbstractMPSSite *  siteOnLeft) { this->SiteOnLeft = siteOnLeft; };
   inline void SetRightSite (AbstractMPSSite *  siteOnRight) { this->SiteOnRight = siteOnRight; };   
   
   inline long int GetVectorOneSiteIndice(int leftIndice, int rightIndice, int physicalIndice) const;
   inline long int GetVectorTwoSiteIndice(int leftIndice, int rightIndice, int physicalIndice) const;
   inline void DecodeVectorOneSiteIndice(long int vectorIndice, int& leftIndice, int & rightIndice, int& physicalIndice) const;
   
   inline void DecodeVectorTwoSiteIndice(long int vectorIndice, int& leftIndice, int & rightIndice, int& physicalIndice) const;
   
};




inline void AbstractMPSSite::SetRightDimension(int bondDimensionRight)
{
  this->BondDimensionRight = bondDimensionRight;
  this->SiteOnRight->BondDimensionLeft = bondDimensionRight;
}


inline void AbstractMPSSite::SetLeftDimension(int bondDimensionLeft)
{
  this->BondDimensionLeft = bondDimensionLeft;
  this->SiteOnLeft->BondDimensionRight = bondDimensionLeft;
}

inline void AbstractMPSSite::SetBondDimension(int bondDimensionLeft, int bondDimensionRight)
{
  this->BondDimensionLeft = bondDimensionLeft;
  this->BondDimensionRight = bondDimensionRight;
}


inline long int AbstractMPSSite::GetVectorOneSiteIndice(int leftIndice, int rightIndice, int physicalIndice) const
{
  return this->BondDimensionRight * (this->BondDimensionLeft * physicalIndice + leftIndice) + rightIndice;
}

inline long int AbstractMPSSite::GetVectorTwoSiteIndice(int leftIndice, int rightIndice, int physicalIndice) const
{
 return physicalIndice + this->SquarePhysicalDimension * (leftIndice + rightIndice * this->BondDimensionLeft);
}


inline void AbstractMPSSite::DecodeVectorOneSiteIndice(long int vectorIndice, int & leftIndice, int & rightIndice, int& physicalIndice) const
{
   rightIndice = vectorIndice%this->BondDimensionRight;
   vectorIndice/=this->BondDimensionRight;
   leftIndice =  vectorIndice % this->BondDimensionLeft;
   physicalIndice = vectorIndice / this->BondDimensionLeft;
}

inline void AbstractMPSSite::DecodeVectorTwoSiteIndice(long int vectorIndice, int& leftIndice, int & rightIndice, int& physicalIndice) const
{
   physicalIndice = vectorIndice%this->SquarePhysicalDimension;
   vectorIndice /= this->SquarePhysicalDimension;
   leftIndice =  vectorIndice % this->BondDimensionLeft;
   rightIndice = vectorIndice / this->BondDimensionLeft;
}


#endif
