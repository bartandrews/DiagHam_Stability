#include "AbstractMPOperatorOBC.h"
#include "Tensor/Tensor3.h"
#include "Matrix/RealSymmetricMatrix.h"
#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;

AbstractMPOperatorOBC::AbstractMPOperatorOBC()
{
  this->NbrNonZeroElements = 0;
  this->IndexValues = 0;
  this->PhysicalDimension = 0;
  this->MPOBondDimension = 0;
}

AbstractMPOperatorOBC::~AbstractMPOperatorOBC()
{
}


// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractMPOperatorOBC::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->HilbertSpace = hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractMPOperatorOBC::GetHilbertSpace ()
{
  return this->HilbertSpace;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractMPOperatorOBC::GetHilbertSpaceDimension ()
{
  if(IDMRGFlag)
    return  this->SiteRight->GetBondDimensionRight()* this->SiteLeft->GetBondDimensionLeft()*this->PhysicalDimension *this->PhysicalDimension;
  return  this->Site->GetBondDimensionRight()* this->Site->GetBondDimensionLeft()*this->PhysicalDimension;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractMPOperatorOBC::GetTwoSitesHilbertSpaceDimension()
{
  return  this->SiteRight->GetBondDimensionRight()* this->SiteLeft->GetBondDimensionLeft()*this->PhysicalDimension *this->PhysicalDimension;
}
 
// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractMPOperatorOBC::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// set site to be acted on
//
// site = pointer to the siteto use 

void AbstractMPOperatorOBC::SetSite (AbstractMPSSite* site)
{
 this->Site = site;
}


// set site to be acted on
//
// site = pointer to the siteto use 

void AbstractMPOperatorOBC::SetSiteLeftAndRight (AbstractMPSSite* siteLeft,AbstractMPSSite* siteRight)
{
 this->SiteLeft = siteLeft;
 this->SiteRight = siteRight;
}




// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored
RealVector& AbstractMPOperatorOBC::LowLevelMultiply(RealVector& vSource, RealVector& vDestination)
{ 
  return this->LowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored
ComplexVector& AbstractMPOperatorOBC::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{ 
//cout <<"ComplexVector& AbstractMPOperatorOBC::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)"<<endl;
  return this->LowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}



// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored
RealVector& AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelMultiplyTwoSites(vSource, vDestination, 0, this->GetTwoSitesHilbertSpaceDimension());
}


// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored
ComplexVector& AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelMultiplyTwoSites(vSource, vDestination, 0, this->GetTwoSitesHilbertSpaceDimension());
}



// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

 RealSymmetricMatrix& AbstractMPOperatorOBC::GetTwoSitesHamiltonian (RealSymmetricMatrix & M)
{
cout <<"using low-level undefined  RealSymmetricMatrix& AbstractMPOperatorOBC::GetTwoSitesHamiltonian (RealSymmetricMatrix & M)" <<endl;
  return M;  
}

// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

HermitianMatrix& AbstractMPOperatorOBC::GetTwoSitesHamiltonian (HermitianMatrix & M)
{
cout <<"using low-level undefined  HermitianMatrix& AbstractMPOperatorOBC::GetTwoSitesHamiltonian (HermitianMatrix & M)" <<endl;
  return M;  
}



void AbstractMPOperatorOBC::ComputeL(Tensor3<double> & L)
{
  cout <<"using low-level undefined  void AbstractMPOperatorOBC::ComputeL(Tensor3<double> & L)" <<endl;
}


void AbstractMPOperatorOBC::ComputeR(Tensor3<double> & R)
{
 cout <<"using low-level undefined  void AbstractMPOperatorOBC::ComputeR(Tensor3<double> & R)" <<endl;
}




RealVector& AbstractMPOperatorOBC::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent)
{
 if(this->IDMRGFlag)
  return this->LowLevelMultiplyTwoSites(vSource,vDestination, firstComponent, nbrComponent); 
else 
 return this->LowLevelMultiplyOneSite(vSource,vDestination, firstComponent, nbrComponent);
}

ComplexVector& AbstractMPOperatorOBC::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				       int firstComponent, int nbrComponent)
{
//cout <<"ComplexVector& AbstractMPOperatorOBC::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination,     int firstComponent, int nbrComponent)"<<endl;
 if(this->IDMRGFlag)
  return this->LowLevelMultiplyTwoSites(vSource,vDestination, firstComponent, nbrComponent); 
else 
 return this->LowLevelMultiplyOneSite(vSource,vDestination, firstComponent, nbrComponent);
}



void AbstractMPOperatorOBC::ComputeL(Tensor3<Complex> & L)
{
  cout <<"using low-level undefined  void AbstractMPOperatorOBC::ComputeL(Tensor3<Complex> & L)" <<endl;
}


void AbstractMPOperatorOBC::ComputeR(Tensor3<Complex> & R)
{
 cout <<"using low-level undefined  void AbstractMPOperatorOBC::ComputeR(Tensor3<Complex> & R)" <<endl;
}




// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractMPOperatorOBC::LowLevelMultiplyOneSite(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent)
{
   cout <<"using low-level undefined RealVector& AbstractMPOperatorOBC::LowLevelMultiplyOneSite(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)" <<endl;
 return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractMPOperatorOBC::LowLevelMultiplyOneSite(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
   cout <<"using low-level undefined ComplexVector& AbstractMPOperatorOBC::LowLevelMultiplyOneSite(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)" <<endl;
 return vDestination;
}




// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination,  int firstComponent, int nbrComponent)
{
 cout <<"using low-level undefined RealVector& AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)" <<endl;

 return vDestination;
}

ComplexVector& AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(ComplexVector& vSource, ComplexVector& vDestination,  int firstComponent, int nbrComponent)
{
 cout <<"using low-level undefined ComplexVector& AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(ComplexVector& vSource, ComplexVector& vDestination,  int firstComponent, int nbrComponent)" <<endl;

 return vDestination;
}
