////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of particle on torus density operator                //
//                                                                            //
//                        last modification : 11/12/2002                      //
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


#include "config.h"
#include "Operator/ParticleOnTorusDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
using std::cout;
using std::endl;


// constructor from default datas
//
// particle = hilbert space associated to the particles
// index = index of the density operator

ParticleOnTorusDensityOperator::ParticleOnTorusDensityOperator(ParticleOnTorus* particle, int index, int qx, double ratio)
: ParticleOnSphereDensityOperator(particle, index)
{
  this->Particle= (ParticleOnTorus*) (particle->Clone());
  this->OperatorIndex = index;
  this->OperatorIndexDagger = index;
  this->Ratio = ratio;
  this->Qx = qx;
  this->CalculateFormFactor();
  if (fabs(Imag(FormFactor))> 1e-12)
    cout << "Form factor in ParticleOnTorusDensityOperator is complex"<<std::endl;
  else
    cout << "Form factor in ParticleOnTorusDensityOperator is real"<<std::endl;
  }

// constructor when dealing with two different Hilbert spaces
//
// particle = hilbert space associated to the right hand state (target space has to be fixed to the hilbert space associated to the left hand state)
// indexDagger = index of the creation operator that is part of the density operator
// index = index of the annihilation operator that is part of the density operator
 
ParticleOnTorusDensityOperator::ParticleOnTorusDensityOperator(ParticleOnTorus* particle, int indexDagger, int index, int qx, double ratio)
: ParticleOnSphereDensityOperator(particle, indexDagger, index)
{
  this->Particle= (ParticleOnTorus*) (particle->Clone());
  this->OperatorIndexDagger = indexDagger;
  this->OperatorIndex = index;
  this->Ratio = ratio;
  this->Qx = qx;
  this->CalculateFormFactor();
  if (fabs(Imag(FormFactor))> 1e-12)
    cout << "Form factor in ParticleOnTorusDensityOperator is complex"<<std::endl;
  else
    cout << "Form factor in ParticleOnTorusDensityOperator is real"<<std::endl;
}

// copy constructor
//
// oper = operator to copy
  
ParticleOnTorusDensityOperator::ParticleOnTorusDensityOperator(ParticleOnTorusDensityOperator& oper)
: ParticleOnSphereDensityOperator(oper)
{
  this->Particle = (ParticleOnTorus*) (oper.Particle->Clone());
  this->OperatorIndexDagger = oper.OperatorIndexDagger;
  this->OperatorIndex = oper.OperatorIndex;
  this->Ratio = oper.Ratio;
  this->Qx = oper.Qx;
  this->FormFactor = oper.FormFactor;
}

// destructor
//

ParticleOnTorusDensityOperator::~ParticleOnTorusDensityOperator()
{
  delete this->Particle;
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnTorusDensityOperator::Clone ()
{
  return new ParticleOnTorusDensityOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnTorus*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnTorusDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnTorusDensityOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element
//
// no form factors implemented

Complex ParticleOnTorusDensityOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  double Element = 0.0;
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Dim = firstComponent + nbrComponent;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (int i = firstComponent; i < Dim; ++i)
	    {
	      Element += V1[i] * V2[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (int i = firstComponent; i < Dim; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		Element += V1[TmpIndex] * V2[i] * TmpCoefficient;
	    }
	}
    }
  else
    {
      long Dim = firstComponent + nbrComponent;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (long i = firstComponent; i < Dim; ++i)
	    {
	      Element += V1[i] * V2[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (long i = firstComponent; i < Dim; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		Element += V1[TmpIndex] * V2[i] * TmpCoefficient;
	    }
	}
     }
  return Complex(this->FormFactor*Element);
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnTorusDensityOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								 int firstComponent, int nbrComponent)
{
  if (fabs(Imag(FormFactor))> 1e-12)
  {
    std::cerr << "Error: calling real multiplication routine in ParticleOnTorusDensityOperator for a complex valued operator"<<endl;
    std::exit(1);
  }
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Last = firstComponent + nbrComponent;;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      vDestination[i] += vSource[i] * Norm(this->FormFactor) * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		vDestination[TmpIndex] +=  vSource[i] * TmpCoefficient * Norm(this->FormFactor);
	    }
	}
    }
  else
    {
      long Last = ((long) firstComponent) + ((long) nbrComponent);
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      vDestination[i] += vSource[i] * Norm(this->FormFactor) * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  long TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		vDestination[TmpIndex] +=  vSource[i] * TmpCoefficient * Norm(this->FormFactor);
	    }
	}
    }
  return vDestination;
}
  
// multiply a set of vectors by the current operator for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* ParticleOnTorusDensityOperator::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									 int firstComponent, int nbrComponent)
{
  if (fabs(Imag(FormFactor))> 1e-12)
  {
    std::cerr << "Error: calling real multiplication routine in ParticleOnTorusDensityOperator for a complex valued operator"<<endl;
    std::exit(1);
  }
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Last = firstComponent + nbrComponent;;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      double TmpCoefficient = Norm(this->FormFactor) * this->Particle->AdA(i, this->OperatorIndex);
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][i] += vSources[k][i] * TmpCoefficient;
		}
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		{
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][TmpIndex] +=  vSources[k][i] * TmpCoefficient * Norm(this->FormFactor);
		    }
		}
	    }
	}
    }
  else
    {
      long Last = ((long) firstComponent) + ((long) nbrComponent);
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      double TmpCoefficient = Norm(this->FormFactor) * this->Particle->AdA(i, this->OperatorIndex);
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][i] += vSources[k][i] * TmpCoefficient;
		}
	    }
	}
      else
	{
	  long TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		{
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][TmpIndex] +=  vSources[k][i] * TmpCoefficient * Norm(this->FormFactor);
		    }
		}
	    }
	}
    }
  return vDestinations;
}


// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element
//
// no form factors implemented

Complex ParticleOnTorusDensityOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  Complex Element = 0.0;
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Dim = firstComponent + nbrComponent;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (int i = firstComponent; i < Dim; ++i)
	    {
	      Element += Conj(V1[i]) * V2[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  int FullDim = this->Particle->GetTargetHilbertSpaceDimension();
	  for (int i = firstComponent; i < Dim; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if ((TmpIndex != FullDim) && (TmpCoefficient != 0.0))
		Element += Conj(V1[TmpIndex]) * V2[i] * TmpCoefficient;
	    }
	}
    }
  else
    {
      long Dim = firstComponent + nbrComponent;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (long i = firstComponent; i < Dim; ++i)
	    {
	      Element += Conj(V1[i]) * V2[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  long TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (long i = firstComponent; i < Dim; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		Element += Conj(V1[TmpIndex]) * V2[i] * TmpCoefficient;
	    }
	}
     }
  return Element*FormFactor;
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnTorusDensityOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								    int firstComponent, int nbrComponent)
{
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  int Last = firstComponent + nbrComponent;;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      vDestination[i] += vSource[i] * this->FormFactor * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  int Last = firstComponent + nbrComponent;;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		vDestination[TmpIndex] += vSource[i] * TmpCoefficient * FormFactor;
	    }
	}
    }
  else
    {
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  long Last = ((long) firstComponent) + ((long) nbrComponent);
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      vDestination[i] += vSource[i] * this->FormFactor * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  long TmpIndex;
	  double TmpCoefficient = 0.0;
	  long Last = ((long) firstComponent) + ((long) nbrComponent);
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		vDestination[TmpIndex] += vSource[i] * TmpCoefficient * FormFactor;
	    }
	}
    }
  return vDestination;
}
  
// multiply a set of vectors by the current operator for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ParticleOnTorusDensityOperator::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									    int firstComponent, int nbrComponent)
{
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Last = firstComponent + nbrComponent;;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      Complex TmpCoefficient = this->FormFactor * this->Particle->AdA(i, this->OperatorIndex);
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][i] += vSources[k][i] * TmpCoefficient;
		}
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		{
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][TmpIndex] +=  vSources[k][i] * TmpCoefficient * FormFactor;
		    }
		}
	    }
	}
    }
  else
    {
      long Last = ((long) firstComponent) + ((long) nbrComponent);
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      Complex TmpCoefficient = this->FormFactor * this->Particle->AdA(i, this->OperatorIndex);
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][i] += vSources[k][i] * TmpCoefficient;
		}
	    }
	}
      else
	{
	  long TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		{
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][TmpIndex] +=  vSources[k][i] * TmpCoefficient * FormFactor;
		    }
		}
	    }
	}
    }
  return vDestinations;
}

// calculate the form factor for this fourier component
//
// return value = form factor

void ParticleOnTorusDensityOperator::CalculateFormFactor()
{  
  int nbrFluxQuanta=this->Particle->GetNbrOrbitals();
  
  double Lx=sqrt(2.0*M_PI*nbrFluxQuanta/this->Ratio);
  const double invLx = 1.0 / Lx;
  double Ly=sqrt(2.0*M_PI*nbrFluxQuanta*this->Ratio);
  const double invLy = 1.0 / Ly;
  
  const double kxFactor=2.0*M_PI*invLx;
  const double kyFactor=2.0*M_PI*invLy;
  const double qxFactor=kxFactor;
  
  int ky = this->OperatorIndex;
  int qx = this->Qx;
  int qy = (this->OperatorIndexDagger - this->OperatorIndex);
  
  int m1=(int)(-qy*kyFactor/Lx);
  int m2=m1+1;
  Complex numerator=0.0, denominator=0.0, exp1Plus=0, exp1Minus=0, exp2=0.0;
  
  while (exp1Plus>-34.5 || exp1Minus>-34.5) //while terms are machine significant
  {
    exp1Plus=(Complex(0.5*(m1*Lx+kyFactor*ky+kyFactor*(ky+qy)),-qxFactor*qx)*Complex(0.5*(m1*Lx+kyFactor*ky+kyFactor*(ky+qy)),-qxFactor*qx)
	      -0.5*(m1*Lx+kyFactor*(ky+qy))*(m1*Lx+kyFactor*(ky+qy))-0.5*kyFactor*ky*kyFactor*ky); //numerator exponent with positive m
    exp1Minus=(Complex(0.5*(m2*Lx+kyFactor*ky+kyFactor*(ky+qy)),-qxFactor*qx)*Complex(0.5*(m2*Lx+kyFactor*ky+kyFactor*(ky+qy)),-qxFactor*qx)
	      -0.5*(m2*Lx+kyFactor*(ky+qy))*(m2*Lx+kyFactor*(ky+qy))-0.5*kyFactor*ky*kyFactor*ky); //numerator exponent with negative m

    numerator+=exp(exp1Plus)+exp(exp1Minus); //sum symmetrically out from the origin
    
    --m1;
    ++m2;
    
  }
  
  m2=0;
  while (exp2>-34.5) //while terms are machine significant
  {
    exp2=-0.25*m2*m2*Lx*Lx; //denominator exponent
    
    if (m2==0)
      denominator+=exp(exp2); //only one term at the origin
    else
      denominator+=2.0*exp(exp2); //sum symmetrically out from the origin
     
    m2++;
  }
  
  this->FormFactor = numerator/denominator;
  
}


