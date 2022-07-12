////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//           class implementing generalized Halperin wave functions on the sphere          //
//                                                                            //
//                        last modification : 02/11/2007                      //
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
#include "ExtendedHalperinWavefunction.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//

ExtendedHalperinWavefunction::ExtendedHalperinWavefunction()
{
  this->NbrParticles = 0;  
  this->NbrParticlesPerLayer = 0;  
}


// constructor
//
// nbrParticles= number of particles per layer
// k = Jastrow factors between same species
// m = Jastrow factors between different species
// p = power inside Cauchy-determinant
// q = power inside Cauchy-permanent
// r = power of Cauchy-determinant
// s = power of Cauchy-permanent
// t = power of overall Pfaffian factor
// moveJastrowInside = flag to indicate whether Jastrow factors should be moved inside Cauchy determinant
ExtendedHalperinWavefunction::ExtendedHalperinWavefunction(int nbrParticles, int k, int m, int p, int q, int r, int s, int t, int u, int v, int b, bool moveJastrowInside )
{  
  if (nbrParticles&1)
    {
      cout << "Error: This implementation of Halperin Wavefunctions requires an even number of particles!"<<endl;
      exit(1);
    }
  if (nbrParticles<=0)
    {
      cout << "Error: The number of particles in ExtendedHalperinWavefunction has to be positive!"<<endl;
      exit(1);
    }
  this->NbrParticles = nbrParticles;
  this->NbrParticlesPerLayer = nbrParticles/2;
  this->P=p;
  this->Q=q;
  this->R=r;
  this->S=s;
  this->T=t;
  this->U=u;
  this->V=v;
  // cout << "P="<<P<<" Q="<<Q<<" R="<<R<<" S="<<S<<endl;
  if ( (moveJastrowInside) && (P!=0) && (R==1))
    {
      this->K_outside=(k%2);
      this->K_inside=(k/2)*2;
      this->M_outside=m%2;
      this->M_inside=(m/2)*2;
      this->JastrowInside=true;
    }
  else
    {
      this->K_outside=k;
      this->K_inside=0;
      this->M_outside=m;
      this->M_inside=0;
      this->JastrowInside=false;
    }
  if (((P!=0)||(JastrowInside))&&(R!=0)) HaveDeterminant=true;
  else  HaveDeterminant=false;
  
  this->DeterminantNorm=1.0;
  this->PermanentNorm=1.0;
  this->JastrowNorm=1.0;

#ifdef USE_LAPACK_CFCB
  this->Matrix = new ComplexLapackDeterminant(this->NbrParticlesPerLayer);
#else
  this->Matrix = new ComplexMatrix(this->NbrParticlesPerLayer,this->NbrParticlesPerLayer);
#endif
  this->Matrix2 = new ComplexMatrix(this->NbrParticlesPerLayer,this->NbrParticlesPerLayer);
  
  if (this->T!=0)
    this->PfaffianFactor = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  else
    this->PfaffianFactor = NULL;

  if (this->U!=0)
    this->PfaffianFactorIntraUp = new ComplexSkewSymmetricMatrix(this->NbrParticlesPerLayer);
  else
    this->PfaffianFactorIntraUp = NULL;

  if (this->V!=0)
    this->PfaffianFactorIntraDown = new ComplexSkewSymmetricMatrix(this->NbrParticlesPerLayer);
  else
    this->PfaffianFactorIntraDown = NULL;

  if (this->B!=0)
    this->PfaffianFactorInter = new ComplexSkewSymmetricMatrix(this->NbrParticlesPerLayer);
  else
    this->PfaffianFactorInter = NULL;

  
  this->J11 = new Complex[this->NbrParticlesPerLayer];
  this->J12 = new Complex[this->NbrParticlesPerLayer];
  this->J21 = new Complex[this->NbrParticlesPerLayer];
  this->J22 = new Complex[this->NbrParticlesPerLayer];
  
  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
  
  this->JastrowFactorElements = new Complex*[this->NbrParticles];
  for (int i=0; i< this->NbrParticles; ++i)
    JastrowFactorElements[i]= new Complex[NbrParticles];

}

// copy constructor
//
// function = reference on the wave function to copy

ExtendedHalperinWavefunction::ExtendedHalperinWavefunction(const ExtendedHalperinWavefunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrParticlesPerLayer = function.NbrParticlesPerLayer;
  this->K_outside=function.K_outside;
  this->K_inside=function.K_inside;
  this->M_outside=function.M_outside;
  this->M_inside=function.M_inside;
  this->P=function.P;
  this->Q=function.Q;
  this->R=function.R;
  this->S=function.S;
  this->T=function.T;
  this->U=function.U;  
  this->V=function.V;
  this->B=function.B;
  this->HaveDeterminant=function.HaveDeterminant;
  this->DeterminantNorm=function.DeterminantNorm;
  this->PermanentNorm=function.PermanentNorm;
  this->JastrowNorm=function.JastrowNorm;
  this->JastrowInside=function.JastrowInside;

#ifdef USE_LAPACK_CFCB
  this->Matrix = new ComplexLapackDeterminant(this->NbrParticlesPerLayer);
#else
  this->Matrix = new ComplexMatrix(this->NbrParticlesPerLayer,this->NbrParticlesPerLayer);
#endif
  this->Matrix2 = new ComplexMatrix(this->NbrParticlesPerLayer,this->NbrParticlesPerLayer);
  
  if (this->T!=0)
    this->PfaffianFactor = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  else
    this->PfaffianFactor = NULL;

  if (this->U!=0)
    this->PfaffianFactorIntraUp = new ComplexSkewSymmetricMatrix(this->NbrParticlesPerLayer);
  else
    this->PfaffianFactorIntraUp = NULL;

  if (this->V!=0)
    this->PfaffianFactorIntraDown = new ComplexSkewSymmetricMatrix(this->NbrParticlesPerLayer);
  else
    this->PfaffianFactorIntraDown = NULL;

  if (this->B!=0)
    this->PfaffianFactorInter = new ComplexSkewSymmetricMatrix(this->NbrParticlesPerLayer);
  else
    this->PfaffianFactorInter = NULL;


  this->J11 = new Complex[this->NbrParticlesPerLayer];
  this->J12 = new Complex[this->NbrParticlesPerLayer];
  this->J21 = new Complex[this->NbrParticlesPerLayer];
  this->J22 = new Complex[this->NbrParticlesPerLayer];
  
  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];

  this->JastrowFactorElements = new Complex*[this->NbrParticles];
  for (int i=0; i< this->NbrParticles; ++i)
    JastrowFactorElements[i]= new Complex[NbrParticles];

}

// destructor
//

ExtendedHalperinWavefunction::~ExtendedHalperinWavefunction()
{
  if (this->NbrParticlesPerLayer!=0)
    {
      for (int i=0; i< this->NbrParticles; ++i)
	delete [] JastrowFactorElements[i];
      delete [] JastrowFactorElements;
      delete [] SpinorUCoordinates;
      delete [] SpinorVCoordinates;
      delete [] J11;
      delete [] J12;
      delete [] J21;
      delete [] J22;
      delete Matrix;
      delete Matrix2;
      if (this->T!=0)
	delete [] PfaffianFactor;
      if (this->U!=0)
	delete [] PfaffianFactorIntraUp;
      if (this->V!=0)
	delete [] PfaffianFactorIntraDown;
      if (this->B!=0)
	delete [] PfaffianFactorInter;
    
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* ExtendedHalperinWavefunction::Clone ()
{
  return new ExtendedHalperinWavefunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex ExtendedHalperinWavefunction::operator()(RealVector& x)
{
  double s,c;
  // Convert particle positions into spinor form:
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
      this->SpinorUCoordinates[i].Im = this->SpinorUCoordinates[i].Re;
      this->SpinorUCoordinates[i].Re *= (c=cos(0.5 * x[1 + (i << 1)]));
      this->SpinorUCoordinates[i].Im *= -(s=sin(0.5 * x[1 + (i << 1)]));
      this->SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
      this->SpinorVCoordinates[i].Im = this->SpinorVCoordinates[i].Re;
      this->SpinorVCoordinates[i].Re *= c;
      this->SpinorVCoordinates[i].Im *= s;
      //cout << "U["<<i<<"]="<<SpinorUCoordinates[i]<<", "<< "V["<<i<<"]="<<SpinorVCoordinates[i]<<endl;
    }
  this->EvaluateTables();
  Complex result=1.0, tmp;
  // calculate Jastrow part:
  if (K_outside != 0)
    {
        tmp=1.0;
	for (int i=1; i<NbrParticlesPerLayer; ++i)
	  for (int j=0; j<i; ++j)
	    {
	      tmp *= JastrowNorm*JastrowFactorElements[i][j];
	      tmp *= JastrowNorm*JastrowFactorElements[i+NbrParticlesPerLayer][j+NbrParticlesPerLayer];
	    }
	result *= pow (tmp,K_outside);
    }
  if (M_outside != 0)
    {
      tmp=1.0;
      for (int i=0; i<NbrParticlesPerLayer; ++i)
	for (int j=NbrParticlesPerLayer; j<2*NbrParticlesPerLayer; ++j)
	  tmp *= JastrowNorm*JastrowFactorElements[i][j];
      result *= pow (tmp,M_outside);
    }

  // calculate Cauchy determinant
  if (HaveDeterminant)
    {      
      for (int i=0; i<NbrParticlesPerLayer; ++i)
	{
	  for (int j=NbrParticlesPerLayer; j<2*NbrParticlesPerLayer; ++j)
	    {
	      tmp=DeterminantNorm;
	      for (int p=P; p>0; --p) tmp*=JastrowFactorElements[i][j];
	      for (int p=P; p<0; ++p) tmp/=JastrowFactorElements[i][j];
	      for (int k=0;k<K_inside; k+=2)
		{
		  tmp*= J11[i];
		  tmp*= J22[j-NbrParticlesPerLayer];
		}
	      for (int m=0;m<M_inside; m+=2)
		{
		  tmp*= J12[i];
		  tmp*= J21[j-NbrParticlesPerLayer];
		}
	      //cout << " " <<tmp;
	      // initialize Cauchy determinant 
#ifdef USE_LAPACK_CFCB
	      Matrix->SetMatrixElement(i,j-NbrParticlesPerLayer,Real(tmp), Imag(tmp));
#else
	      (*Matrix)[i].Re(j-NbrParticlesPerLayer) = Real(tmp);
	      (*Matrix)[i].Im(j-NbrParticlesPerLayer) = Imag(tmp);
#endif
	    }
	  //cout << endl;
	}            
      tmp = Matrix->Determinant();
      //cout << "Calculated Determinant " << tmp << endl;      
      for (int r=this->R; r>0; --r) result*=tmp;
      for (int r=this->R; r<0; ++r) result/=tmp;
    }
    // calculate Cauchy Permanent
  if ((Q!=0)&&(S!=0))
    {
      for (int i=0; i<NbrParticlesPerLayer; ++i)
	{
	  for (int j=NbrParticlesPerLayer; j<2*NbrParticlesPerLayer; ++j)
	    {
	      tmp=PermanentNorm;
	      for (int q=Q; q>0; --q) tmp*=JastrowFactorElements[i][j];
	      for (int q=Q; q<0; ++q) tmp/=JastrowFactorElements[i][j];
	      // initialize Cauchy determinant 
	      (*Matrix2)[i].Re(j-NbrParticlesPerLayer) = Real(tmp);
	      (*Matrix2)[i].Im(j-NbrParticlesPerLayer) = Imag(tmp);
	    }
	}
      tmp = Matrix2->Permanent();
      for (int s=this->S; s>0; --s) result*=tmp;
      for (int s=this->S; s<0; ++s) result/=tmp;
    }
  //Calculate Pfaffian(s)
  if (T!=0)
    {
      for (int i=1; i<NbrParticles; ++i)
	for (int j=0; j<i; ++j)
	  PfaffianFactor->SetMatrixElement(i,j,0.5*M_PI/JastrowFactorElements[i][j]);

      tmp=PfaffianFactor->Pfaffian(); 	
      for (int t=this->T; t>0; --t) result*=tmp;
      for (int t=this->T; t<0; ++t) result/=tmp;
    }
  if (U!=0)
    {
     for (int i=1; i<NbrParticlesPerLayer; ++i)
	  for (int j=0; j<i; ++j)
	     PfaffianFactorIntraUp->SetMatrixElement(i,j,0.5*M_PI/JastrowFactorElements[i][j]);

      tmp=PfaffianFactorIntraUp->Pfaffian(); 	
      for (int u=this->U; u>0; --u) result*=tmp;
      for (int u=this->U; u<0; ++u) result/=tmp;
    }
  if (V!=0)
    {
     for (int i=NbrParticlesPerLayer+1; i<2*NbrParticlesPerLayer; ++i)
	  for (int j=NbrParticlesPerLayer; j<i; ++j)
	     PfaffianFactorIntraDown->SetMatrixElement(i-NbrParticlesPerLayer,j-NbrParticlesPerLayer,0.5*M_PI/JastrowFactorElements[i][j]);

      tmp=PfaffianFactorIntraDown->Pfaffian(); 	
      for (int v=this->V; v>0; --v) result*=tmp;
      for (int v=this->V; v<0; ++v) result/=tmp;
    }
  if (B!=0)
    {
     for (int i=0; i<NbrParticlesPerLayer; ++i)
	  for (int j=NbrParticlesPerLayer; j<2*NbrParticlesPerLayer; ++j)
	     PfaffianFactorInter->SetMatrixElement(i,j-NbrParticlesPerLayer,0.5*M_PI/JastrowFactorElements[i][j]);

      tmp=PfaffianFactorInter->Pfaffian(); 	
      for (int b=this->B; b>0; --b) result*=tmp;
      for (int b=this->B; b<0; ++b) result/=tmp;
    }

  return result;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex ExtendedHalperinWavefunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  // Import from spinors
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = uv.Re(2*i);
      this->SpinorUCoordinates[i].Im = uv.Im(2*i);
      this->SpinorVCoordinates[i].Re = uv.Re(2*i+1);
      this->SpinorVCoordinates[i].Im = uv.Im(2*i+1);
    }
  this->EvaluateTables();
  Complex result=1.0, tmp;
  // calculate Jastrow part:
  if (K_outside != 0)
    {
        tmp=1.0;
	for (int i=1; i<NbrParticlesPerLayer; ++i)
	  for (int j=0; j<i; ++j)
	    {
	      tmp *= JastrowNorm*JastrowFactorElements[i][j];
	      tmp *= JastrowNorm*JastrowFactorElements[i+NbrParticlesPerLayer][j+NbrParticlesPerLayer];
	    }
	result *= pow (tmp,K_outside);
    }
  if (M_outside != 0)
    {
      tmp=1.0;
      for (int i=0; i<NbrParticlesPerLayer; ++i)
	for (int j=NbrParticlesPerLayer; j<2*NbrParticlesPerLayer; ++j)
	  tmp *= JastrowNorm*JastrowFactorElements[i][j];
      result *= pow (tmp,M_outside);
    }

  // calculate Cauchy determinant
  if (HaveDeterminant)
    {      
      for (int i=0; i<NbrParticlesPerLayer; ++i)
	{
	  for (int j=NbrParticlesPerLayer; j<2*NbrParticlesPerLayer; ++j)
	    {
	      tmp=DeterminantNorm;
	      for (int p=P; p>0; --p) tmp*=JastrowFactorElements[i][j];
	      for (int p=P; p<0; ++p) tmp/=JastrowFactorElements[i][j];
	      for (int k=0;k<K_inside; k+=2)
		{
		  tmp*= J11[i];
		  tmp*= J22[j-NbrParticlesPerLayer];
		}
	      for (int m=0;m<M_inside; m+=2)
		{
		  tmp*= J12[i];
		  tmp*= J21[j-NbrParticlesPerLayer];
		}
	      //cout << " " <<tmp;
	      // initialize Cauchy determinant 
#ifdef USE_LAPACK_CFCB
	      Matrix->SetMatrixElement(i,j-NbrParticlesPerLayer,Real(tmp), Imag(tmp));
#else
	      (*Matrix)[i].Re(j-NbrParticlesPerLayer) = Real(tmp);
	      (*Matrix)[i].Im(j-NbrParticlesPerLayer) = Imag(tmp);
#endif
	    }
	  //cout << endl;
	}            
      tmp = Matrix->Determinant();
      //cout << "Calculated Determinant " << tmp << endl;      
      for (int r=this->R; r>0; --r) result*=tmp;
      for (int r=this->R; r<0; ++r) result/=tmp;
    }
    // calculate Cauchy Permanent
  if ((Q!=0)&&(S!=0))
    {
      for (int i=0; i<NbrParticlesPerLayer; ++i)
	{
	  for (int j=NbrParticlesPerLayer; j<2*NbrParticlesPerLayer; ++j)
	    {
	      tmp=PermanentNorm;
	      for (int q=Q; q>0; --q) tmp*=JastrowFactorElements[i][j];
	      for (int q=Q; q<0; ++q) tmp/=JastrowFactorElements[i][j];
	      // initialize Cauchy determinant 
	      (*Matrix2)[i].Re(j-NbrParticlesPerLayer) = Real(tmp);
	      (*Matrix2)[i].Im(j-NbrParticlesPerLayer) = Imag(tmp);
	    }
	}
      tmp = Matrix2->Permanent();
      for (int s=this->S; s>0; --s) result*=tmp;
      for (int s=this->S; s<0; ++s) result/=tmp;
    }
  //Calculate Pfaffian(s)
  if (T!=0)
    {
      for (int i=1; i<NbrParticles; ++i)
	for (int j=0; j<i; ++j)
	  PfaffianFactor->SetMatrixElement(i,j,0.5*M_PI/JastrowFactorElements[i][j]);

      tmp=PfaffianFactor->Pfaffian(); 	
      for (int t=this->T; t>0; --t) result*=tmp;
      for (int t=this->T; t<0; ++t) result/=tmp;
    }
  if (U!=0)
    {
     for (int i=1; i<NbrParticlesPerLayer; ++i)
	  for (int j=0; j<i; ++j)
	     PfaffianFactorIntraUp->SetMatrixElement(i,j,0.5*M_PI/JastrowFactorElements[i][j]);

      tmp=PfaffianFactorIntraUp->Pfaffian(); 	
      for (int u=this->U; u>0; --u) result*=tmp;
      for (int u=this->U; u<0; ++u) result/=tmp;
    }
  if (V!=0)
    {
     for (int i=NbrParticlesPerLayer+1; i<2*NbrParticlesPerLayer; ++i)
	  for (int j=NbrParticlesPerLayer; j<i; ++j)
	     PfaffianFactorIntraDown->SetMatrixElement(i-NbrParticlesPerLayer,j-NbrParticlesPerLayer,0.5*M_PI/JastrowFactorElements[i][j]);

      tmp=PfaffianFactorIntraDown->Pfaffian(); 	
      for (int v=this->V; v>0; --v) result*=tmp;
      for (int v=this->V; v<0; ++v) result/=tmp;
    }
  if (B!=0)
    {
     for (int i=0; i<NbrParticlesPerLayer; ++i)
	  for (int j=NbrParticlesPerLayer; j<2*NbrParticlesPerLayer; ++j)
	     PfaffianFactorInter->SetMatrixElement(i,j-NbrParticlesPerLayer,0.5*M_PI/JastrowFactorElements[i][j]);

      tmp=PfaffianFactorInter->Pfaffian(); 	
      for (int b=this->B; b>0; --b) result*=tmp;
      for (int b=this->B; b<0; ++b) result/=tmp;
    }

  return result;
}


// change the normalization of the funtion by a multiplicative factor
// factor = factor to be multiplied
void ExtendedHalperinWavefunction::Renormalize(double factor)
{
  int TotalJastrowTerms = NbrParticlesPerLayer*(K_outside*(NbrParticlesPerLayer-1)+M_outside*NbrParticlesPerLayer);
  this->JastrowNorm*= pow(factor,(double)1.0/TotalJastrowTerms);
}


// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void ExtendedHalperinWavefunction::AdaptNorm(RealVector& x)
{
  double det;
  int TotalJastrowTerms = NbrParticlesPerLayer*(K_outside*(NbrParticlesPerLayer-1)+M_outside*NbrParticlesPerLayer);
  bool actualHaveDeterminant = HaveDeterminant;
  int actualS = this->S;
  // switch off Cauchy-determinant and -permanent part
  this->HaveDeterminant=false;
  this->S=0; 
  det=Norm((*this)(x));
  while ((det<.1)||(det>50.0))
    {
      cout <<"Nj'="<< this->JastrowNorm << " det="<<det<<endl;
      if (det>1e300) 
	this->JastrowNorm*= pow((double)1.0e-300,(double)1.0/TotalJastrowTerms);
      else if (det==0.0) 
	this->JastrowNorm*= pow((double)1.0e300,(double)1.0/TotalJastrowTerms);
      else 
	this->JastrowNorm*= pow(det,(double)-1.0/TotalJastrowTerms);
      det=Norm((*this)(x));
      cout <<"Nj'="<< this->JastrowNorm << " det="<<det<<endl;
    }
  // switch determinant back on, and renormalize again:
  this->HaveDeterminant=actualHaveDeterminant;
  if (this->HaveDeterminant)
    {
      det=Norm((*this)(x));
      while ((det<.1)||(det>50.0))
	{
	  cout <<"Nd'="<< this->DeterminantNorm << " det="<<det<<endl;
	  if (det>1e300) 
	    this->DeterminantNorm*= pow((double)1.0e-300,1.0/((double)this->R*this->NbrParticlesPerLayer));
	  else if (det==0.0) 
	    this->DeterminantNorm*= pow((double)1.0e300,1.0/((double)this->R*this->NbrParticlesPerLayer));
	  else 
	    this->DeterminantNorm*= pow(det,-1.0/((double)this->R*this->NbrParticlesPerLayer));
	  det=Norm((*this)(x));
	  cout <<"Nd'="<< this->DeterminantNorm << " det="<<det<<endl;
	}
    }
  // switch determinant back on, and renormalize again:
  this->S=actualS;
  if ((this->Q)&&(this->S))
    {
      det=Norm((*this)(x));
      while ((det<.1)||(det>50.0))
	{
	  cout <<"Np'="<< this->PermanentNorm << " per="<<det<<endl;
	  if (det>1e300) 
	    this->PermanentNorm*= pow((double)1.0e-300,1.0/((double)this->S*this->NbrParticlesPerLayer));
	  else if (det==0.0) 
	    this->PermanentNorm*= pow((double)1.0e300,1.0/((double)this->S*this->NbrParticlesPerLayer));
	  else 
	    this->PermanentNorm*= pow(det,-1.0/((double)this->S*this->NbrParticlesPerLayer));
	  det=Norm((*this)(x));
	  cout <<"Np'="<< this->PermanentNorm << " per="<<det<<endl;
	}
    }
}


// normalize the wave-function over an average number of MC positions

void ExtendedHalperinWavefunction::AdaptAverageMCNorm(int thermalize, int average)
{
  int TotalJastrowTerms = NbrParticlesPerLayer*(K_outside*(NbrParticlesPerLayer-1)+M_outside*NbrParticlesPerLayer);
  ParticleOnSphereCollection * Particles = new ParticleOnSphereCollection(2*this->NbrParticlesPerLayer);
  cout << "AdaptNorm 1st time:"<<endl;
  this->AdaptNorm(Particles->GetPositions());
  Complex TmpMetropolis, TrialValue = (*this)(Particles->GetPositions());  
  double PreviousSamplingAmplitude = SqrNorm(TrialValue);
  double CurrentSamplingAmplitude = PreviousSamplingAmplitude;
  int NextCoordinates=0;
  // do some MC moves: accept or reject move according to probability |Psi_new|^2  / |Psi_old|^2
  for (int i = 0; i < thermalize; ++i)
    {
      Particles->Move(NextCoordinates);
      TmpMetropolis = (*this)(Particles->GetPositions());
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((Particles->GetRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      NextCoordinates = (int) (((double) NbrParticlesPerLayer) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticlesPerLayer) --NextCoordinates;      
    }
  cout << "AdaptNorm 2nd time:"<<endl;
  this->AdaptNorm(Particles->GetPositions());
  double SumTrialValues=0.0;
  for (int i = 0; i < average; ++i)
    {
      Particles->Move(NextCoordinates);
      TmpMetropolis = (*this)(Particles->GetPositions());
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((Particles->GetRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      NextCoordinates = (int) (((double) NbrParticlesPerLayer) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticlesPerLayer) --NextCoordinates;            
      SumTrialValues+=Norm(TmpMetropolis);
    }  
  if (TotalJastrowTerms!=0)
    this->JastrowNorm*= pow(SumTrialValues/average,(double)-1.0/TotalJastrowTerms);
  else
    {
      this->DeterminantNorm*= pow(SumTrialValues/average,(double)-1.0/this->NbrParticlesPerLayer);
    }
  TmpMetropolis = (*this)(Particles->GetPositions());
  cout << "Final adaptation: fct value = "<<TmpMetropolis<<endl;

  // testing:
  SumTrialValues=0.0;
  double SumSqrTrialValues=0.0;
  for (int i = 0; i < average; ++i)
    {
      Particles->Move(NextCoordinates);
      TmpMetropolis = (*this)(Particles->GetPositions());
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((Particles->GetRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      NextCoordinates = (int) (((double) NbrParticlesPerLayer) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticlesPerLayer) --NextCoordinates;            
      SumTrialValues+=Norm(TmpMetropolis);
      SumSqrTrialValues+=SqrNorm(TmpMetropolis);
    }


  TmpMetropolis = (*this)(Particles->GetPositions());
  cout << "Final adaptation: fct value = "<<TmpMetropolis<<endl;
  cout << "Final adaptation: ave value = "<<SumTrialValues/average<<endl;
  cout << "                  sqr value = "<<SumSqrTrialValues/average<<endl;

  
  delete Particles;
}


// this is the main part of the calculation of the paired wavefunction:
void ExtendedHalperinWavefunction::EvaluateTables()
{
  Complex Tmp;  

  for (int i=0;i<this->NbrParticles;++i)
    {
      for(int j=0;j<i;++j)
	{
	  Tmp = (this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]);
	  JastrowFactorElements[i][j] = Tmp;
	  JastrowFactorElements[j][i] = -Tmp;
	}
      JastrowFactorElements[i][i] = 0.0;
    }

  if (JastrowInside)
    {
      for (int i=0;i<this->NbrParticlesPerLayer;++i)
	{
	  J12[i]=1.0;
	  J21[i]=1.0;
	  for(int j=0;j<this->NbrParticlesPerLayer;++j)
	    {
	      J12[i] *= JastrowFactorElements[i][j+this->NbrParticlesPerLayer];
	      J21[i] *= -JastrowFactorElements[j+this->NbrParticlesPerLayer][i];
	    }
	  J11[i]=1.0;
	  J22[i]=1.0;
	  for(int j=0;j<i;j++)
	    {
	      J11[i] *= JastrowFactorElements[i][j];
	      J22[i] *= JastrowFactorElements[i+this->NbrParticlesPerLayer][j+this->NbrParticlesPerLayer];
	    }
	  for(int j=i+1;j<NbrParticlesPerLayer;j++)
	    {
	      J11[i] *= JastrowFactorElements[i][j];
	      J22[i] *= JastrowFactorElements[i+this->NbrParticlesPerLayer][j+this->NbrParticlesPerLayer];
	    }
	}
    }

}
