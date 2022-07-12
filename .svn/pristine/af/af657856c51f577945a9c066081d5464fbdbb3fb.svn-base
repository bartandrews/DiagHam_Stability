
#ifndef _TENSOR3_H
#define _TENSOR3_H
#include "config.h"

#include<iostream>
#include "GeneralTools/GarbageFlag.h"
#include "Matrix/Matrix.h"
#include "MathTools/Complex.h"
using std::cout;
using std::endl;

class Matrix;

template <class T> class Tensor3
{
 private :
  unsigned int FirstDimension;
  unsigned int SecondDimension;
  unsigned int ThirdDimension;
  T * TensorElements;
  GarbageFlag Flag;
  
 public :
  
  Tensor3();
  Tensor3(int firstDimension,int secondDimension, int thirdDimension); 
  Tensor3(int firstDimension,int secondDimension, int thirdDimension,bool initiateFlag);
  // copy constructor (without duplicating datas)
  //
  Tensor3(const Tensor3 & tensor3); 
  ~Tensor3();
  
  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  Tensor3 & operator = (const Tensor3 & tensor3);
  
  void PrintTensor();
  
  T & operator()(unsigned int index1, unsigned int index2, unsigned int index3);
  T & operator[] (unsigned int index);  
  template <int X,int Y> inline void Contract(Tensor3<T> * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent);

  virtual void ContractWithMatrixOnThirdAndFirstIndices(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent);
  virtual void ContractWithMatrixOnThirdAndSecondIndices(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent);
  virtual void ContractWithMatrixOnSecondAndFirstIndices(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent);
  virtual void ContractWithMatrixOnSecondAndSecondIndices(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent);
  virtual void ContractWithMatrixOnFirstAndFirstIndices(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent);
  virtual void ContractWithMatrixOnFirstAndSecondIndices(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent);
};

template <class T>
Tensor3<T>::Tensor3()
{
  this->FirstDimension = 0;
  this->SecondDimension = 0;
  this->ThirdDimension = 0;
  this->TensorElements = 0;
}

template <class T>
Tensor3<T>::Tensor3(int firstDimension,int secondDimension, int thirdDimension)
{
  this->FirstDimension = firstDimension;
  this->SecondDimension = secondDimension;
  this->ThirdDimension = thirdDimension;
  this->TensorElements = new T [this->FirstDimension* this->SecondDimension* this->ThirdDimension];
  this->Flag.Initialize();
}

template <class T>
Tensor3<T>::Tensor3(int firstDimension,int secondDimension, int thirdDimension, bool initiateFlag)
{
  this->FirstDimension = firstDimension;
  this->SecondDimension = secondDimension;
  this->ThirdDimension = thirdDimension;
  this->TensorElements = new T [this->FirstDimension * this->SecondDimension * this->ThirdDimension];
  this->Flag.Initialize();
  if (initiateFlag)
    {
      for (int i = 0; i< this->FirstDimension* this->SecondDimension* this->ThirdDimension ; i++)
	this->TensorElements[i] = 0;
    }
}


template <class T>
Tensor3<T>::~Tensor3()
{
 if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
     delete [] this->TensorElements;
}

template <class T>
Tensor3<T>::Tensor3(const Tensor3<T> & tensor3)
{
  this->FirstDimension = tensor3.FirstDimension;
  this->SecondDimension = tensor3.SecondDimension;
  this->ThirdDimension = tensor3.ThirdDimension;
  this->TensorElements = tensor3.TensorElements;
  this->Flag = tensor3.Flag;
}


// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

template <class T>
Tensor3<T> & Tensor3<T>::operator = (const Tensor3<T> & tensor3) 
{
  this->FirstDimension = tensor3.FirstDimension;
  this->SecondDimension = tensor3.SecondDimension;
  this->ThirdDimension = tensor3.ThirdDimension;
  this->TensorElements = tensor3.TensorElements;
  this->Flag = tensor3.Flag;
  return *this;
}

template <class T>
void Tensor3<T>::ContractWithMatrixOnThirdAndFirstIndices(Tensor3<T> * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
    for(int i=firstComponent; i< firstComponent + nbrComponent ; i++)
    {
      for(int k=0; k< this->SecondDimension; k++)
	{
	  for(int l=0; l< sourceMatrix->GetNbrColumn() ; l++)
	    {
           T & Tmp = (*this)(i,k,l);
           T Tmp1;
	  for(int p=0; p< this->ThirdDimension; p++)
	    {
       sourceMatrix->GetMatrixElement(p, l, Tmp1);
       Tmp += (*sourceTensor)(i,k,p) * Tmp1;
}
}
}
}
}  

template <class T>
void Tensor3<T>::ContractWithMatrixOnThirdAndSecondIndices(Tensor3<T> * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
    for(int i=firstComponent; i< firstComponent + nbrComponent ; i++)
    {
      for(int k=0; k< this->SecondDimension; k++)
	{
	  for(int l=0; l< sourceMatrix->GetNbrRow() ; l++)
	    {
           T & Tmp = (*this)(i,k,l);
           T Tmp1;
	  for(int p=0; p< this->ThirdDimension; p++)
	    {
       sourceMatrix->GetMatrixElement(p, l, Tmp1);
       Tmp += (*sourceTensor)(i,k,p) * Tmp1;
}
}
}
}
}  


template <class T>
void Tensor3<T>::ContractWithMatrixOnSecondAndFirstIndices(Tensor3<T> * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
    for(int i=firstComponent; i< firstComponent + nbrComponent ; i++)
    {
 for(int p=0; p< this->ThirdDimension; p++)
	    {
	  for(int l=0; l< sourceMatrix->GetNbrColumn() ; l++)
	    {
           T & Tmp = (*this)(i,l,p);
           T Tmp1;
	       for(int k=0; k< this->SecondDimension; k++)
	{
       sourceMatrix->GetMatrixElement(k, l, Tmp1);
       Tmp += (*sourceTensor)(i,k,p) * Tmp1;
}
}
}
}
}  

template <class T>
void Tensor3<T>::ContractWithMatrixOnSecondAndSecondIndices(Tensor3<T> * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
    T Tmp1;
    for(int i=firstComponent; i< firstComponent + nbrComponent ; i++)
    {

     for(int p=0; p< this->ThirdDimension; p++)
     {

	  for(int l=0; l< sourceMatrix->GetNbrRow() ; l++)
	    {

           T & Tmp = (*this)(i,l,p);
           T Tmp1;
      for(int k=0; k< this->SecondDimension; k++)
	{
       sourceMatrix->GetMatrixElement(l, k, Tmp1);
       Tmp += (*sourceTensor)(i,k,p) * Tmp1;
}
}
}
}
}  


template <class T>
void Tensor3<T>::ContractWithMatrixOnFirstAndFirstIndices(Tensor3<T> * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
      T Tmp1;
      for(int k=0; k< this->SecondDimension; k++)
	{
       for(int p=0; p< this->ThirdDimension; p++)
	{
	  for(int l=firstComponent; l< firstComponent + nbrComponent; l++)
	    {
           T & Tmp = (*this)(l,k,p);
  
    for(int i=0; i< sourceMatrix->GetNbrRow(); i++)
    {
           sourceMatrix->GetMatrixElement(i, l, Tmp1);
           Tmp += (*sourceTensor)(i,k,p) * Tmp1;
    }
  }
}
}
}  

template <class T>
void Tensor3<T>::ContractWithMatrixOnFirstAndSecondIndices(Tensor3<T> * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
     T Tmp1;
      for(int k=0; k< this->SecondDimension; k++)
	{
       for(int p=0; p< this->ThirdDimension; p++)
	{
	  for(int l=0; l< sourceMatrix->GetNbrRow() ; l++)
	    {
           T & Tmp = (*this)(l,k,p);
  
    for(int i=firstComponent; i< firstComponent + nbrComponent ; i++)
    {
           sourceMatrix->GetMatrixElement(l, i, Tmp1);
           Tmp += (*sourceTensor)(i,k,p) * Tmp1;
    }
  }
}
}

}



template <class T>
inline T & Tensor3<T>::operator [] (unsigned int index) 
{ 
  return TensorElements[index];
} 


template <class T>
inline T & Tensor3<T>::operator() (unsigned int index1, unsigned int index2, unsigned int index3) 
{ 
  return TensorElements[index1 + index2*this->FirstDimension+this->FirstDimension*this->SecondDimension*index3];
} 

template <class T>
void Tensor3<T>::PrintTensor()
{
  for(int i=0; i< this->FirstDimension; i++)
    {
      for(int k=0; k< this->SecondDimension; k++)
	{
	  for(int p=0; p< this->ThirdDimension; p++)
	    {
           if ((*this)(i,k,p) != 0)
	      cout <<i << " "<<k<< " "<<p<<" " << (*this)(i,k,p) << endl;
	    }
	}
    }
}

template <class T> template <int X,int Y>
inline void Tensor3<T>::Contract(Tensor3<T> * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
cout <<"impossible Contraction"<<endl;
}

template <> template <>
inline void Tensor3<double>::Contract<0,0>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnFirstAndFirstIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}

template <> template <>
inline void Tensor3<double>::Contract<0,1>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnFirstAndSecondIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}


template <> template <>
inline void Tensor3<double>::Contract<1,0>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnSecondAndFirstIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}

template <> template <>
inline void Tensor3<double>::Contract<1,1>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnSecondAndSecondIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}


template <> template <>
inline void Tensor3<double>::Contract<2,0>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnThirdAndFirstIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}

template <> template <>
inline void Tensor3<double>::Contract<2,1>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnThirdAndSecondIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}

template <> template <>
inline void Tensor3<Complex>::Contract<0,0>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnFirstAndFirstIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}

template <> template <>
inline void Tensor3<Complex>::Contract<0,1>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnFirstAndSecondIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}


template <> template <>
inline void Tensor3<Complex>::Contract<1,0>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnSecondAndFirstIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}

template <> template <>
inline void Tensor3<Complex>::Contract<1,1>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnSecondAndSecondIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}


template <> template <>
inline void Tensor3<Complex>::Contract<2,0>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnThirdAndFirstIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}

template <> template <>
inline void Tensor3<Complex>::Contract<2,1>(Tensor3 * sourceTensor,Matrix * sourceMatrix, int  firstComponent, int nbrComponent)
{
  ContractWithMatrixOnThirdAndSecondIndices(sourceTensor,sourceMatrix,firstComponent,nbrComponent);
}


#endif
