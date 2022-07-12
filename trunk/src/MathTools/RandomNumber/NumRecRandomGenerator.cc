#include "NumRecRandomGenerator.h"

#include <cmath>
#include <cstdio>

// defines valid for the whole file:
#define EPS 1.2e-14
#define RNMX (1.0-EPS)


// defines for GetRealRandomNumber()
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NDIV (1+(IM-1)/NTAB)



using std::log;
using std::sqrt;


/* ###################################### global variable definitions ###################################### */ 

NumRecRandomGenerator::NumRecRandomGenerator(long seed)
{
  this->idum=-1;
  this->iy=0;
  this->iset=0;
  this->GetRealRandomNumber(); //get one random number -> if array iv was not initialized, now it will be
  this->NbrGeneratedNumbers = 1ul;
  if (seed!=-1) this->SetSeed(seed);  
}


NumRecRandomGenerator::NumRecRandomGenerator(const NumRecRandomGenerator& generator)
{  
  this->idum=generator.idum;
  this->iy=generator.iy;
  this->iset=generator.iset;
  for (int i=0; i<NTAB; ++i)
    this->iv[i]=generator.iv[i];
  this->gset=generator.gset;
  this->NbrGeneratedNumbers = generator.NbrGeneratedNumbers;
}

NumRecRandomGenerator::~NumRecRandomGenerator()
{
}


// clone random number generator 
//
// return value = clone of the random number generator

AbstractRandomNumberGenerator* NumRecRandomGenerator::Clone ()
{
  return new NumRecRandomGenerator(*this);
}


void NumRecRandomGenerator::SetSeed(const unsigned long& seed) 
{
  double tmp;
  this->idum = seed;
  for (int i=0;i<NTAB;i++)
    tmp=this->GetRealRandomNumber();        //fill the whole tab with new random numbers;
}

// get real random number between 0 and 1
//
// return value = random number
double NumRecRandomGenerator::GetRealRandomNumber()
{
  ++this->NbrGeneratedNumbers;
  int j;
  long k;
  double temp;
  if (idum<= 0 || iy==0) {
    if  (-(idum) < 1)  idum=1;
    else idum = -(idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(idum)/IQ;
      idum=IA*(idum-k*IQ)-IR*k;
      if (idum<0) idum += IM;
      if (j < NTAB) iv[j]= idum;
    }
    iy=iv[0];
  }
  k=(idum)/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum <0)  idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]= idum;
  if ((temp=AM*iy) > RNMX ) return RNMX;
  else return temp;
}


