#include "SimplePermutations.h"
#include "FactorialCoefficient.h"

// default constructor
// nbrElements = number of elements
SimplePermutations::SimplePermutations(int nbrElements)
{
  this->NbrElements=nbrElements;
  FactorialCoefficient Number;
  Number.SetToOne();
  Number.FactorialMultiply(nbrElements);
  this->NbrPermutations = Number.GetIntegerValue();
  this->MyPermutation=new unsigned[nbrElements];
  this->Permutations=new unsigned*[NbrPermutations];
  for (int i=0; i<NbrPermutations; ++i)
    this->Permutations[i]=new unsigned[nbrElements];
  CurrentPermutation=0;
  for (int i=0; i<nbrElements; ++i)
    this->MyPermutation[i]=i;
  this->Permute(0,nbrElements,1);
  this->CleanUp=true;
}

// destructor
SimplePermutations::~SimplePermutations()
{
  if (CleanUp)
    {
      for (int i=0; i<NbrPermutations; ++i)
	delete [] this->Permutations[i];
      delete [] this->Permutations;
    }
  delete [] this->MyPermutation;
}


// check out Permutations to an external method
// if this method is called, the caller is responsible for cleaning up
// if necessary, a new copy is created
unsigned** SimplePermutations::CheckOutPermutations()
{
  if (this->CleanUp == true)
    {
      this->CleanUp = false;
      return this->Permutations;
    }
  else
    {
      unsigned **Result = new unsigned*[NbrPermutations];
      for (int i=0; i<NbrPermutations; ++i)
	{
	  Result[i]=new unsigned[NbrElements];
	  for (int j=0; j<NbrElements; ++j)
	    Result[i][j] = this->Permutations[i][j];
	}
      return Result;
    }
}




void SimplePermutations::Swap(const int i, const int j, int &sign)
{
  //  cout << "Swap ("<<i<<", "<<j<<" sign: " << sign << endl;
  unsigned t;
  t = MyPermutation[i];
  MyPermutation[i] = MyPermutation[j];
  MyPermutation[j] = t;
  sign ^= 1;
} // Swap


void SimplePermutations::RotateLeft(const int start, const int n, int &sign)
{
  //cout << "RotateLeft ("<<start<<", "<<n<<" sign: " << sign << endl;
  int tmp = MyPermutation[start];
  for (int i = start; i < n-1; i++) {
    MyPermutation[i] = MyPermutation[i+1];
  }
  MyPermutation[n-1] = tmp;
  sign ^= (((n-start)^1)&1);
} // RotateLeft


void SimplePermutations::Permute(const int start, const int n, int sign)
{
  for (int i=0; i<NbrElements; ++i)
    this->Permutations[CurrentPermutation][i]=MyPermutation[i];
  ++CurrentPermutation;
  //cout << "Permute ("<<start<<", "<<n<<" sign: " << sign << endl;
  if (start < n) {
    int i, j;
    for (i = n-2; i >= start; i--) {
      for (j = i + 1; j < n; j++) {
	this->Swap(i, j, sign);
	Permute(i+1, n, sign);
      } // for j
      RotateLeft(i, n, sign);
    } // for i
  }
} // Permute
