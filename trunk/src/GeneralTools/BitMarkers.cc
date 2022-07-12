#include "BitMarkers.h"

// standard constructor - create empty array
//
BitMarkers::BitMarkers()
{
  this->Length=0;
  this->NbrWords=0;
  this->MarkerArray=NULL;
}


// constructor - create array of given length
// length = length in bits
BitMarkers::BitMarkers(unsigned length)
{
  this->Length=length;
  this->NbrWords=(length+31)/32;
  this->MarkerArray=new unsigned[Length];
  for (unsigned w=0; w<NbrWords; ++w)
    this->MarkerArray[w]=0x0u;
}


// standard destructor
//
BitMarkers::~BitMarkers()
{
  if (this->MarkerArray!=NULL)
    delete[] this->MarkerArray;
}

// copy constructor
// do full copy
BitMarkers::BitMarkers(const BitMarkers &markers)
{
  this->Length = markers.Length;
  this->NbrWords = markers.NbrWords;
  this->MarkerArray=new unsigned[NbrWords];
  for (unsigned w=0; w<NbrWords; ++w)
    this->MarkerArray[w]=markers.MarkerArray[w];
}

// assignment operator
BitMarkers& BitMarkers::operator = (const BitMarkers &markers)
{
  if (this->NbrWords!=markers.NbrWords)
    {
      delete[] this->MarkerArray;
      this->MarkerArray=new unsigned[markers.NbrWords];
      this->NbrWords= markers.NbrWords;
    }
  this->Length = markers.Length;
  for (unsigned w=0; w<NbrWords; ++w)
    this->MarkerArray[w]=markers.MarkerArray[w];
  return *this;
}


