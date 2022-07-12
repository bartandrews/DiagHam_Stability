#include "PermutationElement.h"

PermutationElement::PermutationElement()
{
  this->Multiplicity=1;
}

PermutationElement::PermutationElement(SmallIntegerArray value, unsigned long multiplicity)
{
  this->Value = value;
  this->Multiplicity = multiplicity;
}

PermutationElement::PermutationElement(const PermutationElement &per)
{
  this->Value = per.Value;
  this->Multiplicity = per.Multiplicity;
}

PermutationElement::~PermutationElement()
{
}


// assignment operator
//
// per = permutation to assign
// return value = reference on current permutation
PermutationElement& PermutationElement::operator = (const PermutationElement& per)
{
  this->Value = per.Value;
  this->Multiplicity = per.Multiplicity;
  return *this;
}

// multiply on multiplicity
PermutationElement& PermutationElement::operator *= (const unsigned long mult)
{
  this->Multiplicity*=mult;
  return *this;
}


bool operator == (const PermutationElement& a1, const PermutationElement& a2)
{
  return (a1.Value==a2.Value);
}

bool operator < (const PermutationElement& a1,const PermutationElement& a2)
{
  return (a1.Value<a2.Value);
}

bool operator > (const PermutationElement& a1,const PermutationElement& a2)
{
  return (a1.Value>a2.Value);
}

ostream& operator << ( ostream &Str, PermutationElement PE)
{
  Str << "[" << PE.Value <<"] x "<<PE.Multiplicity;
  return Str;
}
