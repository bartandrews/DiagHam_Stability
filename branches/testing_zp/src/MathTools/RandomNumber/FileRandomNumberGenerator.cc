////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of random number generator from a external generated file       //
//                                                                            //
//                        last modification : 08/04/2008                      //
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
#include "MathTools/RandomNumber/FileRandomNumberGenerator.h"
#include "GeneralTools/Endian.h"

#include <iostream>


using std::ios;
using std::cout;
using std::endl;


// constructor
//
// inputFile = name of the file where random numbers are stored
// nbrRandomNumbers = number of random numbers to read from inputFile
// seekPosition = position of random numbers to read

FileRandomNumberGenerator::FileRandomNumberGenerator(char* inputFile, unsigned long nbrRandomNumbers, unsigned long seekPosition)
{
  this->NbrRandomNumbers = nbrRandomNumbers;
  this->RealRandomNumbers = new double [this->NbrRandomNumbers];
  this->IntegerRandomNumbers = new unsigned long [this->NbrRandomNumbers];  
  this->Flag.Initialize();
  this->Position = 0ul;
  this->iset=0;
  ifstream File;
  File.open(inputFile, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << inputFile << endl;
      this->NbrRandomNumbers = 0;
      return;
    }  
  unsigned long TmpNbrRandomNumbers = 0;
  ReadLittleEndian(File, TmpNbrRandomNumbers);
  ReadLittleEndian(File, this->MaxInteger);
  if ((seekPosition + this->NbrRandomNumbers) >  TmpNbrRandomNumbers)
    {
      this->NbrRandomNumbers = TmpNbrRandomNumbers - seekPosition;
    }
  File.seekg((sizeof(unsigned long) + sizeof(double)) * seekPosition, ios::cur);
  cout << this->NbrRandomNumbers << endl;
  for (unsigned long i = 0; i < this->NbrRandomNumbers; ++i)
     {
       ReadLittleEndian(File, this->IntegerRandomNumbers[i]);
       ReadLittleEndian(File, this->RealRandomNumbers[i]);
     }
   File.close(); 
}

// constructor from another random generator, used to generate a file of random numbers
//
// generator = generator to use to generate the random numbers
// nbrRandomNumbers = number of random numbers to generate and to write into outputFile
// outputFile = name of the file where random numbers will be stored

FileRandomNumberGenerator::FileRandomNumberGenerator(AbstractRandomNumberGenerator* generator, int nbrRandomNumbers, char* outputFile)
{
  this->NbrRandomNumbers = nbrRandomNumbers;
  this->RealRandomNumbers = new double [this->NbrRandomNumbers];
  this->IntegerRandomNumbers = new unsigned long [this->NbrRandomNumbers];
  this->MaxInteger = generator->GetMaxInteger();
  double InvMax = 1.0 / ((double) generator->GetMaxInteger());
  this->iset=0;
  for (unsigned long i = 0; i < this->NbrRandomNumbers; ++i)
    {
      this->IntegerRandomNumbers[i] = generator->GetIntegerRandomNumber();
      this->RealRandomNumbers[i] = ((double) this->IntegerRandomNumbers[i]) * InvMax;      
    }  
  this->Position = 0;
  this->Flag.Initialize();
  ofstream File;
  File.open(outputFile, ios::binary | ios::out);
  WriteLittleEndian(File, this->NbrRandomNumbers);
  WriteLittleEndian(File, this->MaxInteger);
  for (unsigned long i = 0; i < this->NbrRandomNumbers; ++i)
    {
      WriteLittleEndian(File, this->IntegerRandomNumbers[i]);
      WriteLittleEndian(File, this->RealRandomNumbers[i]);
    }
  File.close();
}

// copy constructor
//
// generator = generator to copy

FileRandomNumberGenerator::FileRandomNumberGenerator(const FileRandomNumberGenerator& generator)
{
  this->Flag = generator.Flag;
  this->RealRandomNumbers = generator.RealRandomNumbers;
  this->IntegerRandomNumbers = generator.IntegerRandomNumbers;
  this->NbrRandomNumbers = generator.NbrRandomNumbers;
  this->Position = generator.Position;
  this->MaxInteger = generator.MaxInteger;
}

// destructor
//

FileRandomNumberGenerator::~FileRandomNumberGenerator()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete[] this->RealRandomNumbers;
      delete[] this->IntegerRandomNumbers;
    }
}

// clone random number generator 
//
// return value = clone of the random number generator

AbstractRandomNumberGenerator* FileRandomNumberGenerator::Clone ()
{
  return new FileRandomNumberGenerator(*this);
}

