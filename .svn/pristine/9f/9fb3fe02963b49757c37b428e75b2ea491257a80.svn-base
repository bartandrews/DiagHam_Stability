#include "AbstractZDensityProfile.h"
#include "InfiniteWellDensityProfile.h"
#include "FangHowardDensityProfile.h"
#include "TabulatedDensityProfile.h"
#include "InfiniteWellBilayerDensityProfile.h"
#include "InfiniteWellMixedDensityProfile.h"

#include <iostream>

using std::cout;
using std::endl;

AbstractZDensityProfile::~AbstractZDensityProfile()
{
}



// a static class function to return an actual DensityProfile object of some type
AbstractZDensityProfile* AbstractZDensityProfile::CreateZDensityProfile (char *type, double width)
{
  int ProfileType=0;
  if (strlen(type)<3)
    ProfileType=atoi(type);  
  switch(ProfileType)
    {
    case AbstractZDensityProfile::InfiniteWell:
      return new InfiniteWellDensityProfile(width);
      break;
    case AbstractZDensityProfile::FangHoward:
      return (AbstractZDensityProfile*) new FangHowardDensityProfile(width);
      break;
    case AbstractZDensityProfile::TabulatedProfile:
      return (AbstractZDensityProfile*) new TabulatedDensityProfile(type, width);
      break;
    case AbstractZDensityProfile::InfiniteWellExc:
      return new InfiniteWellDensityProfile(width,1);
      break;
    case AbstractZDensityProfile::InfiniteWellBilayerLeft:
      return new InfiniteWellBilayerDensityProfile(width,-1);
      break;
    case AbstractZDensityProfile::InfiniteWellBilayerRight:
      return new InfiniteWellBilayerDensityProfile(width,+1);
      break;
    case AbstractZDensityProfile::InfiniteWellMixed:
      return new InfiniteWellMixedDensityProfile(width);
      break;
    default:
      cout << "This type of Density Profile is not defined, yet"<<endl;
      return 0;
      break;
    } 
}

char *AbstractZDensityProfile::DensityProfileName(char *type)
{
  char * buffer = new char[1000];
  int ProfileType=0;
  if (strlen(type)<3)
    ProfileType=atoi(type);  
  switch(ProfileType)
    {
    case AbstractZDensityProfile::InfiniteWell:
      sprintf(buffer,"Infinite Well Potential");
      break;
    case AbstractZDensityProfile::FangHoward:
      sprintf(buffer,"Fang Howard Potential");
      break;
    case AbstractZDensityProfile::TabulatedProfile:
      sprintf(buffer,"tabulated file: %s",type);
      break;
    case AbstractZDensityProfile::InfiniteWellExc:
      sprintf(buffer,"Infinite Well Potential, 1st excited state");
      break;
    case AbstractZDensityProfile::InfiniteWellBilayerLeft:
      sprintf(buffer,"Infinite Well Potential, 1/sqrt(2)[ |0> - |1> ] state");
      break;
    case AbstractZDensityProfile::InfiniteWellBilayerRight:
      sprintf(buffer,"Infinite Well Potential, 1/sqrt(2)[ |0> + |1> ] state");
      break;
    case AbstractZDensityProfile::InfiniteWellMixed:
      sprintf(buffer,"Infinite Well Mixed Potential");
      break;
    default:
      sprintf(buffer,"Unknown");
      break;
    } 
  char *rst = new char[strlen(buffer)+1];
  strcpy(rst,buffer);
  delete [] buffer;
  return rst;

}
