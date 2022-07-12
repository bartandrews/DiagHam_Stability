////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            DarkRay version  0.05                           //
//                                                                            //
//                  Copyright (C) 1998-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of option manager                         //
//                                                                            //
//                        last modification : 19/05/2004                      //
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
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/SystemTools.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/OptionManager.h"
#include "Options/Options.h"

#include <iostream>
#include <string.h>
#include <cstdlib>
#include <sys/time.h>


using std::ostream;
using std::endl;
using std::cout;

#define QUOTEME(x) #x


// constructor
//
// groupName = group full name

OptionManager::OptionManager(const char* programName, const char* programVersion, const char* programAdditionalInformations)
{
  this->ProgramName = new char [strlen(programName) + 1];
  strcpy (this->ProgramName, programName);      
  if (programVersion == 0)
    {
      this->ProgramVersion = 0;
    }
  else
    {
      this->ProgramVersion = new char [strlen(programVersion) + 1];
      strcpy (this->ProgramVersion, programVersion);      
    }
  if (programAdditionalInformations == 0)
    {
      this->ProgramAdditionalInformations = 0;
    }
  else
    {
      this->ProgramAdditionalInformations = new char [strlen(programAdditionalInformations) + 1];
      strcpy (this->ProgramAdditionalInformations, programAdditionalInformations);      
    }

  this->GlobalGroup  = new OptionGroup ("global options");
  /*
  (*GlobalGroup) += new SingleStringOption  ('\n', "save-command", "save detailed settings of all options to file");
  (*GlobalGroup) += new SingleStringOption  ('\n', "load-command", "load detailed option parameters from file");
  */
#ifdef HAVE_GLOBAL_COMMAND_LOG
  char Buffer[1024]=GLOBAL_COMMAND_LOG;
  char *TmpC=new char[1024];
  sprintf (TmpC,"append command line to file [overrides default '%s']", Buffer );
  (*GlobalGroup) += new SingleStringOption  ('\n', "append-cmdline", TmpC);
  delete [] TmpC;
#else
  (*GlobalGroup) += new SingleStringOption  ('\n', "append-cmdline", "append command line to file");
#endif
  (*GlobalGroup) += new SingleStringInternalOption  ('\n', "repeat-cmdline", "load command line from file [format: filename.line]");
#ifdef HAVE_GLOBAL_COMMAND_LOG
  (*GlobalGroup) += new BooleanOption ('\n',"cmdlog-off","Do not log command line");
#else
  (*GlobalGroup) += new BooleanInternalOption ('\n',"cmdlog-off","Do not log command line");
#endif
}

// destructor
//

OptionManager::~OptionManager()
{
  delete[] this->ProgramName;
  if (this->ProgramVersion != 0)
    delete[] this->ProgramVersion;
  if (this->ProgramAdditionalInformations != 0)
    delete[] this->ProgramAdditionalInformations;
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  while ((TmpGroup = IterGroup()))
    {
      delete *TmpGroup;
    }
  delete GlobalGroup;
}

// add an option group to the manager
// 
// group = pointer to the option group to add 
// return value = reference on the current option manager

OptionManager& OptionManager::operator += (OptionGroup* group)
{
  this->Groups += group;
  return *this;
}

// get option from its name
//
// optionName = string containing option name
// return value = poitner to the option if it has been found, 0 either

AbstractOption* OptionManager::operator[] (const char* optionName)
{
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  AbstractOption* TmpOption = 0;
  while ((TmpOption == 0) && ((TmpGroup = IterGroup())))
    {
      TmpOption = (**TmpGroup)[optionName];
    }
  if (TmpOption==0)
    TmpOption = (*GlobalGroup)[optionName];
  return TmpOption; 
}

// get an option group from its name
//
// optionGroupName = string containing option group name
// return value = pointer to the option group if it has been found, 0 either

OptionGroup* OptionManager::GetOptionGroup(const char* optionGroupName){
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  while ((TmpGroup = IterGroup()))
    {
      if ((*TmpGroup)->IsGroupName(optionGroupName) == true)
	return (*TmpGroup);
    }
  if (GlobalGroup->IsGroupName(optionGroupName) == true)
    return (GlobalGroup);
  return 0;  
}

// Proceed running options from command line arguments
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// output = reference on output stream used to display errors
// return value = true if proceeding succeded, false if an error occurs

bool OptionManager::ProceedOptions (char** argumentValues, int nbrArgument, ostream& output)
{
  int Pos = 1;
  int Inc;
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;

  // check for repeating of stored command line
  AbstractOption* RepeatOption = (*this)["repeat-cmdline"];
  while (Pos < nbrArgument)
    {
      Inc = RepeatOption->ReadOption(argumentValues, nbrArgument, Pos);
      if (Inc == -1)
	{
	  RepeatOption->PrintError(output);
	  return false; 
	}
      if (Inc > 0)
	break;
      ++Pos;
    }
  if (this->GetString("repeat-cmdline"))
    {
      char **Arguments;
      int NbrArguments = SplitLine(this->GetString("repeat-cmdline"), Arguments, ',');
      int NbrLine=1;
      if (NbrArguments==2)
	{
	  NbrLine=strtod(Arguments[1],NULL);
	}
      char *CommandFile = Arguments[0];
      char *Line =  GetLineFromFile (CommandFile, NbrLine);
      for (int i=0; i<NbrArguments; ++i)
	delete [] Arguments[i];
      delete [] Arguments;
      char **Columns;
      int NbrColumns = SplitLine(Line, Columns, '\t');
      if (NbrColumns != 6)
	{
	  cout << "Problem with number of columns: command line not identified"<<endl;
	  exit(1);
	}
      NbrArguments = SplitLine(Columns[6], Arguments, ' ');
      // test same executable
      char *Path,*StoredFileName;
      ExtractPathAndFileName (Arguments[0], Path, StoredFileName);
      if (Path!=NULL) delete [] Path;
      char *CurrentFileName;
      ExtractPathAndFileName (argumentValues[0], Path, CurrentFileName);
      if (Path!=NULL) delete [] Path;
      if (strcmp(StoredFileName,CurrentFileName)!=0)
	{
	  cout << "Could not read command, as executables do not coincide. You can execute the command manually:"<<endl;
	  cout << Columns[6];
	  exit(-1);
	}
      cout << "Repeating the following command:"<<endl;
      cout << Columns[6];
      return this->ProceedOptions (Arguments, NbrArguments, output);
    }
  // evaluate other options
  Pos=1;
  while (Pos < nbrArgument)
    {
      Inc = 0;
      IterGroup.DefineList(this->Groups);
      while ((Inc == 0) && ((TmpGroup = IterGroup())) )
	{
	  Inc = (*TmpGroup)->ReadOption(argumentValues, nbrArgument, Pos);
	  if (Inc == -1)
	    {
	      (*TmpGroup)->PrintError(output);
	      return false; 
	    }
	}
      if (Inc == 0)
	{
	  Inc = this->GlobalGroup->ReadOption(argumentValues, nbrArgument, Pos);
	  if (Inc == -1)
	    {
	      (*TmpGroup)->PrintError(output);
	      return false; 
	    }
	  if (Inc == 0)
	    {
	      output << "unknown option " <<  argumentValues[Pos] << endl;
	      return false;
	    }
	}
      Pos += Inc;
    }
  // proceed global options:

  /*
  if (this->GetString("save-command")) // save detailed command line options to file
    {
      cout << "option '--save-command' not implemented, yet"<<endl;
      exit(1);
    }
  if (this->GetString("load-command")) // load detailed list of command line options from file
    {
      cout << "option '--load-command' not implemented, yet"<<endl;
      exit(1);
    }
  */
#ifdef HAVE_GLOBAL_COMMAND_LOG
  if (this->GetBoolean("cmdlog-off")==false)
    {
      char *AppendFile;
      if (this->GetString("append-cmdline"))
	{
	  AppendFile = new char[strlen(this->GetString("append-cmdline"))+1];
	  strcpy (AppendFile,this->GetString("append-cmdline"));
	}
      else
	{
	  char Buffer[1024]=GLOBAL_COMMAND_LOG;
	  AppendFile = new char[1024];
	  strcpy (AppendFile,Buffer);
	}
#else
  if (this->GetString("append-cmdline"))
    {
      char *AppendFile = new char[strlen(this->GetString("append-cmdline"))+1];
      strcpy (AppendFile,this->GetString("append-cmdline"));
#endif
      ofstream File;
      ifstream TestFile;
      TestFile.open(AppendFile, ios::in);
      long Count;
      if (TestFile.is_open())
	{
	  TestFile.close();
	  Count=GetFileNbrLines(AppendFile);
	  File.open(AppendFile, ios::app);
	}
      else
	{
	  File.open(AppendFile, ios::out);
	  File << "#item\ttime\thost\tcmd-line"<<endl;
	  Count=1;
	}
      timeval LaunchTime;
      gettimeofday (&(LaunchTime), 0);
      
      //Find the current time as a string
      time_t CurrentTime = time(0);
      //convert it to tm
      tm Now=*localtime(&CurrentTime);
      char TimeStr[BUFSIZ]={0};
      //Format string determines the conversion specification's behaviour
      const char Format[]="%m/%d/%y-%X";
      
      //strftime - converts date and time to a string
      if (strftime(TimeStr, sizeof(TimeStr)-1, Format, &Now)<=0)
	{
	  std::cout<<"Could not convert time";
	  exit(1);
	}

      // alternative:
      // #include <unistd.h>
      // char *getcwd(char *buf, size_t size); 
      
      char *CmdString = new char[1024];
      // get full executable path and machine name
      sprintf(CmdString,"which %s",argumentValues[0]);
      char *Executable = GetLineOutputFromSystemCommand(CmdString);
      
      sprintf(CmdString,"hostname");
      char *Host = GetLineOutputFromSystemCommand(CmdString);

      sprintf(CmdString,"pwd");
      char *Pwd = GetLineOutputFromSystemCommand(CmdString);

      //LaunchTime.tv_sec
      File << Count << "\t" << TimeStr << "\t" << Host << "\t" << Pwd << "\t" << Executable << "\t";
      for (Pos=0; Pos<nbrArgument-1; ++Pos)
	File << argumentValues[Pos] << " ";
      File << argumentValues[Pos] << endl;
      delete [] Executable;
      delete [] CmdString;
      delete [] Pwd;
      delete [] Host;
      delete [] AppendFile;
    }
  return true;
}

// ProceedOptions and test if help should be displayed
//
// argumentValues = string array of arguments
// nbrArgument = number of arguments in argumentValues array
// output = reference on output stream used to display errors  
void OptionManager::StandardProceedings(char** argumentValues, int nbrArgument, ostream& output)
{
  if (this->ProceedOptions(argumentValues, nbrArgument, output) == false)
    {
      cout << "see man page for option syntax or type " << this->ProgramName << " -h" << endl;
      exit(-1);
    }
  
  if (this->GetBoolean("help") == true)
    {
      this->DisplayHelp (output);
      exit(0);
    }
}

// print the options and their values in the current group
//  
// output = reference on output stream;
// shortVersion = true if return only option code and the option value, false if return option description in addition
// comment = if different from the null character, add it in front of each line
// return value = reference on current output stream

ostream& OptionManager::DisplayOption (ostream& output, bool shortVersion, char comment)
{
  if (shortVersion)
    {
      if (comment != '\0')
	output << comment << " ";
      output << this->ProgramName << "  ";
      ListIterator<OptionGroup*> IterGroup(this->Groups);
      OptionGroup** TmpGroup;
      while ((TmpGroup = IterGroup()))
	{
	  (*TmpGroup)->DisplayOption(output, shortVersion, comment);
	}
      return output;
    }
  else
    {
      if ((this->ProgramVersion != 0) || (this->ProgramAdditionalInformations != 0))
	{
	  if (comment != '\0')
	    output << comment << " ";
	  output << this->ProgramName;
	  if (this->ProgramVersion != 0)
	    {
	      output << ", version " << this->ProgramVersion << endl;
	    }
	  else
	    {
	      output << endl;
	    }
	  if (this->ProgramAdditionalInformations != 0)
	    {
	      if (comment != '\0')
		output << comment << " ";
	      output << this->ProgramAdditionalInformations << endl;
	    }
	}  
      if (comment != '\0')
	output << comment << " ";
      output << endl;
      if (comment != '\0')
	output << comment << " ";
      output << "Options:" << endl;
      if (comment != '\0')
	output << comment << " ";
      cout << endl;
      ListIterator<OptionGroup*> IterGroup(this->Groups);
      OptionGroup** TmpGroup;
      while ((TmpGroup = IterGroup()))
	{
	  (*TmpGroup)->DisplayOption(output, shortVersion, comment);
	  if (comment != '\0')
	    output << comment << " ";	  
	  cout << endl;
	}
      return output;
    }
}

// print help concerning current option group
//
// output = reference on output stream;
// return value = reference on current output stream

ostream& OptionManager::DisplayHelp (ostream& output)
{
  if ((this->ProgramVersion != 0) || (this->ProgramAdditionalInformations != 0))
    {
      output << this->ProgramName;
      if (this->ProgramVersion != 0)
	{
	   output << ", version " << this->ProgramVersion << endl;
	}
      else
	{
	  output << endl;
	}
      if (this->ProgramAdditionalInformations != 0)
	{
	  output << this->ProgramAdditionalInformations << endl;
	}
      output << endl;      
    }  
  output << "Usage: " << this->ProgramName << " [options] " << endl;
  output << endl << "Options:" << endl << endl;
  ListIterator<OptionGroup*> IterGroup(this->Groups);
  OptionGroup** TmpGroup;
  while ((TmpGroup = IterGroup()))
    {
      (*TmpGroup)->DisplayHelp(output) << endl;
    }
  GlobalGroup->DisplayHelp(output) << endl;
  return output;
}

// dump some of the options into a formatted string
//
// format = string describing the format to use (each symbol %optionname% is replaced by the value associated to the option referred as optionname)
// return value = formatted string (deallocation has to be done manually, 0 if an error occured)

char* OptionManager::GetFormattedString (const char* format)
{
  char* TmpScratch = new char [strlen(format)];
  int StringLength = 0;
  List<char*> Values;
  char* TmpFormat = (char*)format;
  char** TmpValue;
  AbstractOption* TmpOption;
  while ((*TmpFormat) != '\0')
    {
      while (((*TmpFormat) != '%') && ((*TmpFormat) != '\0'))
	{
	  ++TmpFormat;
	  ++StringLength;
	}
      if ((*TmpFormat) != '\0')
	{
	  ++TmpFormat;
	  int TmpSize = 0;
	  while (((*TmpFormat) != '%') && ((*TmpFormat) != '\0'))
	    {
	      TmpScratch[TmpSize] = (*TmpFormat);
	      ++TmpSize;
	      ++TmpFormat;
	    }
	  if ((*TmpFormat) == '\0')
	    {
	      ListIterator<char*> IterValues(Values);
	      while ((TmpValue = IterValues()))
		{
		  delete[] (*TmpValue);
		} 
	      delete[] TmpScratch;
	      return 0l;
	    }
	  TmpScratch[TmpSize] = '\0';
	  TmpOption = (*this)[TmpScratch];
	  if (TmpOption == 0)
	    {
	      ListIterator<char*> IterValues(Values);
	      while ((TmpValue = IterValues()))
		{
		  delete[] (*TmpValue);
		} 
	      delete[] TmpScratch;
	      return 0l;
	    }
	  else
	    {
	      char* TmpString = TmpOption->GetAsAString();
	      StringLength += strlen (TmpString);
	      Values += TmpString;
	    }
	  ++TmpFormat;
	}
    }
  TmpFormat = (char*)format;
  ListIterator<char*> IterValues(Values);
  char* TmpFormattedString = new char [StringLength + 1];
  char* TmpFormattedString2 = TmpFormattedString;
  while ((*TmpFormat) != '\0')
    {
      while (((*TmpFormat) != '%') && ((*TmpFormat) != '\0'))
	{
	  (*TmpFormattedString2) = (*TmpFormat);
	  ++TmpFormattedString2;
	  ++TmpFormat;
	}
      if ((*TmpFormat) != '\0')
	{
	  ++TmpFormat;
	  while ((*TmpFormat) != '%')
	    {
	      ++TmpFormat;	  
	    }
	  ++TmpFormat;
	  TmpValue = IterValues();
	  strcpy (TmpFormattedString2, *TmpValue);
	  TmpFormattedString2 += strlen(*TmpValue);
	  delete[] (*TmpValue);
	}
    }
  delete[] TmpScratch;
  (*TmpFormattedString2) = '\0';
  return TmpFormattedString;
}


//accessor routine for Boolean value
bool OptionManager::GetBoolean(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTBoolean)
	return ((BooleanOption*)OptionPointer)->GetBoolean();
      else
	{
	  cout << "Boolean value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//accessor routine for Double value
double OptionManager::GetDouble(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTDouble)
	return ((SingleDoubleOption*)OptionPointer)->GetDouble();
      else
	{
	  cout << "Double value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//accessor routine for Multiple Double value
double* OptionManager::GetDoubles(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTDoubles)
	return ((MultipleDoubleOption*)OptionPointer)->GetDoubles();
      else
	{
	  cout << "MultipleDouble value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

// alternative accessor routine for Multiple Double value
double* OptionManager::GetDoubles(const char *optionName, int & length)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTDoubles)
	{
	  length = ((MultipleDoubleOption*)OptionPointer)->GetLength();
	  return ((MultipleDoubleOption*)OptionPointer)->GetDoubles();
	}
      else
	{
	  cout << "MultipleDouble value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//test routine for Multiple Double value
bool OptionManager::HasDoubles(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTDoubles)
	return (((MultipleDoubleOption*)OptionPointer)->GetLength()>0);
      else
	{
	  cout << "MultipleDouble value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}


//accessor routine for Integer value
long OptionManager::GetInteger(const char *optionName)
  {
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTInteger)
	return ((SingleIntegerOption*)OptionPointer)->GetInteger();
      else
	{
	  cout << "Integer value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//accessor routine for Multiple Integer value
int* OptionManager::GetIntegers(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTIntegers)
	return ((MultipleIntegerOption*)OptionPointer)->GetIntegers();
      else
	{
	  cout << "MultipleInteger value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

// alternative accessor routine for Multiple Integer value
int* OptionManager::GetIntegers(const char *optionName, int & length)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTIntegers)
	{
	  length = ((MultipleIntegerOption*)OptionPointer)->GetLength();
	  return ((MultipleIntegerOption*)OptionPointer)->GetIntegers();
	}
      else
	{
	  cout << "MultipleInteger value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//test routine for Multiple Double value
bool OptionManager::HasIntegers(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTIntegers)
	return (((MultipleIntegerOption*)OptionPointer)->GetLength()>0);
      else
	{
	  cout << "MultipleInteger value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}


//accessor routine for String value
char* OptionManager::GetString(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTString)
	return ((SingleStringOption*)OptionPointer)->GetString();
      else
	{
	  cout << "String value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//accessor routine for Multiple String value
char** OptionManager::GetStrings(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTStrings)
	return ((MultipleStringOption*)OptionPointer)->GetStrings();
      else
	{
	  cout << "MultipleString value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

// alternative accessor routine for Multiple String value
char** OptionManager::GetStrings(const char *optionName, int & length)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTStrings)
	{
	  length = ((MultipleStringOption*)OptionPointer)->GetNbrStrings();
	  return ((MultipleStringOption*)OptionPointer)->GetStrings();
	}
      else
	{
	  cout << "MultipleString value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}

//test routine for Multiple String value
bool OptionManager::HasStrings(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTStrings)
	return (((MultipleStringOption*)OptionPointer)->GetNbrStrings()>0);
      else
	{
	  cout << "MultipleString value of Option '"<<optionName<<"' was requested, but is of different type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}


//test routine for Multiple String value
bool OptionManager::HasValues(const char *optionName)
{
  AbstractOption* OptionPointer = (*this)[optionName];
  if (OptionPointer!=0)
    {
      if (OptionPointer->GetOptionType() == AbstractOption::OTStrings)
	return (((MultipleStringOption*)OptionPointer)->GetNbrStrings()>0);
      else if (OptionPointer->GetOptionType() == AbstractOption::OTIntegers)
	return (((MultipleIntegerOption*)OptionPointer)->GetLength()>0);
      else if (OptionPointer->GetOptionType() == AbstractOption::OTDoubles)
	return (((MultipleDoubleOption*)OptionPointer)->GetLength()>0);
      else
	{
	  cout << "HasValues() function of Option '"<<optionName<<"' was requested, but is not of multiple value type!" << endl;
	  exit(-1);
	}
    }
  else
    {
      cout << "Option '"<<optionName<<"' was requested, but is not implemented!" << endl;
      exit(-1);
    }
}
