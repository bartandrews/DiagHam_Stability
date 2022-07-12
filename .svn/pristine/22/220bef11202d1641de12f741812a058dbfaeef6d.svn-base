#include "config.h"
#include "GenericSignalHandler.h"
#include <signal.h>

#include <iostream>
using std::cout;
using std::endl;

// Flag indicating whether error handler was initialized
volatile bool GenericSignalHandler::SignalHandlerInitialized = false;

// current number of signal handlers
volatile sig_atomic_t GenericSignalHandler::NbrSignalHandlers=0;

GenericSignalHandler** GenericSignalHandler::DiagHamSignalHandlers=NULL;


// standard constructor
// signum = number of signal to be caught
// iterateOnRelease = when signal is released, execute previous default signal handler
GenericSignalHandler::GenericSignalHandler(int signum, bool iterateOnRelease)
{
  ++NbrSignalHandlers;
  this->SignalNumber=signum;
  this->IteratePreviousSignalOnRelease=iterateOnRelease;
  this->SignalHandlerSet=false;
  this->SignalDeferred=0;
  this->SignalPending=0;
  this->SignalHandlerSet = false;
  if (!SignalHandlerInitialized)
    {
      // initialize signal handler fields
      DiagHamSignalHandlers = new GenericSignalHandler*[NBRDIAGHAMSIGNALHANDLERS];
      for (int i=0; i<NBRDIAGHAMSIGNALHANDLERS; ++i)
	DiagHamSignalHandlers[i]=NULL;
      SignalHandlerInitialized=true;
    }
  else if (DiagHamSignalHandlers[signum]!=NULL)
    {
      cout << "Currently, only a single signal handler can be defined for each signal"<<endl;
      exit(-1);
    }
  DiagHamSignalHandlers[signum]=this;
}

// destructor
GenericSignalHandler::~GenericSignalHandler()
{
  --NbrSignalHandlers;
  signal(this->SignalNumber,SIG_DFL);
  DiagHamSignalHandlers[this->SignalNumber]=NULL;
}


// function that will be executed when a signal is caught (and deferred)
void GenericSignalHandler::FunctionOnSignal()
{
  cout << "Deferring signal no. "<<SignalNumber<<" in GenericSignalHandler"<<endl;
}
  
// function that will be executed when the signal is released
void GenericSignalHandler::FunctionOnRelease()
{
  cout << "Executing function for deferred signal no. "<<SignalNumber<<endl;
  if (IteratePreviousSignalOnRelease)
    cout << "The standard signal will be executed in sequel."<<endl;
}

  
// function called when a signal is caught
void GenericSignalHandler::CatchSignal(int signum)
{
  this->SignalPending = signum;
  this->FunctionOnSignal();
}

// routine to call to block Signal
void GenericSignalHandler::StartToDeferSignal()
{
  //  cout << "Activating deferred signal no. "<<SignalNumber<<endl;
  if (this->SignalHandlerSet == false)
    {
      signal(this->SignalNumber, &DiagHamSignalHandler);
      this->SignalHandlerSet=true;
    }
  SignalDeferred++;
}

// routine to call to process pending Signal using build-in function
void GenericSignalHandler::ProcessDeferredSignal()
{
  //cout << "Processing signals status: "<<SignalDeferred<<", pending: "<<SignalPending<<endl;
  SignalDeferred--;
  if ((SignalDeferred == 0) && (SignalPending != 0))
    {
      this->FunctionOnRelease();
      if (IteratePreviousSignalOnRelease)
	raise (SignalPending);
    }
  SignalHandlerSet = false;
  signal(this->SignalNumber, SIG_DFL); // register default SIGTERM handler
}


// routine to call to clear pending Signal
// revertHandler = flag indicating whether the standard-handler shall be put in place
void GenericSignalHandler::ClearDeferredSignal(bool revertHandler)
{
  //  cout << "Clearing signals status: "<<SignalDeferred<<", pending: "<<SignalPending<<endl;
  SignalPending=0;
  if (revertHandler)
    {
      SignalDeferred--;
      SignalHandlerSet = false;
      signal(this->SignalNumber, SIG_DFL); // register default SIGTERM handler
    }
}


/* ===== GLOBAL FUNCTIONS  =====  */

// signal handling routines for deferring exit in  with disk storage
// actual signal handler
void DiagHamSignalHandler(int signum)
{
  if ((GenericSignalHandler::DiagHamSignalHandlers[signum]!=NULL)&&(GenericSignalHandler::DiagHamSignalHandlers[signum]->IsDeferred()))
    {
      GenericSignalHandler::DiagHamSignalHandlers[signum]->CatchSignal(signum);
    }
  else
    {
      // reiterate original system-wide term-signal here
      signal(signum, SIG_DFL); // register default SIGTERM handler
      raise(signum);
    }
}
