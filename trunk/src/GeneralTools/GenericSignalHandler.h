////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of abstract Lanczos algorithm                   //
//                                                                            //
//                        last modification : 30/04/2001                      //
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


#ifndef GENERICSIGNALHANDLER_H
#define GENERICSIGNALHANDLER_H


#include "config.h"

#include <cstdlib>
#include <signal.h>


#define NBRDIAGHAMSIGNALHANDLERS 34


class GenericSignalHandler
{

 protected:
  // code of signal that is being handled
  int SignalNumber;

  // flag indicating whether signal shall be passed on, or deleted upon deletion of this objectt
  volatile bool IteratePreviousSignalOnRelease;

  // function that will be executed when a signal is caught (and deferred)
  virtual void FunctionOnSignal();

  // function that will be executed when the signal is released
  virtual void FunctionOnRelease();
  
  /* If this flag is nonzero, don't handle the signal right away. */
  volatile sig_atomic_t SignalPending;
  
  /* This is nonzero if a signal arrived and was not handled. */
  volatile sig_atomic_t SignalDeferred;

  /* Flag indicating whether error handler was set */
  volatile bool SignalHandlerSet;

  /* Flag indicating whether global array of error handlers was initialized */
  static volatile bool SignalHandlerInitialized;

  // total number of signal handlers
  static volatile sig_atomic_t NbrSignalHandlers;

  
 public:

  // standard constructor
  // signum = number of signal to be caught
  // iterateOnRelease = when signal is released, execute previous default signal handler
  GenericSignalHandler(int signum, bool iterateOnRelease=true);

  // destructor
  ~GenericSignalHandler();

  // start to defer the signal
  void StartToDeferSignal();

  // process any pending deferred signals
  void ProcessDeferredSignal();

  // routine to call to clear pending Signal
  // revertHandler = flag indicating whether the standard-handler shall be put in place
  //
  void ClearDeferredSignal(bool revertHandler=true);

  // accessor routine to request whether signal is pending
  bool HavePendingSignal(){return (SignalPending!=0);}

  // function called when a signal is caught
  void CatchSignal(int signum);

  // accessor routines for deferred signal flag
  sig_atomic_t IsDeferred(){return this->SignalDeferred;}
  //sig_atomic_t& SignalPending(){return this->SignalDeferred;}

  static GenericSignalHandler** DiagHamSignalHandlers;

};



/* ===== GLOBAL FUNCTIONS  =====  */

// signal handling routines for deferring exit in  with disk storage
// actual signal handler
void DiagHamSignalHandler(int signum);


#endif  // SIGNALMANAGER_H


/* =========================================================================== */
/* for convenience, the following includes a list of unix signals
SIGHUP	1	Exit	Hangup
SIGINT	2	Exit	Interrupt
SIGQUIT	3	Core	Quit
SIGILL	4	Core	Illegal Instruction
SIGTRAP	5	Core	Trace/Breakpoint Trap
SIGABRT	6	Core	Abort
SIGEMT	7	Core	Emulation Trap
SIGFPE	8	Core	Arithmetic Exception
SIGKILL	9	Exit	Killed
SIGBUS	10	Core	Bus Error
SIGSEGV	11	Core	Segmentation Fault
SIGSYS	12	Core	Bad System Call
SIGPIPE	13	Exit	Broken Pipe
SIGALRM	14	Exit	Alarm Clock
SIGTERM	15	Exit	Terminated
SIGUSR1	16	Exit	User Signal 1
SIGUSR2	17	Exit	User Signal 2
SIGCHLD	18	Ignore	Child Status
SIGPWR	19	Ignore	Power Fail/Restart
SIGWINCH	20	Ignore	Window Size Change
SIGURG	21	Ignore	Urgent Socket Condition
SIGPOLL	22	Ignore	Socket I/O Possible
SIGSTOP	23	Stop	Stopped (signal)
SIGTSTP	24	Stop	Stopped (user)
SIGCONT	25	Ignore	Continued
SIGTTIN	26	Stop	Stopped (tty input)
SIGTTOU	27	Stop	Stopped (tty output)
SIGVTALRM	28	Exit	Virtual Timer Expired
SIGPROF	29	Exit	Profiling Timer Expired
SIGXCPU	30	Core	CPU time limit exceeded
SIGXFSZ	31	Core	File size limit exceeded
SIGWAITING	32	Ignore	All LWPs blocked
SIGLWP	33	Ignore	Virtual Interprocessor Interrupt for Threads Library
SIGAIO	34	Ignore	Asynchronous I/O
   =========================================================================== */
