////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of network communicator                      //
//                                                                            //
//                        last modification : 13/04/2003                      //
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


#ifndef NETWORKCOMMUNICATOR_H
#define NETWORKCOMMUNICATOR_H


#include "config.h"


class NetworkCommunication
{

 private:

  // client flag 
  bool ClientFlag;

  // tcp/ip port 
  int Port;

  // distant name (if any)
  char* DistantName;
  // local name (if any)
  char* LocalName;

#ifdef __NETWORK__
  // socket file descriptor for the server
  int ServerSocketFileDescriptor;
  // socket file descriptor for the client
  int ClientSocketFileDescriptor;

  // structure storing details of the server 
  sockaddr_in ServerAddress;
  // structure storing details of the client 
  sockaddr_in ClientAddress;

#endif  

  // job id counter
  int JobID;
  
  // reserved code describing each communication chunk (ranging from 1 to 255)
  enum CommunicationCode
    {
      // ask if the distant client is alive
      IsAlive = 0x00000001,
      // return that local server is alive
      ImAlive = 0x00000002,
      // close communication
      Close = 0x00000003
      // start a new job
      NewJob = 0x00000004
      // end a job
      EndJob = 0x00000005
    };

 public:
  
  // constructor for 
  //
  // port = port to be opened by the server 
  NetworkCommunication(int port);
  
  // destructor
  //
  ~NetworkCommunication();
  
  // try to enable communication with a client
  // 
  // client name = fully qualified client name
  // return value = true if communications are allowed
  bool EnableClientCommunication(char* clientName);

  // try to enable communication with server
  // 
  // server name = fully qualified client name
  // return value = true if communications are allowed
  bool EnableServerCommunication(char* serverName);
 
  // close communication
  //
  // return value = true if communications have been closed
  bool CloseCommunication();

  // wait for a new job
  //
  // return value = id of the newly available job (-1 if an error occurs)
  int WaitJob ();

};

#endif
