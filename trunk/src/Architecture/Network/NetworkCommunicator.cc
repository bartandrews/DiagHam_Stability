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


#include "config.h"
#include "Architecture/Network/NetworkCommunicator.h"

#ifdef __NETWORK__
#include <unistd.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif

// constructor for 
//
// port = port to be opened by the server 

NetworkCommunicator::NetworkCommunicator(int port)
{
  this->Port = port;
}
  
// destructor
//

NetworkCommunicator::~NetworkCommunicator()
{
}

// try to enable communication with a client
// 
// return value = true if communications are allowed

bool NetworkCommunicator::EnableClientCommunication()
{
#ifdef __NETWORK__
  this->ServerAddress()
  this->ServerAddress.sin_family = AF_INET;
  this->ServerAddress.sin_port = htons(this->Port);
  this->ServerAddress.sin_addr.s_addr = INADDR_ANY;
  bzero(&(this->ServerAddress.sin_zero), 8);
  if (bind(this->ServerSocketFileDescriptor, (sockaddr*) &(this->ServerAddress), sizeof(sockaddr)) == -1)
    {
#ifdef __DEBUG__
      cout << "error bind" << endl;
#endif
      return false;
    }
  if (listen(this->ServerSocketFileDescriptor, this->ServerAddress.sin_port) == -1)
    {
#ifdef __DEBUG__
      cout << "error listen" << endl;
#endif
      return false;
    }
  int SinSize = sizeof (sockaddr_in);
  if ((this->ClientSocketFileDescriptor = accept(this->ServerSocketFileDescriptor, 
						 (sockaddr*) &(this->ClientAddress), (socklen_t*) &SinSize)) == -1)
    {
#ifdef __DEBUG__
      cout << "error listen" << endl;
#endif
      return false;
    }      
#ifdef __DEBUG__
  cout << "connection from" << inet_ntoa(this->ClientAddress.sin_addr) << endl;
#endif
  u_int32 Tag = NetworkCommunicator::IsAlive;
  if ((send(this->ClientSocketFileDescriptor, (char*) &Tag, sizeof(u_int32), 0) == -1))
    {
#ifdef __DEBUG__
      cout << "send error" << endl;
#endif
      return false;
    }
  int BufferSize;
  if (((BufferSize = recv(this->ClientSocketFileDescriptor, (char*) &Tag, 
			  sizeof(u_int32), 0)) == -1) || (BufferSize != sizeof(u_int32)))
    {
#ifdef __DEBUG__
      cout << "receive error" << endl;
#endif
      return false;
    }
  if (Tag == NetworkCommunicator::ImAlive)
    {
      this->ClientFlag = false;
      return true;
    }
  else
    return false;
#else
  return false;
#endif
}

// try to enable communication with server
// 
// server name = fully qualified client name
// return value = true if communications are allowed

bool NetworkCommunicator::EnableServerCommunication(char* serverName)
{
#ifdef __NETWORK__
  hostent* Server;
  if ((Server = gethostbyname(serverName)) == 0)
    {
#ifdef __DEBUG__
      cout << "host " << serverName << " does not exist" << endl;
#endif
      return false;
    }
  if ((SocketFileDescriptor = socket(AF_INET, SOCK_STREAM, 0)) == -1)
    {
#ifdef __DEBUG__
      cout << "error socket" << endl;
#endif
      return false;
    }
  this->ServerAddress.sin_family = AF_INET;
  this->ServerAddress.sin_port = htons(3910);
  this->ServerAddress.sin_addr = *((in_addr*) Server->h_addr);
  bzero(&(this->ServerAddress.sin_zero), 8);
  if (connect(this->ServerSocketFileDescriptor, (sockaddr*) &(this->ServerAddress), sizeof (sockaddr)) == -1)
    {
#ifdef __DEBUG__
      cout << "connect error" << endl;
#endif
      return false;
    }
  return true;
#else
  return false;
#endif
} 

// close communication
//
// return value = true if communications have been closed

bool NetworkCommunicator::CloseCommunication()
{
}

// wait for a new job
//
// return value = id of the newly available job (-1 if an error occurs)
int NetworkCommunicator::WaitJob ();
{
  if (this->ClientFlag == true)
    {
    }
  else
    {
    }
}

