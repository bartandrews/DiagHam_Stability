////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//         Copyright (C) 2007 Gunnar Moeller               //
//                                                                            //
//                                                                            //
//           class for storing the history of a Monte-Carlo overlap calculation  //
//                                                                            //
//                        last modification : 28/05/2007                      //
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

#ifndef ABSTRACTMCHISTORYDATA_H
#define ABSTRACTMCHISTORYDATA_H

class AbstractMCHistoryData
{
 protected:
  unsigned MCHistoryDataType;
  
 public:
  enum
    {
      SomeMCHistoryDataType = 0x01u
    };

  virtual ~AbstractMCHistoryData(){};
  
  // Accessor routine to test for Type of additional data:
  unsigned GetHistoryDataType() {return this->MCHistoryDataType;}

  // Request for written size of additional data
  virtual unsigned GetSizeOnDisk() {return 0;}
  // wrapper routines for writing and reading History Data
  virtual bool WriteMCHistoryData() {return false;}
  virtual bool ReadMCHistoryData() {return false;}
};

#endif // ABSTRACTMCHISTORYDATA_H
