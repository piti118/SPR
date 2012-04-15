//------------------------------------------------------------------------
// File and Version Information:
//      $Id: SprRWFactory.hh,v 1.1 2007-10-05 20:03:09 narsky Exp $
//
// Description:
//      Class SprRWFactory :
//         Makes writers and readers of the requested type.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2007              California Institute of Technology
//
//------------------------------------------------------------------------
 
#ifndef _SprRWFactory_HH
#define _SprRWFactory_HH

#include "StatPatternRecognition/SprAbsWriter.hh"
#include "StatPatternRecognition/SprAbsReader.hh"

#include "StatPatternRecognition/SprAsciiWriter.hh"
#include "StatPatternRecognition/SprSimpleReader.hh"

#ifdef SPRROOTTUPLE
#include "StatPatternRecognition/SprRootWriter.hh"
#include "StatPatternRecognition/SprRootReader.hh"
#endif

#include <iostream>

class SprPreFilter;


struct SprRWFactory
{
  enum DataType { Ascii, Root };

  static SprAbsReader* makeReader(DataType dataType, 
				  int mode,
				  SprPreFilter* filter=0) {
    if(      dataType == Ascii )
      return new SprSimpleReader(mode,filter);
#ifdef SPRROOTTUPLE
    else if( dataType == Root )
      return new SprRootReader(filter);
#endif
    else {
      std::cerr << "Unknown reader type requested." << std::endl;
      return 0;
    }
  }


  static SprAbsWriter* makeWriter(DataType dataType, const char* label) {
    if(      dataType == Ascii )
      return new SprAsciiWriter(label);
#ifdef SPRROOTTUPLE
    else if( dataType == Root )
      return new SprRootWriter(label);
#endif
    else {
      std::cout << "Unknown writer type requested. " 
		<< "Returning Ascii writer by default." << std::endl;
      return new SprAsciiWriter(label);
    }
  }

};

#endif
