//------------------------------------------------------------------------
// File and Version Information:
//      $Id: SprEmptyFilter.hh,v 1.5 2008-01-30 21:27:38 narsky Exp $
//
// Description:
//      Class SprEmptyFilter :
//         No filter applied.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2005              California Institute of Technology
//
//------------------------------------------------------------------------
 
#ifndef _SprEmptyFilter_HH
#define _SprEmptyFilter_HH

#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprClass.hh"

#include <vector>

struct SprPoint;
class SprData;


class SprEmptyFilter : public SprAbsFilter
{
public:
  virtual ~SprEmptyFilter() {}

  SprEmptyFilter(const SprData* data, bool ownData=false) 
    : SprAbsFilter(data,ownData) {}

  SprEmptyFilter(const SprData* data, 
		 const std::vector<SprClass>& classes, 
		 bool ownData=false) 
    : SprAbsFilter(data,classes,ownData) {}

  SprEmptyFilter(const SprData* data, 
		 const std::vector<double>& weights,
		 bool ownData=false) 
    : SprAbsFilter(data,weights,ownData) {}

  SprEmptyFilter(const SprData* data, 
		 const std::vector<SprClass>& classes, 
		 const std::vector<double>& weights,
		 bool ownData=false) 
    : SprAbsFilter(data,classes,weights,ownData) {}

  SprEmptyFilter(const SprEmptyFilter& filter) 
    : SprAbsFilter(filter) {}

  SprEmptyFilter(const SprAbsFilter* filter)
    : SprAbsFilter(*filter) {}

  // specific reset
  bool reset() { return true; }

  // accept or reject a point
  bool pass(const SprPoint* p) const { return true; }
};

#endif

