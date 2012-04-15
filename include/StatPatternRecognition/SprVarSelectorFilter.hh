//------------------------------------------------------------------------
// File and Version Information:
//      $Id: SprVarSelectorFilter.hh,v 1.2 2008-01-30 21:27:38 narsky Exp $
//
// Description:
//      Class SprVarSelectorFilter :
//         Chooses a subset of input variables.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2007              California Institute of Technology
//
//------------------------------------------------------------------------
 
#ifndef _SprVarSelectorFilter_HH
#define _SprVarSelectorFilter_HH

#include "StatPatternRecognition/SprAbsFilter.hh"

#include <vector>
#include <set>
#include <string>

struct SprPoint;
class SprData;


class SprVarSelectorFilter : public SprAbsFilter
{
public:
  virtual ~SprVarSelectorFilter() {}

  SprVarSelectorFilter(const SprData* data, 
		       bool ownData=false) 
    : SprAbsFilter(data,ownData) {}

  SprVarSelectorFilter(const SprVarSelectorFilter& filter) 
    : SprAbsFilter(filter), vars_(filter.vars_) {}

  SprVarSelectorFilter(const SprAbsFilter* filter)
    : SprAbsFilter(*filter) {}

  // specific reset
  bool reset() { 
    vars_.clear();
    return true;
  }

  // reimplements filtering to reduce the number of input variables
  bool filter();

  // accept or reject a point
  bool pass(const SprPoint* p) const { return true; }

  // Local methods
  bool chooseVars(const std::vector<std::string>& vars) {
    vars_ = vars;
  }
  bool chooseVars(const std::set<std::string>& vars);

private:
  std::vector<std::string> vars_;
};

#endif

