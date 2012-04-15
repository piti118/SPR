// File and Version Information:
//      $Id: SprTrainedRangeBooster.hh,v 1.1 2008-01-03 20:51:58 narsky Exp $
//
// Description:
//      Class SprTrainedRangeBooster :
//         Classifies a new event by using 0-1 response of the first C-1
//         registered classifiers and the full continuous response
//         of the last C-th classifier. It is the user responsibility
//         to make sure that all supplied trained classifiers are normalized
//         to [0,1] interval by applying useNormalized() method if necessary.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2008              California Institute of Technology
//
//------------------------------------------------------------------------
 
#ifndef _SprTrainedRangeBooster_HH
#define _SprTrainedRangeBooster_HH

#include "StatPatternRecognition/SprTrainedBagger.hh"

#include <vector>
#include <utility>
#include <string>


class SprTrainedRangeBooster : public SprTrainedBagger
{
public:
  virtual ~SprTrainedRangeBooster() {}

  SprTrainedRangeBooster(const std::vector<
			 std::pair<const SprAbsTrainedClassifier*,bool> >& 
			 trained);

  SprTrainedRangeBooster(const SprTrainedRangeBooster& other);

  SprTrainedRangeBooster* clone() const {
    return new SprTrainedRangeBooster(*this);
  }

  /*
    Classifier name.
  */
  std::string name() const { return "RangeBooster"; }

  /*
    Classifier response for a data point. 
    Works only for problems with two categories, e.g., signal and background.
  */
  double response(const std::vector<double>& v) const;
};

#endif

