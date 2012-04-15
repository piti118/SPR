//$Id: SprRangeBooster.hh,v 1.3 2008-05-08 19:57:43 narsky Exp $
//
// Description:
//      Class SprRangeBooster :
//         Implements a boosting algorithm that reduces weights of
//         events outside the specified signal efficiency range.
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
 
#ifndef _SprRangeBooster_HH
#define _SprRangeBooster_HH

#include "StatPatternRecognition/SprBagger.hh"
#include "StatPatternRecognition/SprTrainedRangeBooster.hh"

#include <string>
#include <map>

class SprAbsFilter;
class SprAbsTrainedClassifier;
struct SprPoint;


class SprRangeBooster : public SprBagger
{
public:
  struct IndexValueWeight {
    IndexValueWeight()
      : index0(-1), index1(-1), value(0), weight(0) {}

    IndexValueWeight(int i0, int i1, double v, double w)
      : index0(i0), index1(i1), value(v), weight(w) {}

    int index0;
    int index1;
    double value;
    double weight;
  };

  virtual ~SprRangeBooster() {}

  SprRangeBooster(SprAbsFilter* data);

  SprRangeBooster(SprAbsFilter* data, 
		  unsigned cycles, 
		  double signalFraction,
		  double epsilon,
		  double threshold,
		  bool discrete=false);

  /*
    Classifier name.
  */
  std::string name() const { return "RangeBooster"; }

 /*
    Trains classifier on data. Returns true on success, false otherwise.
  */
  bool train(int verbose=0);

  /*
    Reset to untrained state.
  */
  bool reset();

  /*
    Make a trained classifier.
  */
  SprTrainedRangeBooster* makeTrained();

private:
  /*
    Returns -1 on error, 0 on success if can continue training,
    and 1 if reweighting is skipped for this training set.
   */
  int reweight(SprAbsTrainedClassifier* t, 
	       const SprAbsFilter* testData,
	       SprAbsFilter* trainData,
	       const SprAbsFilter* origData,
	       std::map<const SprPoint*,IndexValueWeight>& responses,
	       double acceptSignalWeight,
	       int verbose);

  // signal-like fraction of input events in the focus of this classifier
  double signalFraction_;

  // reweighting factor by which all events outside the range are multiplied
  double epsilon_;

  // weight threshold below which events will be discarded
  double threshold_;

  // notification of stopping to reweight events
  bool notified_;
};

#endif
