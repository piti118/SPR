// File and Version Information:
//      $Id: SprMultiClassLearner.hh,v 1.7 2008-05-08 19:57:43 narsky Exp $
//
// Description:
//      Class SprMultiClassLearner :
//          Implements a multi class learning algorithm described
//          in Allwein, Schapire and Singer, 
//          J. of Machine Learning Research 2000,
//          "Reducing multiclass to binary"
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
 
#ifndef _SprMultiClassLearner_HH
#define _SprMultiClassLearner_HH

#include "StatPatternRecognition/SprAbsMultiClassLearner.hh"
#include "StatPatternRecognition/SprTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprMatrix.hh"

#include <vector>
#include <utility>
#include <iostream>
#include <cassert>
#include <string>

class SprAbsFilter;
class SprAbsClassifier;
class SprAbsTrainedClassifier;


class SprMultiClassLearner : public SprAbsMultiClassLearner
{
public:
  enum MultiClassMode { User, OneVsAll, OneVsOne };

  virtual ~SprMultiClassLearner();

  SprMultiClassLearner(SprAbsFilter* data, 
		       SprAbsClassifier* c,
		       const std::vector<int>& classes,
		       const SprMatrix& indicator,
		       MultiClassMode mode=User);

  /*
    Returns classifier name.
  */
  std::string name() const { return "MultiClassLearner"; }

  /*
    Trains classifier on data. Returns true on success, false otherwise.
  */
  bool train(int verbose=0);

  /*
    Reset this classifier to untrained state.
  */
  bool reset();

  /*
    Replace training data.
  */
  bool setData(SprAbsFilter* data);

  /*
    Include non-participating classes in loss evaluation?
  */
  void excludeZeroClasses() { includeZeroIndicator_ = false; }
  void includeZeroClasses() { includeZeroIndicator_ = true; }

  /*
    Set weights for classifier normalization. By default all weights are 1.
    The length of the vector must be equal to the number of columns
    in the indicator matrix.
  */
  void setClassifierWeights(const std::vector<double>& weights);

  /*
    Prints results of training.
  */
  void print(std::ostream& os) const;

  /*
    Make a trained classifier.
  */
  SprTrainedMultiClassLearner* makeTrained();

  /*
    Set parameters.
  */
  void setTrained(const SprMatrix& indicator, 
		  const std::vector<int>& classes,
		  const std::vector<
		  std::pair<const SprAbsTrainedClassifier*,bool> >& trained,
		  const std::vector<double>& weights);

  // Show indicator matrix.
  void printIndicatorMatrix(std::ostream& os) const;
  void indicatorMatrix(SprMatrix& indicator) const {
    indicator = indicator_;
  }

private:
  bool setClasses();
  void destroy();

  MultiClassMode mode_;
  SprMatrix indicator_;
  SprAbsClassifier* trainable_;
  std::vector<std::pair<const SprAbsTrainedClassifier*,bool> > trained_;
  std::vector<double> weights_;// weights for trained classifiers
  bool keepWeights_;// keep user-supplied weights or recompute
  bool includeZeroIndicator_;
};

#endif
