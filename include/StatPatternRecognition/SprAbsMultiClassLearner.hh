// File and Version Information:
//      $Id: SprAbsMultiClassLearner.hh,v 1.5 2008-05-08 19:57:43 narsky Exp $
//
// Description:
//      Class SprAbsMultiClassLearner :
//          Interface for multiclass methods.
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
 
#ifndef _SprAbsMultiClassLearner_HH
#define _SprAbsMultiClassLearner_HH

#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"

#include <iostream>
#include <vector>
#include <string>

class SprAbsFilter;


class SprAbsMultiClassLearner
{
public:
  typedef SprAbsTrainedMultiClassLearner::MissingType MissingType;

  virtual ~SprAbsMultiClassLearner();

  SprAbsMultiClassLearner(SprAbsFilter* data,
			  const std::vector<int>& classes);

  /*
    Classifier name.
  */
  virtual std::string name() const = 0;

  /*
    Trains classifier on data. Returns true on success, false otherwise.
  */
  virtual bool train(int verbose=0) = 0;

  /*
    Reset this classifier to untrained state.
  */
  virtual bool reset() = 0;

  /*
    Replace training data.
  */
  virtual bool setData(SprAbsFilter* data) = 0;

  /*
    Prints results of training.
  */
  virtual void print(std::ostream& os) const = 0;

  /*
    Store training results in a file.
  */
  virtual bool store(const char* filename) const;

  /*
    Make a trained classifier.
  */
  virtual SprAbsTrainedMultiClassLearner* makeTrained() = 0;

  /*
    Responses with missing values.
  */
  bool setDefaultMissing(const SprCut& validRange,
		 const std::vector<MissingType>& defaultMissing);
  void clearDefaultMissing() {
    validRange_.clear();
    defaultMissing_.clear();
  }

  // Get classes.
  void classes(std::vector<int>& classes) const {
    classes = mapper_;
  }

  // Check input class list for consistency.
  bool checkClasses() const;

protected:
  // constructs new training data with replaced missng values
  bool replaceMissing();

  SprAbsFilter* data_;// non-const filter to allow adjustment of weights
  std::vector<int> mapper_;// vector of classes

  // for handling missing values
  SprCut validRange_;
  std::vector<MissingType> defaultMissing_;

  // own data?
  bool ownData_;
};

#endif
