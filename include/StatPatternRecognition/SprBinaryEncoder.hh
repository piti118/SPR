// File and Version Information:
//      $Id: SprBinaryEncoder.hh,v 1.2 2008-05-08 19:57:43 narsky Exp $
//
// Description:
//      Class SprBinaryEncoder :
//         Trains a multi-class learner by encoding multi-class data
//         with D features, K classes and N instances as binary data 
//         with D+1 features, 2 classes and N*K instances. For example,
//         a vector X of class 2 for a problem with 3 classes is encoded
//         as 3 vectors:
//            X,1    X,2    X,3
//         where X,2 is labeled as class 1 and the two others are 
//         labeled as class 0. The disdvantage of this method
//         is the K-fold increase in memory consumption.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2008              California Institute of Technology
//------------------------------------------------------------------------
 
#ifndef _SprBinaryEncoder_HH
#define _SprBinaryEncoder_HH

#include "StatPatternRecognition/SprAbsMultiClassLearner.hh"
#include "StatPatternRecognition/SprTrainedBinaryEncoder.hh"

#include <vector>
#include <iostream>
#include <cassert>
#include <string>

class SprAbsFilter;
class SprEmptyFilter;
class SprAbsClassifier;
class SprAbsTrainedClassifier;


class SprBinaryEncoder : public SprAbsMultiClassLearner
{
public:
  virtual ~SprBinaryEncoder();

  SprBinaryEncoder(SprAbsFilter* data, 
		   SprAbsClassifier* c,
		   const std::vector<int>& classes);

  /*
    Returns classifier name.
  */
  std::string name() const { return "BinaryEncoder"; }

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
    Prints results of training.
  */
  void print(std::ostream& os) const;

  /*
    Make a trained classifier.
  */
  SprTrainedBinaryEncoder* makeTrained();

  /*
    Convert data into binary format. User of this method assumes ownership.
  */
  SprEmptyFilter* convertData() const;

private:
  void destroy();

  SprEmptyFilter* convertedData_;
  SprAbsClassifier* trainable_;
  const SprAbsTrainedClassifier* trained_;
};

#endif
