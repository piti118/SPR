// File and Version Information:
//      $Id: SprTrainedBinarySplit.hh,v 1.5 2008-01-03 20:51:58 narsky Exp $
//
// Description:
//      Class SprTrainedBinarySplit :
//          Interface for trained classifiers.
//          The purpose of this class is to generate response of 
//          a trained classifier on validation or test data.
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
 
#ifndef _SprTrainedBinarySplit_HH
#define _SprTrainedBinarySplit_HH

#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <iostream>
#include <string>
#include <vector>


class SprTrainedBinarySplit : public SprAbsTrainedClassifier
{
public:
  virtual ~SprTrainedBinarySplit() {}

  SprTrainedBinarySplit(unsigned d, const SprCut& inputCut); 

  SprTrainedBinarySplit(const SprTrainedBinarySplit& other);

  /*
    Returns classifier name.
  */
  std::string name() const {
    return "BinarySplit";
  }

  /*
    Make a clone.
  */
  SprTrainedBinarySplit* clone() const {
    return new SprTrainedBinarySplit(*this);
  }

  /*
    Classifier response for a data point. 
    Binary split produces binary response only: 
    0 for background and 1 for signal.
  */
  double response(const std::vector<double>& v) const;

  /*
    Generate code.
  */
  bool generateCode(std::ostream& os) const {
    return false;
  }

  // Change normalization.
  void useStandard() {}
  void useNormalized() {}
  bool normalized() const { return true; }

  /*
    Print out.
  */
  void print(std::ostream& os) const;

  // Local methods.

  // Return input cut.
  SprCut inputCut() const { return inputCut_; }


private:
  unsigned d_;
  SprCut inputCut_;
};

#endif
