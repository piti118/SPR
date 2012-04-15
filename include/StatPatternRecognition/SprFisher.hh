// File and Version Information:
//      $Id: SprFisher.hh,v 1.6 2008-05-08 19:57:43 narsky Exp $
//
// Description:
//      Class SprFisher :
//          Trains linear or quadratic discriminant.
/*
  For description of linear and quadratic discriminant analysis see, e.g.,
    Hastie, Tibshirani and Friedman "The Elements of Statistical Learning"
*/
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
 
#ifndef _SprFisher_HH
#define _SprFisher_HH

#include "StatPatternRecognition/SprAbsClassifier.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprTrainedFisher.hh"

#include "math/SprSymMatrix.hh"
#include "math/SprVector.hh"

#include <string>
#include <iostream>
#include <utility>
#include <string>

class SprAbsFilter;


class SprFisher : public SprAbsClassifier
{
public:
  virtual ~SprFisher() {}

  SprFisher(SprAbsFilter* data, int mode=1); 

  /*
    Classifier name.
  */
  std::string name() const { return "Fisher"; }

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
  SprTrainedFisher* makeTrained();

  /*
    Choose two classes.
  */
  bool setClasses(const SprClass& cls0, const SprClass& cls1) {
    cls0_ = cls0; cls1_ = cls1;
    return true;
  }

  //
  // Local methods for Fisher discriminant.
  //

  //
  // modifiers
  //
  void setMode(int mode) { mode_ = mode; }

  //
  // accessors
  //
  int mode() const { return mode_; }

  std::string charMode() const {
    if(      mode_ == 1 )
      return "linear";
    else if( mode_ == 2 )
      return "quadratic";
    return "unknown";
  }

  void linear(SprVector& v) const { v = linear_; }
  void quadratic(SprSymMatrix& m) const { m = quadr_; }
  double cterm() const { return cterm_; }

private:
  void setClasses();// copies two classes from the filter

  int mode_;// Fisher mode: 1 = linear, 2 = quadratic
  SprClass cls0_;// class indices in the data
  SprClass cls1_;
  SprVector linear_;// linear coefficients
  SprSymMatrix quadr_;// quadratic coefficients
  double cterm_;// const term for proper normalization
};

#endif
