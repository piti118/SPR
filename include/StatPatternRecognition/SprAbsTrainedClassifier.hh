// File and Version Information:
//      $Id: SprAbsTrainedClassifier.hh,v 1.8 2008-01-03 20:51:58 narsky Exp $
//
// Description:
//      Class SprAbsTrainedClassifier :
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
 
#ifndef _SprAbsTrainedClassifier_HH
#define _SprAbsTrainedClassifier_HH

#include "StatPatternRecognition/SprVector.hh"
#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <vector>
#include <iostream>
#include <string>

class SprAbsFilter;


class SprAbsTrainedClassifier
{
public:
  virtual ~SprAbsTrainedClassifier() {}

  SprAbsTrainedClassifier() : cut_(), vars_() {}

  SprAbsTrainedClassifier(const SprAbsTrainedClassifier& other)
    : cut_(other.cut_), vars_(other.vars_) {}

  /*
    Returns classifier name.
  */
  virtual std::string name() const = 0;

  /*
    Make a clone.
  */
  virtual SprAbsTrainedClassifier* clone() const = 0;

  /*
    Classifier response for a data point. 
    Works only for problems with two categories, e.g., signal and background.
  */
  virtual double response(const std::vector<double>& v) const = 0;
  double response(const SprVector& v) const { 
    return this->response(v.std()); 
  }
  double response(const SprPoint* p) const {
    return this->response(p->x_);
  }

  /*
    Generate code.
  */
  virtual bool generateCode(std::ostream& os) const = 0;
  bool storeCode(const char* filename) const;

  /*
    Classifier response. These methods set a cut on the classifier 1D output
    and accept a point if it satisfies these cuts.
  */
  virtual void setCut(const SprCut& cut) { cut_ = cut; }
  virtual void resetCut(const SprCut& cut) { cut_.clear(); }
  virtual SprCut cut() const { return cut_; }
  virtual bool accept(const std::vector<double>& v) const {
    double r = 0;
    return this->accept(v,r);
  }
  virtual bool accept(const SprVector& v) const {
    double r = 0;
    return this->accept(v,r);
  }
  bool accept(const SprPoint* p) const {
    return this->accept(p->x_);
  }
  virtual bool accept(const std::vector<double>& v, double& response) const;
  bool accept(const SprVector& v, double& response) const {
    return this->accept(v.std(),response);
  }
  bool accept(const SprPoint* p, double& response) const {
    return this->accept(p->x_,response);
  }

  /*
    The useNormalized() method forces the classifier output range to
    be in the range from 0 to 1. This method only has an effect for
    classifiers which by default have their output range from -infty
    to +infty (or any other range different from [0,1]). This method
    cannot be used to adjust the output range of a classifier if its
    output by accident does not follow the [0,1] range (for example,
    if the output clumps around 0.5 even though in principle it is
    allowed to extend to 0 and 1.)  useStandard() switches the
    classifier output back to the default range. normalized() tells
    you if the range is currently normalized.
  */
  virtual void useNormalized() = 0;
  virtual void useStandard() = 0;
  virtual bool normalized() const = 0;

  /*
    Variables used for these trained classifier.
    The list of variables can be set optionally.
    It is up to the user to figure out how s/he wants to use the list
    of variables.
    A typical application would be to read the trained classifier
    configuration from a file and then look at the list of variables
    to make sure they make sense.
  */
  void setVars(const std::vector<std::string>& vars) { vars_ = vars; }
  void vars(std::vector<std::string>& vars) const { vars = vars_; }
  unsigned dim() const { return vars_.size(); }

  /*
    Print out.
  */
  bool store(const char* filename) const;
  virtual void print(std::ostream& os) const = 0;
  std::vector<std::string>* getVars() const;
protected:
  SprCut cut_;// cut imposed on the classifier output for accept()
  std::vector<std::string> vars_;
};

#endif
