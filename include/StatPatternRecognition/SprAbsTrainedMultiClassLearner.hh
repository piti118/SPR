// File and Version Information:
//      $Id: SprAbsTrainedMultiClassLearner.hh,v 1.5 2008-04-02 23:36:43 narsky Exp $
//
// Description:
//      Class SprAbsTrainedMultiClassLearner :
//          Interface for trained multiclass methods.
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
 
#ifndef _SprAbsTrainedMultiClassLearner_HH
#define _SprAbsTrainedMultiClassLearner_HH

#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <vector>
#include <utility>
#include <iostream>
#include <string>
#include <map>

class SprAbsFilter;


class SprAbsTrainedMultiClassLearner
{
public:
  typedef std::pair<int,std::vector<double> > MissingType;

  virtual ~SprAbsTrainedMultiClassLearner() {}

  SprAbsTrainedMultiClassLearner();

  SprAbsTrainedMultiClassLearner(const std::vector<int>& classes);

  SprAbsTrainedMultiClassLearner(const SprAbsTrainedMultiClassLearner& other);

  /*
    Returns classifier name.
  */
  virtual std::string name() const = 0;

  /*
    Make a clone.
  */
  virtual SprAbsTrainedMultiClassLearner* clone() const = 0;

  /*
    Classifier response for a data point. 
    Computes loss values for registered classifiers; the map key
    is the class and the value is the corresponding loss.
    The returned integer is the class for which the loss is minimal.
  */
  int response(const std::vector<double>& input,
	       std::map<int,double>& output) const;
  int response(const std::vector<double>& input) const {
    std::map<int,double> output;
    return this->response(input,output);
  }
  int response(const SprPoint* input, 
	       std::map<int,double>& output) const {
    return this->response(input->x_,output);
  }
  int response(const SprPoint* input) const {
    std::map<int,double> output;
    return this->response(input->x_,output);
  }

  /*
    Returns output names associated with map<int,double> 
    returned by response() methods.
  */
  void outputNames(const char* prefix,
		   std::vector<std::string>& names) const;

  // Access to the list of variables.
  void setVars(const std::vector<std::string>& vars) { vars_ = vars; }
  void vars(std::vector<std::string>& vars) const { vars = vars_; }
  unsigned dim() const { return vars_.size(); }

  /*
    Responses with missing values.
  */
  bool setDefaultMissing(const SprCut& validRange,
			 const std::vector<MissingType>& defaultMissing);
  void clearDefaultMissing() {
    validRange_.clear();
    defaultMissing_.clear();
  }

  /*
    Print out.
  */
  bool store(const char* filename) const;
  virtual void print(std::ostream& os) const = 0;

  // returns number of categories
  unsigned nClasses() const { return mapper_.size(); }

  // returns categories
  void classes(std::vector<int>& classes) const {
    classes = mapper_;
  }

protected:
  // names in the absence of missing values
  void outputNamesWithoutMissing(const char* prefix,
				 std::vector<std::string>& names) const;

  // response in the absence of missing values
  virtual int response_one(const std::vector<double>& input,
			   std::map<int,double>& output) const = 0;

  // response with missing values handled for each input class separately
  int response_many(const std::vector<double>& input,
		    std::map<int,double>& output) const;

  // variable list
  std::vector<std::string> vars_;

  // list of integer classes
  std::vector<int> mapper_;

  // for handling missing values
  SprCut validRange_;
  std::vector<MissingType> defaultMissing_;
};

#endif
