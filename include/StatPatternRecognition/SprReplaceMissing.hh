// File and Version Information:
//      $Id: SprReplaceMissing.hh,v 1.1 2008-04-02 23:36:44 narsky Exp $
//
// Description:
//      Class SprReplaceMissing :
//          Replaces values outside of the valid range with supplied defaults
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2008              California Institute of Technology
//------------------------------------------------------------------------
 
#ifndef _SprReplaceMissing_HH
#define _SprReplaceMissing_HH

#include "StatPatternRecognition/SprAbsVarTransformer.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <string>
#include <vector>
#include <utility>
#include <iostream>

class SprAbsFilter;


class SprReplaceMissing : public SprAbsVarTransformer
{
public:
  typedef std::pair<SprClass,std::vector<double> > ClassAndDefaultValues;

  enum Mode { Median, Average };

  virtual ~SprReplaceMissing() {}

  SprReplaceMissing();

  SprReplaceMissing(Mode mode, const SprCut& validRange, bool classBlind);

  SprReplaceMissing(const SprCut& validRange, bool classBlind);

  SprReplaceMissing(const SprReplaceMissing& other);

  // Clone
  SprReplaceMissing* clone() const {
    return new SprReplaceMissing(*this);
  }

  /*
    VarTransformer name.
  */
  std::string name() const { return "ReplaceMissing"; }

  /*
    Computes transformation using input data. 
    Returns true on success, false otherwise.
  */
  bool train(const SprAbsFilter* data, int verbose=0);

  // Applies transformation.
  void transform(const std::vector<double>& in, std::vector<double>& out) 
    const;

  // Applies inverse transformation.
  void inverse(const std::vector<double>& in, std::vector<double>& out) const;

  // Status of the transformer - if returns true, ready to transform.
  bool ready() const {
    return (!oldVars_.empty() && oldVars_.size()==newVars_.size());
  }

  // Returns true if all transformations for individual variables
  // are independent of each other.
  bool allVarsIndependent() const;

  // Eliminates transformation variables which are not in the supplied list.
  bool reduceVars(const std::vector<std::string>& vars);

  // Output
  void print(std::ostream& os) const;

  // Choose class for transform() method.
  bool chooseClass(const SprClass& cls);

  // Set replacement values
  void setReplacement(const std::vector<ClassAndDefaultValues>& replacement) {
    replacement_ = replacement;
  }

  // Get replacement values.
  void replacement(std::vector<ClassAndDefaultValues>& replacement) {
    replacement = replacement_;
  }

private:
  typedef std::vector<std::pair<double,double> > ValuesWithWeights;

  Mode mode_;
  SprCut validRange_;
  bool classBlind_;// true if replacing values do not depend on class
  std::vector<ClassAndDefaultValues> replacement_;
  int chosenClass_;// one chosen class for transform()
};

#endif
