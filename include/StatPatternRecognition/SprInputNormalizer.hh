// File and Version Information:
//      $Id: SprInputNormalizer.hh,v 1.2 2008-02-08 23:22:09 narsky Exp $
//
// Description:
//      Class SprInputNormalizer :
//          Normalizes all input variables in data to
//             (X-<X>)/Sx
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2007              California Institute of Technology
//------------------------------------------------------------------------
 
#ifndef _SprInputNormalizer_HH
#define _SprInputNormalizer_HH

#include "StatPatternRecognition/SprAbsVarTransformer.hh"

#include <string>
#include <vector>
#include <iostream>

class SprAbsFilter;


class SprInputNormalizer : public SprAbsVarTransformer
{
public:
  virtual ~SprInputNormalizer() {}

  SprInputNormalizer();
  SprInputNormalizer(const std::vector<double>& mean,
		     const std::vector<double>& sigma);
  SprInputNormalizer(const SprInputNormalizer& other);

  // Clone
  SprInputNormalizer* clone() const {
    return new SprInputNormalizer(*this);
  }

  /*
    VarTransformer name.
  */
  std::string name() const { return "InputNormalizer"; }

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

private:
  std::vector<double> mean_;
  std::vector<double> sigma_;
};

#endif
