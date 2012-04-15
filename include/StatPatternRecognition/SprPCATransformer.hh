// File and Version Information:
//      $Id: SprPCATransformer.hh,v 1.2 2008-02-08 23:22:09 narsky Exp $
//
// Description:
//      Class SprPCATransformer :
//          Principal Component Analysis
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
 
#ifndef _SprPCATransformer_HH
#define _SprPCATransformer_HH

#include "StatPatternRecognition/SprAbsVarTransformer.hh"
#include "math/SprMatrix.hh"

#include <string>
#include <vector>
#include <utility>
#include <iostream>

class SprAbsFilter;


class SprPCATransformer : public SprAbsVarTransformer
{
public:
  virtual ~SprPCATransformer() {}

  SprPCATransformer();
  SprPCATransformer(const SprMatrix& U, 
		    const std::vector<std::pair<double,int> >& eigenValues); 
  SprPCATransformer(const SprPCATransformer& other);

  // Clone
  SprPCATransformer* clone() const {
    return new SprPCATransformer(*this);
  }

  /*
    VarTransformer name.
  */
  std::string name() const { return "PCA"; }

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
    return (!oldVars_.empty() && oldVars_.size()==U_.num_row());
  }

  // Returns true if all transformations for individual variables
  // are independent of each other.
  bool allVarsIndependent() const { return false; }

  // Eliminates transformation variables which are not in the supplied list.
  bool reduceVars(const std::vector<std::string>& vars) { return true; }

  // Output
  void print(std::ostream& os) const;

private:
  // Transformation matrix.
  // U is orthogonal and its rows give eigenvectors.
  // S = T(U)*D*U, where S is the covariance matrix for the data,
  //   and D is the diagonal matrix of eigenvalues.
  SprMatrix U_;

  // eigenvalues sorted in descending order with indices pointing to data vars
  std::vector<std::pair<double,int> > eigenValues_;
};

#endif
