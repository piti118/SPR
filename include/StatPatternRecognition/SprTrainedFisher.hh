// File and Version Information:
//      $Id: SprTrainedFisher.hh,v 1.8 2008-02-22 23:07:50 narsky Exp $
//
// Description:
//      Class SprTrainedFisher :
//          Returns response of a trained Fisher.
//          The quantity returned is
//    log(f1/f0) = log(N1/N0) -(1/2)*log(|Cov1|/|Cov0|) 
//                 -(1/2)*{m1*Cov1*m1-m0*Cov0*m0}
//                 +x*{Cov1^-1*mu1-Cov0^-1*mu0}
//                 -(1/2)*x*{Cov1^-1-Cov0^-1}*x
//          where the first two lines are a constant term 
//                needed for proper normalization,
//                the 3rd line is the linear term,
//                and the 4th line is the quadratic term.
//          Under assumption Cov1=Cov0, the quadratic term disappears
//          and we obtain the conventional linear Fisher discriminant.
//
//    Notation:
//       N0, N1 - training sample sizes
//       d      - dimensionality of input space
//       f0, f1 - densities
//       m0, m1 - means
//       Cov0, Cov1 - covariance matrices
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
 
#ifndef _SprTrainedFisher_HH
#define _SprTrainedFisher_HH

#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "math/SprVector.hh"
#include "math/SprSymMatrix.hh"

#include <vector>
#include <iostream>
#include <string>


class SprTrainedFisher : public SprAbsTrainedClassifier
{
public:
  virtual ~SprTrainedFisher() {}

  SprTrainedFisher(const SprVector& v, double cterm);

  SprTrainedFisher(const SprVector& v, const SprSymMatrix& m, double cterm);

  SprTrainedFisher(const SprTrainedFisher& other);

  SprTrainedFisher* clone() const {
    return new SprTrainedFisher(*this);
  }

  // inherited methods
  std::string name() const { return "Fisher"; }
  double response(const std::vector<double>& v) const;
  double response(const SprVector& v) const;
  void print(std::ostream& os) const;
  bool generateCode(std::ostream& os) const {
    return false;
  }

  //
  // local methods
  //
  int mode() const {
    if( quadr_.num_row() > 0 ) return 2;
    return 1;
  }
  void linear(SprVector& v) const { v = linear_; }
  void quadratic(SprSymMatrix& m) const { m = quadr_; }
  double cterm() const { return cterm_; }

  /*
    The useStandard flag lets you switch from the standard version to
    the "normalized" version. The standard version returns normal
    L(Q)DA output ranging from -infty to +infty. The normalized
    version applied logit transform 1/[1+exp(-F)] to this output.
  */
  void useStandard();
  void useNormalized();
  bool normalized() const { return !standard_; }

private:
  SprVector linear_;// linear coefficients
  SprSymMatrix quadr_;// quadratic coefficients
  double cterm_;// const term in the Fisher expression
  bool standard_;
};

#endif
