//$Id: SprTrainedLogitR.cc,v 1.7 2008-01-03 20:51:59 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprTrainedLogitR.hh"
#include "StatPatternRecognition/SprTransformation.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprUtils.hh"

#include <iomanip>
#include <cassert>

using namespace std;


SprTrainedLogitR::SprTrainedLogitR(double beta0, const SprVector& beta)
  : 
  SprAbsTrainedClassifier(), 
  beta0_(beta0),
  beta_(beta),
  standard_(false)
{
  this->setCut(SprUtils::lowerBound(0.5));
}


SprTrainedLogitR::SprTrainedLogitR(const SprTrainedLogitR& other)
  :
  SprAbsTrainedClassifier(other),
  beta0_(other.beta0_),
  beta_(other.beta_),
  standard_(other.standard_)
{}


double SprTrainedLogitR::response(const std::vector<double>& v) const
{
  // sanity check
  int size = v.size();
  assert( size == beta_.num_row() );

  // compute linear contribution
  double result = 0;
  for( int i=0;i<size;i++ ) {
    result += v[i] * beta_[i];
  }

  // add const term
  result += beta0_;

  // transform
  if( !standard_ ) 
    result = SprTransformation::logit(result);

  // exit
  return result;
}


double SprTrainedLogitR::response(const SprVector& v) const
{
  // sanity check
  assert( v.num_row() == beta_.num_row() );

  // compute linear contribution
  double result = dot(v,beta_);

  // add const term
  result += beta0_;

  // transform
  if( !standard_ ) 
    result = SprTransformation::logit(result);

  // exit
  return result;
}


void SprTrainedLogitR::print(std::ostream& os) const
{
  os << "Trained LogitR " << SprVersion << endl;
  os << "LogitR dimensionality: " << beta_.num_row();
  os << " Cut: " << cut_.size();
  for( int i=0;i<cut_.size();i++ )
    os << " " << cut_[i].first << " " << cut_[i].second;
  os << endl;
  os << "LogitR response: L = Beta0 + Beta*X" << endl;  
  os << "By default logit transform is applied: L <- 1/[1+exp(-L)]" << endl;
  os << "Beta0: " << beta0_ << endl;
  os << "Vector of Beta Coefficients:" << endl;
  for( int i=0;i<beta_.num_row();i++ )
    os << setw(10) << beta_[i] << " ";
  os << endl;
}
