//$Id: SprTrainedFisher.cc,v 1.9 2008-02-22 23:07:50 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprTrainedFisher.hh"
#include "StatPatternRecognition/SprTransformation.hh"
#include "StatPatternRecognition/SprUtils.hh"

#include <iomanip>
#include <cassert>

using namespace std;


SprTrainedFisher::SprTrainedFisher(const SprVector& v, double cterm)
  : 
  SprAbsTrainedClassifier(), 
  linear_(v), 
  quadr_(), 
  cterm_(cterm),
  standard_(false)
{
  this->setCut(SprUtils::lowerBound(0.5));
}


SprTrainedFisher::SprTrainedFisher(const SprVector& v, 
				   const SprSymMatrix& m, 
				   double cterm)
  : 
  SprAbsTrainedClassifier(), 
  linear_(v), 
  quadr_(m), 
  cterm_(cterm),
  standard_(false)
{
  this->setCut(SprUtils::lowerBound(0.5));
}


SprTrainedFisher::SprTrainedFisher(const SprTrainedFisher& other)
  :
  SprAbsTrainedClassifier(other),
  linear_(other.linear_),
  quadr_(other.quadr_),
  cterm_(other.cterm_),
  standard_(other.standard_)
{}


void SprTrainedFisher::useStandard()
{ 
  standard_ = true;
  this->setCut(SprUtils::lowerBound(0.));
}

void SprTrainedFisher::useNormalized() 
{ 
  standard_ = false; 
  this->setCut(SprUtils::lowerBound(0.5));
}


double SprTrainedFisher::response(const std::vector<double>& v) const
{
  // sanity check
  int size = v.size();
  assert( size == linear_.num_row() );

  // compute linear contribution
  double d = 0;
  for( int i=0;i<size;i++ ) {
    d += v[i] * linear_[i];
  }

  // compute quadratic contribution
  if( this->mode() == 2 ) {
    assert( size == quadr_.num_row() );
    double row = 0;
    for( int i=1;i<size;i++ ) {
      row = 0;
      for( int j=0;j<i;j++ )
	row += quadr_[i][j] * v[j];
      d += v[i] * row;
    }
    d *= 2;
    for( int i=0;i<size;i++ )
      d += v[i]*v[i] * quadr_[i][i];
  }

  // add const term
  d += cterm_;

  // apply transform
  if( !standard_ )
    d = SprTransformation::logit(d);

  // exit
  return d;
}


double SprTrainedFisher::response(const SprVector& v) const
{
  // sanity check
  assert( v.num_row() == linear_.num_row() );

  // compute linear contribution
  double d = dot(v,linear_);

  // compute quadratic contribution
  if( this->mode() == 2 ) {
    assert( v.num_row() == quadr_.num_row() );
    d += dot(v,quadr_*v);
  }

  // add const term
  d += cterm_;

  // apply transform
  if( !standard_ )
    d = SprTransformation::logit(d);

  // exit
  return d;
}


void SprTrainedFisher::print(std::ostream& os) const
{
  os << "Trained Fisher " << SprVersion << endl;
  os << "Fisher dimensionality: " << linear_.num_row();
  os << " Cut: " << cut_.size();
  for( int i=0;i<cut_.size();i++ )
    os << " " << cut_[i].first << " " << cut_[i].second;
  os << endl;
  os << "Fisher response: F = C + T(L)*X + T(X)*Q*X; T is transposition" 
     << endl;
  os << "By default logit transform is applied: F <- 1/[1+exp(-F)]" << endl;
  os << "Fisher order: " << (quadr_.num_row()>0 ? 2 : 1) << endl;
  os << "Const term: " << cterm_ << endl;
  os << "Linear Part:" << endl;
  for( int i=0;i<linear_.num_row();i++ )
    os << setw(10) << linear_[i] << " ";
  os << endl;
  if( quadr_.num_row() > 0 ) {
    os << "Quadratic Part:" << endl;
    for( int i=0;i<quadr_.num_row();i++ ) {
      for( int j=0;j<quadr_.num_col();j++ ) {
        os << setw(10) << quadr_[i][j] << " ";
      }
      os << endl;
    }
  }
}
