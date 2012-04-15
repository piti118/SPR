//$Id: SprTrainedAdaBoost.cc,v 1.10 2008-02-22 23:07:50 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprTrainedAdaBoost.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprTransformation.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <stdio.h>
#include <algorithm>
#include <cmath>

using namespace std;


SprTrainedAdaBoost::SprTrainedAdaBoost(const std::vector<
		       std::pair<const SprAbsTrainedClassifier*,bool> >& 
				       trained, 
				       const std::vector<double>& beta,
				       bool useStandard,
				       AdaBoostMode mode) 
  : 
  SprAbsTrainedClassifier(),
  trained_(trained),
  beta_(beta),
  mode_(mode),
  standard_(useStandard),
  epsilon_(0.01),
  nUsedClassifiers_(0)
{
  assert( trained_.size() == beta_.size() );
  assert( !trained_.empty() );
  if( standard_ )
    this->setCut(SprUtils::lowerBound(0.));
  else
    this->setCut(SprUtils::lowerBound(0.5));
}


SprTrainedAdaBoost::SprTrainedAdaBoost(const SprTrainedAdaBoost& other)
  :
  SprAbsTrainedClassifier(other),
  trained_(),
  beta_(other.beta_),
  mode_(other.mode_),
  standard_(other.standard_),
  epsilon_(other.epsilon_),
  nUsedClassifiers_(other.nUsedClassifiers_)
{
  for( int i=0;i<other.trained_.size();i++ )
    trained_.push_back(pair<const SprAbsTrainedClassifier*,bool>
		       (other.trained_[i].first->clone(),true));
}


void SprTrainedAdaBoost::useStandard()
{ 
  standard_ = true;
  this->setCut(SprUtils::lowerBound(0.));
}

void SprTrainedAdaBoost::useNormalized() 
{ 
  standard_ = false; 
  this->setCut(SprUtils::lowerBound(0.5));
}


double SprTrainedAdaBoost::response(const std::vector<double>& v) const
{
  // determine the number of classifiers to be used
  unsigned nClassifiers = 0;
  if( nUsedClassifiers_ > 0 )
    nClassifiers = ( nUsedClassifiers_<trained_.size() ? 
		     nUsedClassifiers_ : trained_.size() );
  else
    nClassifiers = trained_.size();

  // compute output
  double result = 0;
  if(      mode_==Discrete || mode_==Epsilon ) {
    for( int i=0;i<nClassifiers;i++ ) {
      int out = ( trained_[i].first->accept(v) ? 1 : -1 );
      result += out*beta_[i];
    }
  }
  else if( mode_==Real ) {
    double resp = 0;
    for( int i=0;i<nClassifiers;i++ ) {
      resp = trained_[i].first->response(v);
      resp += (1.-2.*resp)*epsilon_;
      if( resp < SprUtils::eps() ) {
	if( standard_ )
	  return -SprUtils::max();
	else
	  return 0;
      }
      if( resp > 1.-SprUtils::eps() ) {
	if( standard_ )
	  return SprUtils::max();
	else
	  return 1;
      }
      result += SprTransformation::logitHalfInverse(resp)*beta_[i];
    }
  }

  // transform to [0,1] if required
  if( !standard_ )
    result = SprTransformation::logitDouble(result);

  // exit
  return result;
}


void SprTrainedAdaBoost::destroy()
{
  for( int i=0;i<trained_.size();i++ ) {
    if( trained_[i].second )
      delete trained_[i].first;
  }
}


void SprTrainedAdaBoost::print(std::ostream& os) const
{
  assert( beta_.size() == trained_.size() );
  os << "Trained AdaBoost " << SprVersion << endl;
  os << "Classifiers: " << trained_.size();
  os << " Cut: " << cut_.size();
  for( int i=0;i<cut_.size();i++ )
    os << " " << cut_[i].first << " " << cut_[i].second;
  os << endl;
  os << "Mode: " << int(mode_) << "   Epsilon: " << epsilon_ << endl;
  for( int i=0;i<trained_.size();i++ ) {
    char s [200];
    sprintf(s,"Classifier %6i %s Beta: %12.10f",
	    i,trained_[i].first->name().c_str(),beta_[i]);
    os << s << endl;
  }
  os << "Classifiers:" << endl;
  for( int i=0;i<trained_.size();i++ ) {
    os << "Classifier " << i 
       << " " << trained_[i].first->name().c_str() << endl;
    trained_[i].first->print(os);
  }
}


bool SprTrainedAdaBoost::generateCode(std::ostream& os) const 
{ 
  // generate weak classifiers
  for( int i=0;i<trained_.size();i++ ) { 
    string name = trained_[i].first->name();
    os << " // Classifier " << i  
       << " \"" << name.c_str() << "\"" << endl; 
    if( !trained_[i].first->generateCode(os) ) {
      cerr << "Unable to generate code for classifier " << name.c_str() 
	   << endl;
      return false;
    }
    if( i < trained_.size()-1 ) os << endl; 
  }

  // exit
  return true;
} 
