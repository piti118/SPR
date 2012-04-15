//$Id: SprTrainedBagger.cc,v 1.8 2008-01-03 20:51:59 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprTrainedBagger.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <stdio.h>
#include <cassert>

using namespace std;


SprTrainedBagger::SprTrainedBagger(const std::vector<
		   std::pair<const SprAbsTrainedClassifier*,bool> >& 
		   trained, bool discrete) 
  : 
  SprAbsTrainedClassifier(),
  trained_(trained),
  discrete_(discrete),
  nUsedClassifiers_(0)
{
  assert( !trained_.empty() );
  this->setCut(SprUtils::lowerBound(0.5));
}


SprTrainedBagger::SprTrainedBagger(const SprTrainedBagger& other)
  :
  SprAbsTrainedClassifier(other),
  trained_(),
  discrete_(other.discrete_),
  nUsedClassifiers_(other.nUsedClassifiers_)
{
  for( int i=0;i<other.trained_.size();i++ )
    trained_.push_back(pair<const SprAbsTrainedClassifier*,bool>
		       (other.trained_[i].first->clone(),true));
}


double SprTrainedBagger::response(const std::vector<double>& v) const
{
  // determine the number of classifiers to be used
  unsigned nClassifiers = 0;
  if( nUsedClassifiers_ > 0 )
    nClassifiers = ( nUsedClassifiers_<trained_.size() ? 
		     nUsedClassifiers_ : trained_.size() );
  else
    nClassifiers = trained_.size();

  // init
  double r = 0;

  // discrete/continuous
  if( discrete_ ) {
    int out = 0;
    for( int i=0;i<nClassifiers;i++ )
      out += ( trained_[i].first->accept(v) ? 1 : -1 );
    r = out;
    r /= 2.*nClassifiers;
    r += 0.5;
  }
  else {
    for( int i=0;i<nClassifiers;i++ )
      r += trained_[i].first->response(v);
    r /= nClassifiers;
  }

  // exit
  return r;
}


void SprTrainedBagger::destroy()
{
  for( int i=0;i<trained_.size();i++ ) {
    if( trained_[i].second )
      delete trained_[i].first;
  }
}


void SprTrainedBagger::print(std::ostream& os) const
{
  os << "Trained " << this->name() << " " << SprVersion << endl;
  os << "Classifiers: " << trained_.size();
  os << " Cut: " << cut_.size();
  for( int i=0;i<cut_.size();i++ )
    os << " " << cut_[i].first << " " << cut_[i].second;
  os << endl;
  for( int i=0;i<trained_.size();i++ ) {
    os << "Classifier " << i 
       << " " << trained_[i].first->name().c_str() << endl;
    trained_[i].first->print(os);
  }
}


bool SprTrainedBagger::generateCode(std::ostream& os) const 
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


SprTrainedBagger& SprTrainedBagger::operator+=(const SprTrainedBagger& other)
{
  // check vars
  if( vars_.size() != other.vars_.size() ) {
    cerr << "Unable to add Bagger: variable lists do not match." << endl;
    return *this;
  }
  for( int i=0;i<vars_.size();i++ ) {
    if( vars_[i] != other.vars_[i] ) {
      cerr << "Unable to add Bagger: variable lists do not match." << endl;
      cerr << "Variables " << i << ": " 
	   << vars_[i] << " " << other.vars_[i] << endl;
      return *this;
    }
  }

  // check discreteness
  if( discrete_ != other.discrete_ ) {
    cerr << "Unable to add Bagger: discreteness does not match." << endl;
    return *this;
  }

  // add
  for( int i=0;i<other.trained_.size();i++ ) {
    trained_.push_back(pair<const SprAbsTrainedClassifier*,
		       bool>(other.trained_[i].first->clone(),true));
  }
  this->setCut(SprUtils::lowerBound(0.5));

  // exit
  return *this;
}


const SprTrainedBagger operator+(const SprTrainedBagger& l,
				 const SprTrainedBagger& r)
{
  // check variable list
  assert( l.vars_.size() == r.vars_.size() );
  for( int i=0;i<l.vars_.size();i++ )
    assert( l.vars_[i] == r.vars_[i] );

  // add classifiers
  vector<pair<const SprAbsTrainedClassifier*,bool> > trained;
  for( int i=0;i<l.trained_.size();i++ ) {
    trained.push_back(pair<const SprAbsTrainedClassifier*,
		      bool>(l.trained_[i].first->clone(),true));
  }
  
  for( int i=0;i<r.trained_.size();i++ ) {
    trained.push_back(pair<const SprAbsTrainedClassifier*,
		      bool>(r.trained_[i].first->clone(),true));
  }

  // add discrete
  assert( l.discrete_ == r.discrete_ );

  // make bagger and set cut
  SprTrainedBagger newBagger(trained,l.discrete_);
  newBagger.setCut(SprUtils::lowerBound(0.5));

  // exit
  return newBagger;
}
