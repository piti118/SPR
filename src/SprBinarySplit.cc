//$Id: SprBinarySplit.cc,v 1.11 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprBinarySplit.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprAbsTwoClassCriterion.hh"

#include <stdio.h>
#include <iostream>
#include <cassert>
#include <utility>
#include <algorithm>
#include <functional>
#include <cmath>

using namespace std;


struct SBSCmpPairFirst 
  : public binary_function<pair<double,int>,pair<double,int>,bool> {
  bool operator()(const pair<double,int>& l, const pair<double,int>& r)
    const {
    return (l.first < r.first);
  }
};


SprBinarySplit::SprBinarySplit(SprAbsFilter* data, 
			       const SprAbsTwoClassCriterion* crit,
			       unsigned d) 
  : 
  SprAbsClassifier(data), 
  crit_(crit),
  d_(d),
  cls0_(0), 
  cls1_(1), 
  cut_(),
  nSorted_(0),
  sorted0_(),
  sorted1_(),
  division_()
{
  assert( crit_ != 0 );
  assert( d_ < data_->dim() );
  this->setClasses();
}


void SprBinarySplit::setClasses() 
{
  vector<SprClass> classes;
  data_->classes(classes);
  int size = classes.size();
  if( size > 0 ) cls0_ = classes[0];
  if( size > 1 ) cls1_ = classes[1];
  //  cout << "Classes for SprBinarySplit are set to " 
  //       << cls0_ << " " << cls1_ << endl;
}


bool SprBinarySplit::setClasses(const SprClass& cls0, const SprClass& cls1)
{
  cls0_ = cls0; cls1_ = cls1;
  //  cout << "Classes for binary split reset to " 
  //       << cls0_ << " " << cls1_ << std::endl;
  return true;
}


bool SprBinarySplit::reset()
{
  cut_.clear();
  return true;
}


bool SprBinarySplit::setData(SprAbsFilter* data)
{
  assert( data != 0 );
  data_ = data;
  sorted0_.clear();
  sorted1_.clear();
  division_.clear();
  nSorted_ = 0;
  return this->reset();
}


SprTrainedBinarySplit* SprBinarySplit::makeTrained()
{
  if( d_ >= data_->dim() ) return 0;
  return new SprTrainedBinarySplit(d_,cut_);
}

bool SprBinarySplit::train(int verbose)
{
  // if dimension is greater than dimensionality, do nothing
  if( d_ >= data_->dim() ) return true;

  // sort
  if( !this->sort() ) {
    cerr << "SprBinarySplit unable to sort data." << endl;
    return false;
  }

  // init weights
  double wcor0(0), wmis1(0);
  double wmis0(0), wcor1(0);
  for( int j=0;j<sorted0_.size();j++ )
    wmis0 += data_->w(sorted0_[j]);
  for( int j=0;j<sorted1_.size();j++ )
    wcor1 += data_->w(sorted1_[j]);
  assert( wmis0>0 && wcor1>0 );

  // prepare vars
  int i0start(0), i0split(0);
  int i1start(0), i1split(0);
  int ndiv = division_.size();
  vector<double> flo(ndiv), fhi(ndiv);
  bool lbreak = false;// found breaking point

  // loop through divisions
  for( int k=0;k<ndiv;k++ ) {
    double z = division_[k];
    lbreak = false;
    for( i0split=i0start;i0split<sorted0_.size();i0split++ ) {
      if( ((*data_)[sorted0_[i0split]])->x_[d_] > z ) {
	lbreak = true;
	break;
      }
    }
    if( !lbreak ) i0split = sorted0_.size();
    lbreak = false;
    for( i1split=i1start;i1split<sorted1_.size();i1split++ ) {
      if( ((*data_)[sorted1_[i1split]])->x_[d_] > z ) {
	lbreak = true;
	break;
      }
    }
    if( !lbreak ) i1split = sorted1_.size();
    for( int i=i0start;i<i0split;i++ ) {
      double w = data_->w(sorted0_[i]);
      wcor0 += w;
      wmis0 -= w;
    }
    for( int i=i1start;i<i1split;i++ ) {
      double w = data_->w(sorted1_[i]);
      wmis1 += w;
      wcor1 -= w;
    }
    i0start = i0split;
    i1start = i1split;
    flo[k] = crit_->fom(wcor0,wmis0,wcor1,wmis1);
    fhi[k] = crit_->fom(wmis0,wcor0,wmis1,wcor1);
  }
  
  // find minimal point and sign of cut
  vector<double>::iterator ilo = max_element(flo.begin(),flo.end());
  vector<double>::iterator ihi = max_element(fhi.begin(),fhi.end());
  if( *ilo > *ihi ) {
    int k = ilo - flo.begin();
    double z = division_[k];
    cut_ = SprUtils::lowerBound(z);
  }
  else {
    int k = ihi - fhi.begin();
    double z = division_[k];
    cut_ = SprUtils::upperBound(z);
  }
  
  // verbose
  if( verbose > 2 ) {
    cout << "Setting cut on variable " << d_ 
	 << " at " << cut_[0].first << " " << cut_[0].second << endl;
  }
  if( verbose > 3 ) {
    cout << "Low and high: " << *ilo << " " << *ihi << endl;
  }
  
  // exit
  return true;
}


void SprBinarySplit::print(std::ostream& os) const
{
  os << "Trained BinarySplit " << SprVersion << endl;
  os << "Dimension: " << d_ << endl;
  os << "Cut: " << cut_.size() << endl;
  for( int i=0;i<cut_.size();i++ ) {
    char s [200];
    sprintf(s,"%10g %10g",cut_[i].first,cut_[i].second);
    os << s << endl;
  }
}


bool SprBinarySplit::sort()
{
  // check
  if( nSorted_ == data_->size() )
    return true;

  // dump points into vectors
  vector<pair<double,int> > r0, r1;
  vector<double> r(data_->size());
  for( int j=0;j<data_->size();j++ ) {
    const SprPoint* p = (*data_)[j];
    r[j] = (p->x_)[d_];
    if(      p->class_ == cls0_ )
      r0.push_back(pair<double,int>((p->x_)[d_],j));
    else if( p->class_ == cls1_ )
      r1.push_back(pair<double,int>((p->x_)[d_],j));
  }
    
  // there have to be events from both categories
  if( r0.empty() || r1.empty() ) {
    cerr << "One of the categories is empty in the original data." << endl;
    return false;
  }
    
  // sort
  stable_sort(r.begin(),r.end(),less<double>());
  stable_sort(r0.begin(),r0.end(),SBSCmpPairFirst());
  stable_sort(r1.begin(),r1.end(),SBSCmpPairFirst());
    
  // fill out sorted indices
  sorted0_.clear();
  sorted1_.clear();
  sorted0_.resize(r0.size());
  sorted1_.resize(r1.size());
  for( int i=0;i<r0.size();i++ )
    sorted0_[i] = r0[i].second;
  for( int i=0;i<r1.size();i++ )
    sorted1_[i] = r1[i].second;

  // fill out divisions
  division_.clear();
  division_.push_back(SprUtils::min());
  double xprev = r[0];
  for( int k=1;k<r.size();k++ ) {
    double xcurr = r[k];
    if( (xcurr-xprev) > SprUtils::eps() ) {
      division_.push_back(0.5*(xcurr+xprev));
      xprev = xcurr;
    }
  }
  division_.push_back(SprUtils::max());

  // reset nSorted
  nSorted_ = sorted0_.size() + sorted1_.size();

  // exit
  return true;
}
