// $Id: SprVarSelectorFilter.cc,v 1.1 2008-01-03 20:51:59 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprVarSelectorFilter.hh"
#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprData.hh"
#include "StatPatternRecognition/SprCoordinateMapper.hh"

#include <iostream>
#include <algorithm>
#include <memory>

using namespace std;


bool SprVarSelectorFilter::chooseVars(const std::set<std::string>& vars)
{
  // convert set into a vector
  vars_.clear();
  for( set<string>::const_iterator i=vars.begin();i!=vars.end();i++ )
    vars_.push_back(*i);
  return true;
}


bool SprVarSelectorFilter::filter()
{
  // sanity check
  if( vars_.empty() ) {
    cerr << "Variable list for VarSelectorFilter is empty." << endl;
    return false;
  }

  // make sure all requested variables are found in the original data
  vector<string> dataVars;
  data_->vars(dataVars);
  auto_ptr<SprCoordinateMapper> 
    mapper(SprCoordinateMapper::createMapper(vars_,dataVars));
  if( mapper.get() == 0 ) {
    cerr << "VarSelectorFilter cannot make a mapper object." << endl;
    return false;
  }

  // make a new empty copy
  bool ownPoints = true;
  SprData* copy = new SprData("reduced_data",vars_,ownPoints);

  // find index range
  int istart = ( imin_<0 ? 0 : imin_ );
  int iend = ( (imax_>data_->size()||imax_==0) ? data_->size() : imax_ );

  // loop through points and accept
  copyWeights_.clear();
  for( int i=istart;i<iend;i++ ) {
    SprPoint* p = (*data_)[i];
    if( this->category(p) ) {
      vector<double> newX;
      mapper->map(p->x_,newX);
      SprPoint* newP = new SprPoint(p->index_,p->class_,newX);
      copy->uncheckedInsert(newP);
      copyWeights_.push_back(dataWeights_[i]);
    }
  }

  // save copy
  if( ownCopy_ ) delete copy_;
  copy_ = copy;
  ownCopy_ = true;

  // exit
  return true;
}
