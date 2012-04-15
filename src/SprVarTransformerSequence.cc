// $Id: SprVarTransformerSequence.cc,v 1.2 2008-02-08 23:22:09 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprVarTransformerSequence.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprTransformerFilter.hh"

#include <cassert>

using namespace std;


SprVarTransformerSequence::~SprVarTransformerSequence()
{
  for( int i=0;i<transformers_.size();i++ ) {
    if( transformers_[i].second )
      delete transformers_[i].first;
  }
}


SprVarTransformerSequence::SprVarTransformerSequence()
  :
  SprAbsVarTransformer(),
  transformers_()
{}


SprVarTransformerSequence::SprVarTransformerSequence(
   const std::vector<std::pair<SprAbsVarTransformer*,bool> >& transformers)
  :
  SprAbsVarTransformer(),
  transformers_(transformers)
{}
  

SprVarTransformerSequence::SprVarTransformerSequence(const 
						     SprVarTransformerSequence&
						     other)
  :
  SprAbsVarTransformer(other),
  transformers_()
{
  for( int i=0;i<other.transformers_.size();i++ ) {
    transformers_.push_back(pair<SprAbsVarTransformer*,
			    bool>(other.transformers_[i].first->clone(),true));
  }
}


bool SprVarTransformerSequence::allVarsIndependent() const
{
  for( int i=0;i<transformers_.size();i++ )
    if( !transformers_[i].first->allVarsIndependent() ) return false;
  return true;
}


bool SprVarTransformerSequence::reduceVars(const std::vector<std::string>& 
					   vars)
{
  // sanity check
  if( !this->allVarsIndependent() || transformers_.empty() ) return true;

  // reduce vars
  vector<string> allowedVars = vars;
  for( int i=0;i<transformers_.size();i++ ) {
    assert( transformers_[i].first != 0 );
    if( !transformers_[i].first->reduceVars(allowedVars) ) {
      cerr << "Unable to reduce variable list for transformer " 
	   << transformers_[i].first->name().c_str() << endl;
      return false;
    }
    transformers_[i].first->newVars(allowedVars);
  }  

  // reset vars
  transformers_[0].first->oldVars(oldVars_);
  transformers_[transformers_.size()-1].first->newVars(newVars_);

  // exit
  return true;
}


bool SprVarTransformerSequence::train(const SprAbsFilter* data, int verbose)
{
  // init current copy of data
  SprTransformerFilter* currentData = new SprTransformerFilter(data);

  // loop thru transformers
  for( int i=0;i<transformers_.size();i++ ) {
    SprAbsVarTransformer* transformer = transformers_[i].first;
    assert( transformer != 0 );

    // train
    if( !transformer->train(currentData) ) {
      cerr << "Unable to train transformer " 
	   << transformer->name().c_str() << endl;
      return false;
    }

    // transform data
    bool replaceOriginalData = true;
    if( !currentData->transform(transformer,replaceOriginalData) ) {
      cerr << "Unable to transform data with transformer " 
	   << transformer->name().c_str() << endl;
      return false;
    }

    // update current copy
    SprTransformerFilter* temp = new SprTransformerFilter(currentData);
    delete currentData;
    currentData = temp;
  }

  // set variables
  return this->initVars();
}


void SprVarTransformerSequence::transform(const std::vector<double>& in,
					  std::vector<double>& out) const
{
  out.clear();
  vector<double> temp = in;
  for( int i=0;i<transformers_.size();i++ ) {
    transformers_[i].first->transform(temp,out);
    temp = out;
  }
}


void SprVarTransformerSequence::inverse(const std::vector<double>& in,
					std::vector<double>& out) const
{
  out.clear();
  vector<double> temp = in;
  for( int i=transformers_.size()-1;i>=0;i-- ) {
    transformers_[i].first->inverse(temp,out);
    temp = out;
  }
}


bool SprVarTransformerSequence::ready() const
{
  if( oldVars_.empty() || newVars_.empty() ) return false;
  for( int i=0;i<transformers_.size();i++ ) {
    if( !transformers_[i].first->ready() ) return false;
  }
  return true;
}


void SprVarTransformerSequence::print(std::ostream& os) const
{
  // save name
  os << "VarTransformer: " << this->name().c_str() 
     << " " << SprVersion.c_str() << endl;

  // 2nd line contains all transformer names
  os << transformers_.size() << " ";
  for( int i=0;i<transformers_.size();i++ ) {
    assert( transformers_[i].first != 0 );
    os << transformers_[i].first->name().c_str() << " ";
  }
  os << endl;
  
  // print out all transformers
  for( int i=0;i<transformers_.size();i++ )
    transformers_[i].first->printWithVars(os);
}


bool SprVarTransformerSequence::initVars()
{
  if( transformers_.empty() ) return false;
  if( transformers_[0].first == 0 ) return false;
  transformers_[0].first->oldVars(oldVars_);
  if( transformers_[transformers_.size()-1].first == 0 ) return false;
  transformers_[transformers_.size()-1].first->newVars(newVars_);
  if( oldVars_.empty() || newVars_.empty() ) return false;
  return true;
}
