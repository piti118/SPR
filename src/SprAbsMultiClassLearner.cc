//$Id: SprAbsMultiClassLearner.cc,v 1.3 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprAbsMultiClassLearner.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprReplaceMissing.hh"
#include "StatPatternRecognition/SprTransformerFilter.hh"

#include <stdio.h>
#include <cassert>
#include <fstream>

using namespace std;


SprAbsMultiClassLearner::~SprAbsMultiClassLearner()
{
  if( ownData_ ) {
    delete data_;
    ownData_ = false;
  }
}


SprAbsMultiClassLearner::SprAbsMultiClassLearner(
					SprAbsFilter* data,
					const std::vector<int>& classes)
  :
  data_(data),
  mapper_(classes),
  validRange_(),
  defaultMissing_(),
  ownData_(false)
{
  assert( data_ != 0 );
  assert( !mapper_.empty() );
  bool status = this->replaceMissing();
  assert( status );
}


bool SprAbsMultiClassLearner::checkClasses() const
{
  // sanity check
  if( mapper_.size() < 2 ) {
    cerr << "Less than 2 classes are specified." << endl;
    return false;
  }
  
  // make sure all classes are distinct
  for( int i=0;i<mapper_.size();i++ ) {
    for( int j=i+1;j<mapper_.size();j++ ) {
      if( mapper_[i] == mapper_[j] ) {
	cerr << "Elements " << i << " and " << j 
	     << " of the input vector of classes are equal." << endl;
	return false;
      }
    }
  }

  // exit
  return true;
}


bool SprAbsMultiClassLearner::store(const char* filename) const
{
  // open file for output
  string fname = filename;
  ofstream os(fname.c_str());
  if( !os ) {
    cerr << "Unable to open file " << fname.c_str() << endl;
    return false;
  }

  // store into file
  this->print(os);

  // store missing values
  os << "==================================================" << endl;
  os << "Valid Range: " << validRange_.size();
  for( int i=0;i<validRange_.size();i++ )
    os << " " << validRange_[i].first << " " << validRange_[i].second;
  os << endl;
  os << "Default values: " << defaultMissing_.size() << endl;
  for( int i=0;i<defaultMissing_.size();i++ ) {
    os << "Class: " << defaultMissing_[i].first
       << " Size: " << defaultMissing_[i].second.size()
       << " Values:";
    for( int j=0;j<defaultMissing_[i].second.size();j++ )
      os << " " << defaultMissing_[i].second[j];
    os << endl;
  }

  // store variables
  vector<string> vars;
  data_->vars(vars);
  assert( vars.size() == data_->dim() );
  os << "==================================================" << endl;
  os << "Dimensions:" << endl;
  for( int i=0;i<vars.size();i++ ) {
    char s [200];
    sprintf(s,"%5i %40s",i,vars[i].c_str());
    os << s << endl;
  }
  os << "==================================================" << endl;

  // exit
  return true;
}


bool SprAbsMultiClassLearner::setDefaultMissing(
				const SprCut& validRange,
				const std::vector<MissingType>& defaultMissing)
{
  // sanity check
  if( validRange.empty() || defaultMissing.empty() ) return true;

  // set the cut
  validRange_ = validRange;

  // check size of missing values
  int dim = defaultMissing[0].second.size();
  for( int ic=1;ic<defaultMissing.size();ic++ ) 
    assert( defaultMissing[ic].second.size() == dim );

  // set missing values only for included classes
  defaultMissing_.clear();
  for( int ic=0;ic<mapper_.size();ic++ ) {
    int found_jc = -1;
    for( int jc=0;jc<defaultMissing.size();jc++ ) {
      if( mapper_[ic] == defaultMissing[jc].first ) {
	found_jc = jc;
	break;
      }
    }
    if( found_jc < 0 ) continue;
    defaultMissing_.push_back(MissingType(mapper_[ic],
					  defaultMissing[found_jc].second));
  }
    
  // exit
  return true;
}


bool SprAbsMultiClassLearner::replaceMissing()
{
  // sanity check
  if( validRange_.empty() || defaultMissing_.empty() ) return true;

  // make new data
  if( ownData_ ) {
    cerr << "Data for multiclass learner has missing values replaced already." 
	 << " Cannot continue." << endl;
    return false;
  }
  SprTransformerFilter* data = new SprTransformerFilter(data_);

  // replace values
  bool replaceOriginalData = false;
  bool classBlind = false;
  SprReplaceMissing trans(SprReplaceMissing::Median,validRange_,classBlind);
  int nClass = defaultMissing_.size();
  vector<SprReplaceMissing::ClassAndDefaultValues> replacement(nClass);
  for( int ic=0;ic<nClass;ic++ ) {
    replacement[ic].first = defaultMissing_[ic].first;
    replacement[ic].second = defaultMissing_[ic].second;
  }
  trans.setReplacement(replacement);
  if( !data->transform(&trans,replaceOriginalData) ) {
    cerr << "Unable to replace missing values for multiclass learner." << endl;
    delete data;
    return false;
  }

  // reassign data
  data_ = data;
  ownData_ = true;

  // exit
  return true;
}
