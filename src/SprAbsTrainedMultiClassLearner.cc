//$Id: SprAbsTrainedMultiClassLearner.cc,v 1.2 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"

#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <functional>

using namespace std;


struct SATMCLCmpPairSecond
  : public binary_function<pair<const int,double>,
			   pair<const int,double>,bool> {
  bool operator()(const pair<const int,double>& l, 
		  const pair<const int,double>& r) const {
    return (l.second < r.second);
  }
};


SprAbsTrainedMultiClassLearner::SprAbsTrainedMultiClassLearner()
  :
  vars_(),
  mapper_(),
  validRange_(),
  defaultMissing_()
{}


SprAbsTrainedMultiClassLearner::SprAbsTrainedMultiClassLearner(
				       const std::vector<int>& classes)
  :
  vars_(),
  mapper_(classes),
  validRange_(),
  defaultMissing_()
{}


SprAbsTrainedMultiClassLearner::SprAbsTrainedMultiClassLearner(
			       const SprAbsTrainedMultiClassLearner& other)
  :
  vars_(other.vars_),
  mapper_(other.mapper_),
  validRange_(other.validRange_),
  defaultMissing_(other.defaultMissing_)
{}


int SprAbsTrainedMultiClassLearner::response(const std::vector<double>& input,
					     std::map<int,double>& output) 
  const
{
  if( validRange_.empty() || defaultMissing_.empty() )
    return this->response_one(input,output);
  return this->response_many(input,output);
}


int SprAbsTrainedMultiClassLearner::response_many(
			 const std::vector<double>& input,
			 std::map<int,double>& output) const
{
  // init vectors
  vector<MissingType> replaced = defaultMissing_;

  // loop over components to find those outside allowed range
  bool foundInvalid = false;
  for( int d=0;d<input.size();d++ ) {
      bool valid = false;
      double r = input[d];
      for( int j=0;j<validRange_.size();j++ ) {
	if( r>validRange_[j].first && r<validRange_[j].second ) {
	  valid = true;
	  break;
	}
      }
      if( !valid ) {
	foundInvalid = true;
	continue;
      }
      
      // copy valid values
      for( int ic=0;ic<replaced.size();ic++ )
	replaced[ic].second[d] = r;
  }

  // clear responses
  output.clear();
  map<int,double> replaced_output;

  // compute responses with replaced values
  if( foundInvalid ) {
    for( int ic=0;ic<replaced.size();ic++ ) {
      int cls = this->response_one(replaced[ic].second,replaced_output);
      output.insert(pair<const int,double>(replaced[ic].first,
					   replaced_output[cls]));
    }
    map<int,double>::const_iterator iter 
      = min_element(output.begin(),output.end(),SATMCLCmpPairSecond());
    return iter->first;
  }

  // if no invalid values found, fill output with the same value
  int cls = this->response_one(replaced[0].second,replaced_output);
  double value = replaced_output[cls];
  for( int ic=0;ic<replaced.size();ic++ )
    output.insert(pair<const int,double>(replaced[ic].first,value));
  return cls;
}


bool SprAbsTrainedMultiClassLearner::store(const char* filename) const
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
  os << "==================================================" << endl;
  os << "Dimensions:" << endl;
  for( int i=0;i<vars_.size();i++ ) {
    char s [200];
    sprintf(s,"%5i %40s",i,vars_[i].c_str());
    os << s << endl;
  }
  os << "==================================================" << endl;

  // exit
  return true;
}


bool SprAbsTrainedMultiClassLearner::setDefaultMissing(
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


void SprAbsTrainedMultiClassLearner::outputNames(
				    const char* prefix,
				    std::vector<std::string>& names) const
{
  // init
  string sprefix = prefix;
  names.clear();

  // process missing values
  if( defaultMissing_.empty() ) {
    this->outputNamesWithoutMissing(sprefix.c_str(),names);
  }
  else {
    for( int i=0;i<defaultMissing_.size();i++ ) {
      string axis = sprefix;
      char s [200];
      sprintf(s,"%i",defaultMissing_[i].first);
      axis += s;
      names.push_back(axis);
    }
  }
}


void SprAbsTrainedMultiClassLearner::outputNamesWithoutMissing(
		       		      const char* prefix,
				      std::vector<std::string>& names) const
{
  string sprefix = prefix;
  names.clear();
  for( int i=0;i<mapper_.size();i++ ) {
    string axis = sprefix;
    char s [200];
    sprintf(s,"%i",mapper_[i]);
    axis += s;
    names.push_back(axis);
  }
}
