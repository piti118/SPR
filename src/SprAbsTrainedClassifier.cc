//$Id: SprAbsTrainedClassifier.cc,v 1.7 2007-11-07 00:56:14 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprUtils.hh"

#include <stdio.h>
#include <utility>
#include <fstream>

using namespace std;


bool SprAbsTrainedClassifier::accept(const std::vector<double>& v, 
				     double& response) const
{
  response = this->response(v);
  if( cut_.empty() ) return true;
  bool passed = false;
  for( int i=0;i<cut_.size();i++ ) {
    const pair<double,double>& lims = cut_[i];
    if( response>lims.first && response<lims.second ) {
      passed = true;
      break;
    }
  }
  return passed;
}


bool SprAbsTrainedClassifier::store(const char* filename) const
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

std::vector<std::string>* SprAbsTrainedClassifier::getVars() const{
    std::vector<std::string>* toReturn = new std::vector<std::string>(vars_);
    return toReturn;
}


bool SprAbsTrainedClassifier::storeCode(const char* filename) const
{
  // open file for output
  string fname = filename;
  ofstream os(fname.c_str());
  if( !os ) {
    cerr << "Unable to open file " << fname.c_str() << endl;
    return false;
  }

  // store
  return this->generateCode(os);
}
