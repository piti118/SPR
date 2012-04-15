//$Id: SprAbsClassifier.cc,v 1.3 2007-02-05 21:49:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprAbsClassifier.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"

#include <fstream>
#include <string>
#include <vector>

using namespace std;


bool SprAbsClassifier::store(const char* filename) const
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
