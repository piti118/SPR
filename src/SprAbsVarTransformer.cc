// $Id: SprAbsVarTransformer.cc,v 1.2 2008-01-03 20:51:59 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprAbsVarTransformer.hh"

#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;


bool SprAbsVarTransformer::store(const char* filename) const
{
  // open file for output
  string fname = filename;
  ofstream os(fname.c_str());
  if( !os ) {
    cerr << "Unable to open file " << fname.c_str() << endl;
    return false;
  }

  // store into file
  this->printWithVars(os);

  // exit
  return true;
}


void SprAbsVarTransformer::printWithVars(std::ostream& os) const
{
  // print
  this->print(os);

  // store old variables
  os << "==================================================" << endl;
  os << "Old Variables:" << endl;
  for( int i=0;i<oldVars_.size();i++ ) {
    char s [200];
    sprintf(s,"%5i %40s",i,oldVars_[i].c_str());
    os << s << endl;
  }
  os << "==================================================" << endl;

  // store new variables
  os << "==================================================" << endl;
  os << "New Variables:" << endl;
  for( int i=0;i<newVars_.size();i++ ) {
    char s [200];
    sprintf(s,"%5i %40s",i,newVars_[i].c_str());
    os << s << endl;
  }
  os << "==================================================" << endl;
}
