//$Id: SprConfig.cc,v 1.1 2008-05-08 19:57:43 narsky Exp $
#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprConfig.hh"

#include <stdlib.h>
#include <fstream>
#include <cassert>

using namespace std;

const char SprConfig::commentSignal_ = '#';
const char SprConfig::keyValueSeparator_ = '=';


SprConfig::SprConfig(const std::string& filename)
  :
  config_()
{
  bool status = this->read(filename);
  assert( status );
}

SprConfig::SprConfig(const SprConfig& configObject)
  :
  config_(configObject.config_)
{}

std::string SprConfig::getStringValue(const std::string& key, 
				      const std::string defaultValue) const
{
  map<string,string>::const_iterator found = config_.find(key);
  if( found != config_.end() ) 
    return found->second;
  return defaultValue;
}

int SprConfig::getIntValue(const std::string& key, const int defaultValue) const
{
  map<string,string>::const_iterator found = config_.find(key);
  if( found != config_.end() ) 
    return atoi(found->second.c_str());
  return defaultValue;
}

double SprConfig::getDoubleValue(const std::string& key, 
				 const double defaultValue) const
{
  map<string,string>::const_iterator found = config_.find(key);
  if( found != config_.end() ) 
    return atof(found->second.c_str());
  return defaultValue;
}

bool SprConfig::getBoolValue(const std::string& key, 
			     const bool defaultValue) const
{
  map<string,string>::const_iterator found = config_.find(key);
  if( found != config_.end() ) 
    return ( found->second.compare("true")==0 ? true : false );
  return defaultValue;
}

void SprConfig::print(std::ostream& os) const
{
  for ( map<string,string>::const_iterator 
	  it=config_.begin(); it != config_.end(); it++ ){
    os << it->first << " => " << it->second << endl;
  }
}

void SprConfig::trim(std::string& s, const std::string drop )
{
  s.erase(s.find_last_not_of(drop)+1);
  s.erase(0,s.find_first_not_of(drop));
}

/*
void SprConfig::ltrim(std::string& s, const std::string drop )
{
  s.erase(0,s.find_first_not_of(drop));
}
*/

bool SprConfig::keyExists(const std::string& key) const
{
  return ( config_.find(key) != config_.end() );
}


bool SprConfig::read(const std::string& filename)
{
  // open file
  ifstream fin(filename.c_str(),ios::in);
  if( !fin ) {
    cerr << "Error Opening File: " << filename.c_str() << endl;
    return false;
  }

  // read  
  unsigned nLine = 0;
  string line;
  while( getline(fin,line) ){
    nLine++;

    //first strip all comment after #
    string::size_type commentcutoff 
      = line.find_first_of(SprConfig::commentSignal_);
    line = line.substr(0,commentcutoff);

    // trim and see if there is anything there
    trim(line);
    if( line.empty() ) continue;

    // find separator
    string::size_type keycutoff 
      = line.find_first_of(SprConfig::keyValueSeparator_);
    if( keycutoff == string::npos ) {
      cerr << "No key separator found on line " << nLine << endl;
      return false;
    }

    // key
    string key = line.substr(0,keycutoff);
    trim(key);
    if( key.empty() ) {
      cerr << "Key is empty on line " << nLine << endl;
      return false;
    }
      
    // value
    string value = line.substr(keycutoff+1);
    trim(value);
    if( value.empty() ) {
      cerr << "Value is empty on line " << nLine << endl;
      return false;
    }

    // store
    if( !config_.insert(pair<const string,string>(key,value)).second ) {
      cerr << "Unable to insert key and value into config map: " 
	   << key.c_str() << " " << value.c_str() << endl;
      return false;
    }
  }

  // exit
  return true;
}
