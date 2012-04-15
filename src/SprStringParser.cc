//$Id: SprStringParser.cc,v 1.3 2006-11-26 02:04:31 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprStringParser.hh"

#include <stdlib.h>

using namespace std;

void SprStringParser::parseToStrings(const char* input, 
			  std::vector<std::vector<std::string> >& output)
{
  // init
  string s = input;
  output.clear();

  // split into subvectors using ":" as separator
  vector<string> v;
  while( !s.empty() ) {
    string::size_type pos = s.find_first_of(':');
    string s1 = s.substr(0,pos);
    if( pos != string::npos ) 
      s = s.substr(pos+1);
    else
      s.clear();
    s1.erase( s1.find_last_not_of(' ') + 1 );
    s1.erase( 0, s1.find_first_not_of(' ') );
    v.push_back(s1);
  }

  // split each element of subvectors using "," as separator
  output.resize(v.size());
  for( int i=0;i<v.size();i++ ) {
    s = v[i];
    while( !s.empty() ) {
      string::size_type pos = s.find_first_of(',');
      string s1 = s.substr(0,pos);
      if( pos != string::npos ) 
	s = s.substr(pos+1);
      else
	s.clear();
      s1.erase( s1.find_last_not_of(' ') + 1 );
      s1.erase( 0, s1.find_first_not_of(' ') );
      output[i].push_back(s1);
    }
  }
}


void SprStringParser::parseToInts(const char* input, 
				  std::vector<std::vector<int> >& output)
{
  output.clear();
  vector<vector<string> > str_output;
  SprStringParser::parseToStrings(input,str_output);
  output.resize(str_output.size());
  for( int i=0;i<str_output.size();i++ ) {
    output[i].resize(str_output[i].size());
    for( int j=0;j<str_output[i].size();j++ ) {
      output[i][j] = atoi(str_output[i][j].c_str());
    }
  }
}


void SprStringParser::parseToDoubles(const char* input, 
				     std::vector<std::vector<double> >& output)
{
  output.clear();
  vector<vector<string> > str_output;
  SprStringParser::parseToStrings(input,str_output);
  output.resize(str_output.size());
  for( int i=0;i<str_output.size();i++ ) {
    output[i].resize(str_output[i].size());
    for( int j=0;j<str_output[i].size();j++ ) {
      output[i][j] = atof(str_output[i][j].c_str());
    }
  }
}
