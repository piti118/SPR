//$Id: SprClass.cc,v 1.7 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprClass.hh"

#include <stdlib.h>
#include <ios>
#include <algorithm>
#include <sstream>

using namespace std;


bool SprClass::operator==(int cls) const 
{
  if( negate_ ) {
    for( int i=0;i<classes_.size();i++ )
      if( cls == classes_[i] ) return false;
    return true;
  }
  else {
    for( int i=0;i<classes_.size();i++ )
      if( cls == classes_[i] ) return true;
    return false;
  }
  return false;
}


bool SprClass::operator==(const SprClass& other) const 
{
  if( negate_ != other.negate_ ) return false;
  if( classes_.size() != other.classes_.size() ) return false;
  for( int i=0;i<classes_.size();i++ ) {
    if( find(other.classes_.begin(),other.classes_.end(),classes_[i]) 
	== other.classes_.end() ) return false;
  }
  for( int i=0;i<other.classes_.size();i++ ) {
    if( find(classes_.begin(),classes_.end(),other.classes_[i]) 
	== classes_.end() ) return false;
  }
  return true;
}


bool SprClass::checkClasses() const
{
  for( vector<int>::const_iterator iter=classes_.begin();
       iter!=classes_.end();iter++ ) {
    if( find(iter+1,classes_.end(),*iter) != classes_.end() ) {
      cerr << "Class " << *iter << " has been entered twice." << endl;
      return false;
    }
  }
  return true;
}


int SprClass::overlap(const SprClass& other) const
{
  if( negate_ || other.negate_ ) return -1;
  for( int i=0;i<classes_.size();i++ ) {
    if( find(other.classes_.begin(),other.classes_.end(),classes_[i]) 
	!= other.classes_.end() ) return 1;
  }
  for( int i=0;i<other.classes_.size();i++ ) {
    if( find(classes_.begin(),classes_.end(),other.classes_[i]) 
	!= classes_.end() ) return 1;
  }
  return 0;
}


std::string SprClass::toString() const
{
  ostringstream os;
  os << *this;
  return os.str();
}


std::istream& operator>>(std::istream& is, SprClass& cls) {
  // get the whole string
  string strcls;
  is >> strcls;
  if( strcls.empty() ) {
    is.clear(ios_base::failbit);
    return is;
  }

  // get list of classes and negation
  string::size_type left_parenthesis = strcls.find_first_of('(');
  if( left_parenthesis<1 || left_parenthesis==string::npos ) {
    is.clear(ios_base::failbit);
    return is;
  }
  string::size_type right_parenthesis = strcls.find_first_of(')');
  if( right_parenthesis<3 || right_parenthesis==string::npos ) {
    is.clear(ios_base::failbit);
    return is;
  }
  string list_of_classes = strcls.substr( 0, left_parenthesis );
  string negation = strcls.substr( right_parenthesis-1, 
					right_parenthesis );

  // process negation
  int negate = atoi(negation.c_str());
  if( negate!=-1 && negate!=1 ) {
    is.clear(ios_base::failbit);
    return is;
  }

  // process classes
  vector<int> classes;
  string::size_type pos_comma = list_of_classes.find_first_of(',');
  while( pos_comma != string::npos ) {
    int intcls = atoi(list_of_classes.substr(0,pos_comma).c_str());
    classes.push_back(intcls);    
    list_of_classes.erase(0,pos_comma+1);
    pos_comma = list_of_classes.find_first_of(',');
  }

  // make SprClass
  cls.classes_ = classes;
  cls.negate_ = negate;
  if( !cls.checkClasses() ) {
    is.clear(ios_base::failbit);
    return is;
  }

  // exit
  return is;
}
