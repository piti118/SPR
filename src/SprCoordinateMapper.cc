//$Id: SprCoordinateMapper.cc,v 1.3 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprCoordinateMapper.hh"
#include "StatPatternRecognition/SprPoint.hh"

#include <cassert>
#include <algorithm>
#include <iostream>

using namespace std;


SprCoordinateMapper* SprCoordinateMapper::createMapper(
				      const std::vector<std::string>& from,
				      const std::vector<std::string>& to)
{
  vector<unsigned> mapper(from.size());
  for( int i=0;i<from.size();i++ ) {
    vector<string>::const_iterator found = find(to.begin(),to.end(),from[i]);
    if( found == to.end() ) {
      cout << "Warning: CoordinateMapper cannot find variable " 
	   << from[i].c_str() << " among data variables." << endl;
      return 0;
    }
    int d = found - to.begin();
    mapper[i] = d;
  }
  return SprCoordinateMapper::createMapper(mapper);
}


const SprPoint* SprCoordinateMapper::output(const SprPoint* input) 
{
  // sanity check
  if( mapper_.empty() ) return input;
  
  // make new point and copy index+class
  SprPoint* p = new SprPoint;
  p->index_ = input->index_;
  p->class_ = input->class_;
  
  // copy vector elements
  this->map(input->x_,p->x_);
  
  // add to the cleanup list
  toDelete_.push_back(p);
  
  // exit
  return p;
}


void SprCoordinateMapper::map(const std::vector<double>& in,
			      std::vector<double>& out) const
{
  out.clear();
  for( int i=0;i<mapper_.size();i++ ) {
    unsigned d = mapper_[i];
    assert( d < in.size() );
    out.push_back(in[d]);
  }
}


void SprCoordinateMapper::clear() 
{
  for( int i=0;i<toDelete_.size();i++ ) delete toDelete_[i];
  toDelete_.clear();
}
