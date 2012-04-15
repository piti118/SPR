// $Id: SprTrainedDecisionTree.cc,v 1.5 2007-05-14 18:08:08 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprTrainedDecisionTree.hh"

#include <stdio.h>
#include <cassert>

using namespace std;


double SprTrainedDecisionTree::response(const std::vector<double>& v) const
{
  // go through signal nodes
  for( int i=0;i<nodes1_.size();i++ ) {
    bool accepted = true;
    for( SprBox::const_iterator iter=nodes1_[i].begin();
	 iter!=nodes1_[i].end();iter++ ) {
      unsigned d = iter->first;
      assert( d < v.size() );
      if( v[d]<iter->second.first || v[d]>iter->second.second ) {
	accepted = false;
	break;
      }
    }
    if( accepted ) return 1;
  }
  return 0;
}


int SprTrainedDecisionTree::nBox(const std::vector<double>& v) const
{
  // go through signal nodes
  for( int i=0;i<nodes1_.size();i++ ) {
    bool accepted = true;
    for( SprBox::const_iterator iter=nodes1_[i].begin();
	 iter!=nodes1_[i].end();iter++ ) {
      unsigned d = iter->first;
      assert( d < v.size() );
      if( v[d]<iter->second.first || v[d]>iter->second.second ) {
	accepted = false;
	break;
      }
    }
    if( accepted ) return i;
  }
  return -1;
}


void SprTrainedDecisionTree::print(std::ostream& os) const
{
  os << "Trained DecisionTree " << SprVersion << endl;
  os << "Nodes: " << nodes1_.size() << " nodes." << endl;
  for( int i=0;i<nodes1_.size();i++ ) {
    const SprBox& limits = nodes1_[i];
    int size = limits.size();
    os << "Node " << i << " Size " << size << endl;
    for( SprBox::const_iterator iter=limits.begin();
	 iter!=limits.end();iter++ ) {
      char s [200];
      sprintf(s,"Dimension %4i    Limits %15g %15g",
	      iter->first,iter->second.first,iter->second.second);
      os << s << endl;
    }
  }
}
