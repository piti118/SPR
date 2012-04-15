//$Id: SprMultiClassPlotter.cc,v 1.3 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprMultiClassPlotter.hh"
#include "StatPatternRecognition/SprLoss.hh"
#include "StatPatternRecognition/SprAverageLoss.hh"

#include <cassert>
#include <iostream>
#include <utility>
#include <algorithm>

using namespace std;


double SprMultiClassPlotter::multiClassTable(
			      const std::vector<int>& classes,
			      SprAverageLoss* loss,
			      SprClassificationTable& classificationTable,
			      std::map<int,double>& weightInClass,
			      bool normalizePerClass) const
{
  // init
  assert( loss != 0 );
  loss->reset();
  classificationTable.clear();
  weightInClass.clear();
  
  // sanity check
  if( responses_.empty() ) {
    cerr << "No responses have been computed." << endl;
    return 0;
  }
  const int nClasses = classes.size();
  assert( nClasses > 0 );

  // loop over responses
  for( int i=0;i<responses_.size();i++ ) {

    // true class
    int cls = responses_[i].cls;
    if( !classes.empty() && 
	find(classes.begin(),classes.end(),cls)==classes.end() ) continue;

    // find position of the class in the table
    int assigned = responses_[i].assigned;
    vector<int>::const_iterator foundClass = 
      find(classes.begin(),classes.end(),assigned);
    if( foundClass == classes.end() ) continue;
    int iassigned = foundClass - classes.begin();

    // update weights
    double w = responses_[i].weight;
    map<int,double>::iterator weighted = weightInClass.find(cls);
    if( weighted == weightInClass.end() ) {
      pair<map<int,double>::iterator,bool> inserted
	= weightInClass.insert(pair<const int,double>(cls,w));
      assert( inserted.second );
    }
    else
      weighted->second += w;

    // update overall loss
    loss->update(cls,double(assigned),w);

    // indvidual loss mapped to classes
    map<int,vector<double> >::iterator tabulated 
      = classificationTable.find(cls);
    if( tabulated == classificationTable.end() ) {
      pair<map<int,vector<double> >::iterator,bool> inserted
	= classificationTable.insert(pair<const int,
			    vector<double> >(cls,vector<double>(nClasses,0)));
      assert( inserted.second );
      tabulated = inserted.first;
    }
    (tabulated->second)[iassigned] += w;
  }

  // normalize classification table for each class row
  if( normalizePerClass ) {
    for( map<int,vector<double> >::iterator i=classificationTable.begin();
	 i!=classificationTable.end();i++ ) {
      int cls = i->first;
      double wtot = weightInClass[cls];
      assert( wtot > 0 );
      vector<double>& w = i->second;
      for( int j=0;j<w.size();j++ ) w[j] /= wtot;
    }
  }

  // exit
  return loss->value();
}
