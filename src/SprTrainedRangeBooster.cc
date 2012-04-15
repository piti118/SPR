//$Id: SprTrainedRangeBooster.cc,v 1.1 2008-01-03 20:51:59 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprTrainedRangeBooster.hh"

using namespace std;


SprTrainedRangeBooster::SprTrainedRangeBooster(
  const std::vector<std::pair<const SprAbsTrainedClassifier*,bool> >& trained) 
  : 
  SprTrainedBagger(trained,false)
{}


SprTrainedRangeBooster::SprTrainedRangeBooster(
  const SprTrainedRangeBooster& other)
  :
  SprTrainedBagger(other)
{}


double SprTrainedRangeBooster::response(const std::vector<double>& v) const
{
  // determine the number of classifiers to be used
  unsigned nClassifiers = 0;
  if( nUsedClassifiers_ > 0 )
    nClassifiers = ( nUsedClassifiers_<trained_.size() ? 
		     nUsedClassifiers_ : trained_.size() );
  else
    nClassifiers = trained_.size();

  // average of acceptance rate by all classifiers but the last one
  double r = 0;
  for( int i=0;i<nClassifiers-1;i++ )
    r += ( trained_[i].first->accept(v) ? 1 : 0 );
  if( nClassifiers > 1 ) r /= (nClassifiers-1);
  if( r < 0.5 ) return 0;

  // exit
  return trained_[nClassifiers-1].first->response(v);
}


