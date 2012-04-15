//$Id: SprRangeBooster.cc,v 1.2 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprRangeBooster.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprEmptyFilter.hh"
#include "StatPatternRecognition/SprData.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <utility>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <memory>

using namespace std;


struct SRBCmpPairIVWvalue
  : public binary_function<SprRangeBooster::IndexValueWeight,
			   SprRangeBooster::IndexValueWeight,bool> {
  bool operator()(const SprRangeBooster::IndexValueWeight& l, 
		  const SprRangeBooster::IndexValueWeight& r) const {
    return (l.value < r.value);
  }
};


struct SRBCmpPairIVWDvalue
  : public binary_function<SprRangeBooster::IndexValueWeight,
			   double,bool> {
  bool operator()(const SprRangeBooster::IndexValueWeight& l, double r)
    const {
    return (l.value < r);
  }
};


SprRangeBooster::SprRangeBooster(SprAbsFilter* data)
  :
  SprBagger(data), 
  signalFraction_(0.5),
  epsilon_(0.5),
  threshold_(0.1),
  notified_(false)
{}

SprRangeBooster::SprRangeBooster(SprAbsFilter* data, 
				 unsigned cycles, 
				 double signalFraction,
				 double epsilon,
				 double threshold,
				 bool discrete)
  :
  SprBagger(data,cycles,discrete), 
  signalFraction_(signalFraction),
  epsilon_(epsilon),
  threshold_(threshold),
  notified_(false)
{
  assert( signalFraction_ > 0 );
  assert( epsilon_ > 0 );
  cout << "RangeBooster initialized with signalFraction=" 
       << signalFraction_ << " epsilon=" << epsilon_
       << " threshold=" << threshold_ << endl;
}


bool SprRangeBooster::reset()
{
  notified_ = false;
  return SprBagger::reset();
}


bool SprRangeBooster::train(int verbose)
{
  // sanity check
  if( cycles_==0 || trainable_.empty() ) {
    cout << "RangeBooster will exit without training." << endl;
    return this->prepareExit(true);
  }

  // copy input data into temp containers
  SprEmptyFilter origData(data_);// index0, no reweighting
  SprEmptyFilter trainData(data_);// index1, reweighted

  // update responses
  map<const SprPoint*,IndexValueWeight> responses;
  for( int i=0;i<trainData.size();i++ ) {
    const SprPoint* p = trainData[i];
    assert( p == origData[i] );
    map<const SprPoint*,IndexValueWeight>::iterator inserted
      = responses.insert(pair<const SprPoint* const,
			 IndexValueWeight>(p,IndexValueWeight())).first;
    double resp = 0;
    for( int j=0;j<trained_.size();j++ )
      resp += trained_[j].first->response(p);
    double w = 0;
    if( !trained_.empty() ) {
      w = trained_.size();
      resp /= w;
    }
    inserted->second.index0 = i;
    inserted->second.index1 = i;
    inserted->second.value  = resp;
    inserted->second.weight = w;
  }

  // after all betas are filled, do an overall validation
  if( !this->initValBeta() )
    return this->prepareExit(false);

  // if resume training, generate a seed from time of day
  int seed = 0;
  if( cycles_>0 && !trained_.empty() ) seed = -1;

  // compute max signal weight allowed in the region
  double acceptSignalWeight 
    = 0.5 * signalFraction_ * trainData.weightInClass(cls1_);
  if( verbose > 0 ) {
    cout << "Max allowed signal fraction for RangeBooster set to " 
	 << acceptSignalWeight << endl;
  }

  // loop through trainable
  unsigned nCycle = 0;
  unsigned nFailed = 0;
  while( nCycle < cycles_ ) {
    for( int i=0;i<trainable_.size();i++ ) {
      // check cycles
      if( nCycle++ >= cycles_ ) return this->prepareExit((this->nTrained()>0));

      // split the data in 2 equal subsets using original weights
      SprEmptyFilter data1(&trainData);

      // set original weights
      for( int j=0;j<data1.size();j++ ) {
	const SprPoint* p = data1[j];
	map<const SprPoint*,IndexValueWeight>::iterator 
	  found = responses.find(p);
	if( found == responses.end() ) {
	  cerr << "RangeBooster cannot find point " << j
	       << " in the response map." << endl;
	  return this->prepareExit(false);
	}
	data1.setW(j,origData.w(found->second.index0));
      }

      // split
      vector<double> splitWeights;
      bool randomize = true;
      SprData* splitted = data1.split(0.5,splitWeights,randomize,seed);
      if( splitted == 0 ) {
	cerr << "RangeBooster is unable to split input data." << endl;
	return this->prepareExit(false);
      }
      bool ownSplitted = true;
      SprEmptyFilter data2(splitted,splitWeights,ownSplitted);

      // get new classifier
      SprAbsClassifier* c = trainable_[i];
      if( !c->setData(&data1) ) {
	cerr << "Unable to set data for classifier " << i << endl;
	return this->prepareExit(false);
      }
      if( !c->train(verbose-1) ) {
	cerr << "RangeBooster failed to train classifier " << i 
	     << ". Continuing..."<< endl;
	if( ++nFailed >= cycles_ ) {
	  cout << "Exiting after failed to train " << nFailed 
	       << " classifiers." << endl;
	  return this->prepareExit((this->nTrained()>0));
	}
	else
	  continue;
      }

      // register new trained classifier
      SprAbsTrainedClassifier* t = c->makeTrained();
      if( t == 0 ) {
	cerr << "RangeBooster failed to train classifier " << i 
	     << ". Continuing..."<< endl;
	if( ++nFailed >= cycles_ ) {
	  cout << "Exiting after failed to train " << nFailed 
	       << " classifiers." << endl;
	  return this->prepareExit((this->nTrained()>0));
	}
	else
	  continue;
      }
      t->useNormalized();

      // reweight events
      int status = this->reweight(t,&data2,&trainData,&origData,
				  responses,acceptSignalWeight,verbose);
      if(      status < 0 ) {
	cerr << "RangeBooster cannot reweight events." << endl;
	return this->prepareExit(false);
      }

      // add this classifier to the trained list
      trained_.push_back(pair<const SprAbsTrainedClassifier*,bool>(t,true));

      // message
      if( verbose>1 || (nCycle%100)==0 ) {
	cout << "Done cycle " << nCycle << endl;
      }

      // validation
      if( !this->updateValBeta(t,nCycle) )
	return this->prepareExit(false);
    }
  }

  // normal exit
  return this->prepareExit((this->nTrained()>0));
}


int SprRangeBooster::reweight(
		 SprAbsTrainedClassifier* t, 
		 const SprAbsFilter* testData,
		 SprAbsFilter* trainData,
		 const SprAbsFilter* origData,
		 std::map<const SprPoint*,IndexValueWeight>& responses,
		 double acceptSignalWeight,
		 int verbose)
{
  // sanity check
  if( t==0 || testData==0 || trainData==0 || origData==0 ) {
    cerr << "Parameters missing for SprRangeBooster::reweight." << endl;
    return -1;
  }

  // update responses for testData
  // make lists of responses, for signal and background
  vector<IndexValueWeight> signal, bgrnd;
  for( int i=0;i<testData->size();i++ ) {

    // find point in the response container
    const SprPoint* p = (*testData)[i];
    map<const SprPoint*,IndexValueWeight>::iterator 
      found = responses.find(p);
    if( found == responses.end() ) {
      cerr << "RangeBooster cannot find point " << i 
	   << " in the response map." << endl;
      return -1;
    }

    // update response
    int index0 = found->second.index0;
    int index1 = found->second.index1;
    double r = t->response(p);
    double w = testData->w(i);

    // fill out response vector
    if(      p->class_ == cls0_ )
      bgrnd.push_back(IndexValueWeight(index0,index1,r,w));
    else if( p->class_ == cls1_ )
      signal.push_back(IndexValueWeight(index0,index1,r,w));
  }

  // sanity check
  assert( bgrnd.size() == testData->ptsInClass(cls0_) );
  assert( signal.size() == testData->ptsInClass(cls1_) );

  // check total weight
  double wsig = 0;
  for( int i=0;i<signal.size();i++ )
    wsig += signal[i].weight;
  if( wsig <= acceptSignalWeight ) {
    if( !notified_ ) {
      notified_ = true;
      cout << "Weight of available RangeBooster signal data " << wsig
	   << " is less than required gross signal weight " 
	   << acceptSignalWeight << "   Reweighting will not be applied." 
	   << endl;
    }
    return 1;
  }

  // sort the signal response list in descending order
  stable_sort(signal.begin(),signal.end(),not2(SRBCmpPairIVWvalue()));

  // find the divider in response values
  double wdiv = 0;
  double rdiv = 0;
  int sigDiv = 0;
  for( int i=0;i<signal.size();i++ ) {
    wdiv += signal[i].weight;
    if( wdiv > acceptSignalWeight ) {
      if( i == 0 ) {
	rdiv = signal[0].value;
	sigDiv = 1;
      }
      else {
	rdiv = 0.5 * (signal[i-1].value+signal[i].value);
	sigDiv = i;
      }
      break;
    }
  }
  if( sigDiv == 0 ) {
    cout << "RangeBooster cannot find signal weight divider and will exit." 
	 << endl;
    return 1;
  }

  // set cut for the current classifier
  t->setCut(SprUtils::lowerBound(rdiv));

  // partition background vector
  vector<IndexValueWeight>::const_iterator firstLessRdiv 
    = stable_partition(bgrnd.begin(),bgrnd.end(),
		       bind2nd(not2(SRBCmpPairIVWDvalue()),rdiv));
  int bgrDiv = firstLessRdiv - bgrnd.begin();

  // print out
  if( verbose > 1 )
    cout << "sigDiv=" << sigDiv << "   bgrDiv=" << bgrDiv << endl;

  // reduce weights of events that are out of this range
  vector<int> indicesToRemove;
  double factor = epsilon_;
  for( int i=sigDiv;i<signal.size();i++ ) {
    int index = signal[i].index1;
    double w = factor*trainData->w(index);
    if( w < threshold_ )
      indicesToRemove.push_back(index);
    else
      trainData->setW(index,w);
  }
  for( int i=bgrDiv;i<bgrnd.size();i++ ) {
    int index = bgrnd[i].index1;
    double w = factor*trainData->w(index);
    if( w < threshold_ )
      indicesToRemove.push_back(index);
    else
      trainData->setW(index,w);
  }

  // make data to be removed
  if( !indicesToRemove.empty() ) {

    // remove points
    stable_sort(indicesToRemove.begin(),indicesToRemove.end());
    auto_ptr<SprData> dataToRemove(trainData->emptyCopy());
    for( int i=0;i<indicesToRemove.size();i++ ) {
      int index = indicesToRemove[i];
      dataToRemove->uncheckedInsert((*trainData)[index]);
    }
    if( !trainData->fastRemove(dataToRemove.get()) ) {
      cerr << "RangeBooster cannot remove points below threshold." << endl;
      return -1;
    }

    // remap the indices
    for( int i=0;i<trainData->size();i++ ) {
      const SprPoint* p = (*trainData)[i];
      map<const SprPoint*,IndexValueWeight>::iterator found 
	= responses.find(p);
      if( found == responses.end() ) {
	cerr << "RangeBooster cannot find a point for remapping." << endl;
	return -1;
      }
      found->second.index1 = i;
    }// end of remap

    // notify
    if( verbose > 0 ) {
      cout << "RangeBooster removed " << indicesToRemove.size()
	   << " points with weights below " << threshold_ << endl;
    }

  }// end of !indicesToRemove.empty()

  // print out
  if( verbose > 0 ) {
    cout << "After reweighting:   W1=" << trainData->weightInClass(cls1_)
	 << " W0=" << trainData->weightInClass(cls0_)
	 << "    N1=" << trainData->ptsInClass(cls1_)
	 << " N0=" << trainData->ptsInClass(cls0_) << endl;
  }

  // exit
  return 0;
}


SprTrainedRangeBooster* SprRangeBooster::makeTrained()
{
  // sanity check
  if( trained_.empty() ) return 0;

  // make
  SprTrainedRangeBooster* t = 0;
  if( SprEnforceLowMemory ) {
    t = new SprTrainedRangeBooster(trained_);
    for( int i=0;i<trained_.size();i++ ) trained_[i].second = false;
  }
  else {
    // prepare a vector of trained classifiers
    vector<pair<const SprAbsTrainedClassifier*,bool> > trained;
    for( int i=0;i<trained_.size();i++ ) {
      SprAbsTrainedClassifier* c = trained_[i].first->clone();
      trained.push_back(pair<const SprAbsTrainedClassifier*,bool>(c,true));
    }

    // make a trained bagger
    t = new SprTrainedRangeBooster(trained);
  }
  assert( t != 0 );

  // vars
  vector<string> vars;
  data_->vars(vars);
  t->setVars(vars);

  // exit
  return t;
}
