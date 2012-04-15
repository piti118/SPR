//$Id: SprArcE4.cc,v 1.7 2008-02-14 20:21:01 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprArcE4.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprEmptyFilter.hh"
#include "StatPatternRecognition/SprBootstrap.hh"
#include "StatPatternRecognition/SprAbsTwoClassCriterion.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"

#include <cassert>
#include <cmath>
#include <memory>

using namespace std;


SprArcE4::SprArcE4(SprAbsFilter* data, 
		   unsigned cycles, bool discrete)
  : 
  SprBagger(data,cycles,discrete), 
  initialDataWeights_(),
  response_(data->size(),pair<double,double>(0,0))
{ 
  data_->weights(initialDataWeights_);
  cout << "ArcE4 initialized." << endl;
}


bool SprArcE4::setData(SprAbsFilter* data)
{
  // reset base data
  if( !SprBagger::setData(data) ) {
    cerr << "Unable to set data for ArcE4." << endl;
    return false;
  }

  // copy weights
  data_->weights(initialDataWeights_);
  
  // init responses
  response_.clear();
  response_.resize(data_->size(),pair<double,double>(0,0));

  // exit
  return true;
}


bool SprArcE4::train(int verbose)
{
  // sanity check
  if( cycles_==0 || trainable_.empty() ) {
    cout << "ArcE4 will exit without training." << endl;
    return this->prepareExit(true);
  }

  // if resume training, generate a seed from time of day
  if( cycles_>0 && !trained_.empty() ) {
    delete bootstrap_;
    bootstrap_ = new SprBootstrap(data_,-1);
    assert( bootstrap_ != 0 );
  }

  // update responses
  assert( data_->size() == response_.size() );
  if( !trained_.empty() ) {
    for( int i=0;i<data_->size();i++ ) {
      const SprPoint* p = (*data_)[i];
      double resp = 0;
      for( int j=0;j<trained_.size();j++ )
	resp += trained_[j].first->response(p);
      double w = trained_.size();
      resp /= w;
      response_[i].first = resp;
      response_[i].second = w;
    }
  }

  // after all betas are filled, do an overall validation
  if( !this->initValBeta() )
    return this->prepareExit(false);

  // loop through trainable
  unsigned nCycle = 0;
  unsigned nFailed = 0;
  while( nCycle < cycles_ ) {
    for( int i=0;i<trainable_.size();i++ ) {
      // check cycles
      if( nCycle++ >= cycles_ ) return this->prepareExit((this->nTrained()>0));

      // generate replica
      auto_ptr<SprEmptyFilter> temp(bootstrap_->weightedReplica());
      if( temp->size() != data_->size() ) {
	cerr << "Failed to generate bootstrap replica." << endl;
	return this->prepareExit(false);
      }

      // get new classifier
      SprAbsClassifier* c = trainable_[i];
      if( !c->setData(temp.get()) ) {
	cerr << "Unable to set data for classifier " << i << endl;
	return this->prepareExit(false);
      }
      if( !c->train(verbose-1) ) {
	cerr << "ArcE4 failed to train classifier " << i 
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
	cerr << "ArcE4 failed to train classifier " << i 
	     << ". Continuing..."<< endl;
	if( ++nFailed >= cycles_ ) {
	  cout << "Exiting after failed to train " << nFailed 
	       << " classifiers." << endl;
	  return this->prepareExit((this->nTrained()>0));
	}
	else
	  continue;
      }
      trained_.push_back(pair<const SprAbsTrainedClassifier*,bool>(t,true));

      // reweight events
      this->reweight(t);
      if( verbose > 1 ) {
	cout << "After reweighting:   W1=" << data_->weightInClass(cls1_)
	     << " W0=" << data_->weightInClass(cls0_)
	     << "    N1=" << data_->ptsInClass(cls1_)
	     << " N0=" << data_->ptsInClass(cls0_) << endl;
      }

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


bool SprArcE4::prepareExit(bool status)
{
  // restore weights
  data_->setWeights(initialDataWeights_);

  // do basic restore
  return SprBagger::prepareExit(status);
}


void SprArcE4::reweight(const SprAbsTrainedClassifier* t)
{
  unsigned size = data_->size();
  assert( size == initialDataWeights_.size() );
  assert( size == response_.size() );
  for( int i=0;i<size;i++ ) {
    const SprPoint* p = (*data_)[i];

    // update response
    double& resp = response_[i].first;
    double& wresp = response_[i].second;
    resp = wresp*resp + t->response(p);
    wresp += 1.;
    resp /= wresp;

    // reweight
    int cls = -1;
    if(      p->class_ == cls0_ ) 
      cls = 0;
    else if( p->class_ == cls1_ )
      cls = 1;
    if( cls > -1 ) {
      double error = wresp * (resp - cls);
      double w = initialDataWeights_[i] * (1.+pow(fabs(error),4));
      data_->setW(i,w);
    }
  }
}
