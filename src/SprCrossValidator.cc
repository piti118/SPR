//$Id: SprCrossValidator.cc,v 1.9 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprCrossValidator.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprEmptyFilter.hh"
#include "StatPatternRecognition/SprAbsTwoClassCriterion.hh"
#include "StatPatternRecognition/SprAverageLoss.hh"
#include "StatPatternRecognition/SprAbsClassifier.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprAbsMultiClassLearner.hh"
#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprIntegerPermutator.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprFomCalculator.hh"
#include "StatPatternRecognition/SprPlotter.hh"
#include "StatPatternRecognition/SprMultiClassPlotter.hh"

#include <list>
#include <iostream>
#include <cmath>
#include <sstream>

using namespace std;


SprCrossValidator::~SprCrossValidator()
{
  for( int i=0;i<samples_.size();i++ )
    delete samples_[i];
}


SprCrossValidator::SprCrossValidator(const SprAbsFilter* data, 
				     unsigned nPieces)
  : 
  data_(data),
  samples_()
{
  bool status = this->divide(nPieces);
  assert( status );
}


bool SprCrossValidator::divide(unsigned nPieces)
{
  // get data size
  unsigned N = data_->size();

  // get total weight
  vector<double> weights;
  data_->weights(weights);
  double totW = 0;
  for( int i=0;i<weights.size();i++ )
    totW += weights[i];

  // find all classes in input data
  vector<SprClass> classes;
  data_->allClasses(classes);
  assert( !classes.empty() );

  // get weights for each class and perform sanity check
  int nClass = classes.size();  
  vector<double> totWperClass(nClass);
  for( int ic=0;ic<nClass;ic++ )
    totWperClass[ic] = data_->weightInClass(classes[ic]);

  // randomize point indices
  vector<unsigned> index;
  SprIntegerPermutator permu(N);
  if( !permu.sequence(index) ) {
    cerr << "CrossValidator is unable to randomize indices." << endl;
    return false;
  }

  // fill subsamples
  vector<int> pointAssigned(N,0);
  int nstart = 0;
  double accountedW = 0;
  vector<double> accountedWperClass(nClass,0);
  //
  for( int ip=0;ip<nPieces;ip++ ) {
    // last piece
    bool lastPiece = ((ip+1) == nPieces);

    // subsamples own SprData which does not own points
    SprData* sample = data_->emptyCopy();
    SprEmptyFilter* filter = new SprEmptyFilter(sample,classes,true);

    // prepare for loop over events
    weights.clear();
    double weightInPiece = 0;
    bool oneClassAttempted = false;
    vector<double> weightInClass(nClass,0);
    vector<int> classFilled(nClass,0), classAttempted(nClass,0);
    bool canResetStart = true;
    bool endEventProcessing = false;

    // compute max allowed weight per piece
    double maxWperPiece = (totW-accountedW) / (nPieces-ip);
    vector<double> maxWperClass(nClass);
    for( int ic=0;ic<nClass;ic++ ) {
      maxWperClass[ic] = (totWperClass[ic]-accountedWperClass[ic]) 
	/ (nPieces-ip);
    }

    // start event processing
    for( int n=nstart;n<N;n++ ) {
      // end processing? exit event loop
      if( endEventProcessing ) break;

      // do nothing if the point is assigned
      if( pointAssigned[n] == 1 ) continue;

      // get the point
      int ind = index[n];
      SprPoint* p = (*data_)[ind];
      double w = data_->w(ind);

      // check accumulated weights in classes
      for( int ic=0;ic<nClass;ic++ ) {
	// class filled? do not process
	if( classFilled[ic] == 1 ) continue;

	if( p->class_ == classes[ic] ) {
	  // update weights
	  weightInClass[ic] += w;
	  weightInPiece += w;

	  // check total weight
	  if( weightInPiece>maxWperPiece 
	      && oneClassAttempted && !lastPiece ) {
	    if( canResetStart ) nstart = n;
	    endEventProcessing = true;
	    break;
	  }
	  oneClassAttempted = true;

	  // check weight in this class
	  if( weightInClass[ic]>maxWperClass[ic] 
	      && classAttempted[ic]==1 && !lastPiece ) {
	    classFilled[ic] = 1;
	    if( canResetStart ) {
	      nstart = n;
	      canResetStart = false;
	    }
	  }
	  // insert point and update accounted weights
	  else {
	    sample->uncheckedInsert(p);
	    weights.push_back(w);
	    accountedW += w;
	    accountedWperClass[ic] += w;
	    pointAssigned[n] = 1;
	  }
	  classAttempted[ic] = 1;

	  // exit class loop since the correct class has been found
	  break;
	}// end if p->class_ == classes[ic]
      }// end loop over classes
    }// end loop over points

    // add to the list of subsamples or clean up
    if( !filter->empty() ) {
      if( !filter->setWeights(weights) ) {
	cerr << "Unable to set weights for subsample " << ip << endl;
	return false;
      }
      samples_.push_back(filter);
    }
    else
      delete filter;
  }

  // sanity check
  if( samples_.size() < 2 ) {
    cerr << "Not enough samples for cross-validation: " 
	 << samples_.size() << endl;
    return false;
  }

  // exit
  return true;
}


bool SprCrossValidator::validate(const SprAbsTwoClassCriterion* crit,
				 SprAverageLoss* loss,
				 const std::vector<SprAbsClassifier*>& 
				 classifiers,
				 const SprClass& cls0, const SprClass& cls1,
				 std::vector<SprValueWithError>& crossFom,
				 bool integrate,
				 int verbose) const
{
  // print out
  if( verbose > 0 ) {
    cout << "Will cross-validate using " 
	 << samples_.size() << " subsamples: " << endl;
    for( int i=0;i<samples_.size();i++ ) {
      cout << "Subsample " << i 
	   << "  W1=" << samples_[i]->weightInClass(cls1)
	   << "  W0=" << samples_[i]->weightInClass(cls0)
	   << "  N1=" << samples_[i]->ptsInClass(cls1)
	   << "  N0=" << samples_[i]->ptsInClass(cls0) << endl;
    }
  }

  // sanity check
  assert( !classifiers.empty() && !samples_.empty() );

  // make classes
  vector<SprClass> classes(2);
  classes[0] = cls0;
  classes[1] = cls1;

  // make a local copy of data
  SprEmptyFilter data(data_);

  // init
  crossFom.clear();
  crossFom.resize(classifiers.size());

  // loop over classifiers
  for( int ic=0;ic<classifiers.size();ic++ ) {
    SprAbsClassifier* c = classifiers[ic];
    assert( c != 0 );
    double confusion[2][2];
    confusion[0][0]=0.;confusion[0][1]=0.;
    confusion[1][0]=0.;confusion[1][1]=0.; 
    // message
    if( verbose > 0 )
      cout << "Cross-validator processing classifier " << ic << endl;

    // loop over subsamples
    vector<vector<SprPlotter::Response> > responses(samples_.size());
    for( int is=0;is<samples_.size();is++ ) {
      // message
      if( verbose > 0 ) {
	cout << "Cross-validator processing sub-sample " << is 
	     << " for classifier " << ic << endl;
      }

      // remove subsample from training data
      data.clear();
      data.remove(samples_[is]->data());
      vector<SprClass> missing;
      if( !data.checkClasses(classes,missing) ) {
	if( verbose > 0 ) {
	  cout << "Training data for CV set " << is 
	       << " is missing input classes:   ";
	  for( int ic=0;ic<missing.size();ic++ )
	    cout << " " << missing[ic];
	  cout << endl;
	}
	continue;
      }

      // print out
      if( verbose > 1 ) {
	cout << "Will train classifier " << c->name().c_str()
	     << " on a sample: " << endl;
	cout << "  W1=" << data.weightInClass(cls1)
	     << "  W0=" << data.weightInClass(cls0)
	     << "  N1=" << data.ptsInClass(cls1)
	     << "  N0=" << data.ptsInClass(cls0) << endl;
      }

      // reset classifier
      if( !c->setData(&data) ) {
	cerr << "Cross-validator unable to set data for classifier " 
	     << ic << endl;
	return false;
      }

      // train
      if( !c->train(verbose-1) ) {
	cerr << "Unable to train classifier " << ic << endl;
	continue;
      }
      const SprAbsTrainedClassifier* trained = c->makeTrained();
      if( trained == 0 ) {
	cerr << "Cross-validator unable to get trained classifier "
	     << ic << " for subsample " << is << endl;
	continue;
      }

      // loop thru subsample and fill out trained responses
      const SprAbsFilter* sub = samples_[is];
      for( int i=0;i<sub->size();i++ ) {
	const SprPoint* p = (*sub)[i];
	int cls = -1;
	if(      p->class_ == cls0 )
	  cls = 0;
	else if( p->class_ == cls1 )
	  cls = 1;
	else
	  continue;
	double w = sub->w(i);
	SprPlotter::Response resp(cls,w);
	double r = 0;
	int accept = ( trained->accept(p,r) ? 1 : 0 );
        //cout << "Accept: " << accept << "(" << cls << ")" << endl;
	ostringstream ost;
	ost << i;
	resp.set(ost.str().c_str(),r,accept);
	responses[is].push_back(resp);
        confusion[cls][accept]+=w;
      }// end loop over points in this subsample

      // clean up
      delete trained;
    }// end loop over subsamples

    // fill cross-validation FOM
    crossFom[ic] = SprFomCalculator::fom(responses,crit,loss,
    					 integrate,verbose-1);

    cout << "confusion matrix" << endl;
    cout << "0=>0: " << confusion[0][0] << endl; 
    cout << "0=>1: " << confusion[0][1] << endl; 
    cout << "1=>0: " << confusion[1][0] << endl; 
    cout << "1=>1: " << confusion[1][1] << endl; 
    
    // reset classifier to point to the original data
    if( !c->setData(const_cast<SprAbsFilter*>(data_)) ) {
      cerr << "Cross-validator unable to restore data for classifier " 
	   << ic << endl;
      return false;
    }
  }// end loop over classifiers

  // exit
  return true;
}


bool SprCrossValidator::validate(SprAbsMultiClassLearner* mcLearner,
				 SprAverageLoss* loss,
				 SprValueWithError& realLoss,
				 SprClassificationTable& classificationTable,
				 std::map<int,double>& weightInClass,
				 int verbose) const
{
  // get classes
  vector<int> intClasses;
  mcLearner->classes(intClasses);
  assert( !intClasses.empty() );
  vector<SprClass> classes(intClasses.size());
  for( int i=0;i<intClasses.size();i++ )
    classes[i] = intClasses[i];

  // print out
  if( verbose > 0 ) {
    cout << "Will cross-validate using " 
	 << samples_.size() << " subsamples: " << endl;
    for( int i=0;i<samples_.size();i++ ) {
      cout << "Subsample " << i << endl;
      double totW = 0;
      unsigned totN = 0;
      for( int ic=0;ic<classes.size();ic++ ) {
	double w   = samples_[i]->weightInClass(classes[ic]);
	unsigned n = samples_[i]->ptsInClass(classes[ic]);
	totW += w;
	totN += n;
	cout << "    Class " << classes[ic]
	     << "  W=" << w << "  N=" << n << endl;
      }
      cout << "Total W=" << totW << "     Total N=" << totN << endl;
    }
  }

  // sanity check
  assert( mcLearner!=0 && !samples_.empty() );

  // make a local copy of data
  SprEmptyFilter data(data_);

  // loop over subsamples
  vector<SprMultiClassPlotter::Response> responses;
  for( int is=0;is<samples_.size();is++ ) {
    // message
    if( verbose > 0 ) {
      cout << "Cross-validator processing sub-sample " << is 
	   << " for MultiClass learner " << endl;
    }

    // remove subsample from training data
    data.clear();
    data.remove(samples_[is]->data());
    vector<SprClass> missing;
    if( !data.checkClasses(classes,missing) ) {
      if( verbose > 0 ) {
	cout << "Training data for CV set " << is 
	     << " is missing input classes:   ";
	for( int ic=0;ic<missing.size();ic++ )
	  cout << " " << missing[ic];
	cout << endl;
      }
      continue;
    }

    // reset classifier
    if( !mcLearner->setData(&data) ) {
      cerr << "Cross-validator unable to set data for MultiClass learner." 
	   << endl;
      return false;
    }

    // train
    if( !mcLearner->train(verbose-1) ) {
      cerr << "Unable to train MultiClass learner." << endl;
      continue;
    }
    const SprAbsTrainedMultiClassLearner* trained = mcLearner->makeTrained();
    if( trained == 0 ) {
      cerr << "Cross-validator unable to get trained MultiClass learner "
	   << " for subsample " << is << endl;
      continue;
    }

    // loop thru subsample and fill out trained responses
    const SprAbsFilter* sub = samples_[is];
    for( int i=0;i<sub->size();i++ ) {
      const SprPoint* p = (*sub)[i];
      int cls = p->class_;
      if( find(intClasses.begin(),intClasses.end(),cls) 
	  == intClasses.end() ) continue;
      double w = sub->w(i);
      map<int,double> output;
      int r = trained->response(p,output);
      responses.push_back(SprMultiClassPlotter::Response(cls,w,r,output));
    }// end loop over points in this subsample

    // clean up
    delete trained;
  }// end loop over subsamples

  // fill cross-validation FOM
  realLoss = SprFomCalculator::loss(intClasses,
				    responses,
				    loss,
				    classificationTable,
				    weightInClass,
				    verbose-1);

  // reset classifier to point to the original data
  if( !mcLearner->setData(const_cast<SprAbsFilter*>(data_)) ) {
    cerr << "Cross-validator unable to restore data for MultiClass learner." 
	 << endl;
    return false;
  }

  // exit
  return true;
}
