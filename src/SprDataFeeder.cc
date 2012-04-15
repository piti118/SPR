//$Id: SprDataFeeder.cc,v 1.7 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprDataFeeder.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprAbsWriter.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprCoordinateMapper.hh"

#include <iostream>
#include <cassert>
#include <string>

using namespace std;


SprDataFeeder::~SprDataFeeder()
{
  delete mapper_;
  for( int i=0;i<specificMappers_.size();i++ ) 
    delete specificMappers_[i];
  for( int i=0;i<multiSpecificMappers_.size();i++ ) 
    delete multiSpecificMappers_[i];
}


SprDataFeeder::SprDataFeeder(const SprAbsFilter* data,
			     SprAbsWriter* writer,
			     const std::vector<unsigned>& mapper) 
  : 
  data_(data), 
  writer_(writer), 
  mode_(0),
  classifiers_(),
  multiclass_(),
  mapper_(SprCoordinateMapper::createMapper(mapper)),
  specificMappers_(),
  multiSpecificMappers_()
{
  assert( data_ != 0 );
  assert( writer_ != 0 );
  vector<string> axes;
  data_->vars(axes);
  writer_->setAxes(axes);
}


bool SprDataFeeder::addClassifier(const SprAbsTrainedClassifier* c,
				  const char* name,
				  const std::vector<unsigned>& mapper) 
{
  return this->addClassifier(c,name,SprCoordinateMapper::createMapper(mapper));
}


bool SprDataFeeder::addClassifier(const SprAbsTrainedClassifier* c,
				  const char* name,
				  SprCoordinateMapper* mapper) 
{ 
  if( c != 0 ) {
    // sanity check
    if( mode_ == 2 ) {
      cerr << "Unable to add classifier: " 
	   << "DataFeeder is in the multi class mode." << endl;
      return false;
    }
    else
      mode_ = 1;

    // check dimensionality
    if( c->dim()!=data_->dim() && mapper==0 ) {
      cerr << "Dimensionality of classifier " << name 
	   << " does not match dimensionality of input data." << endl;
      return false;
    }

    // add classifier
    classifiers_.push_back(c);
    specificMappers_.push_back(mapper);
    writer_->addAxis(name);
  }

  // exit
  return true;
}


bool SprDataFeeder::addMultiClassLearner(
			      const SprAbsTrainedMultiClassLearner* c,
			      const char* name,
			      const std::vector<unsigned>& mapper) 
{ 
  return this->addMultiClassLearner(c,name,
				    SprCoordinateMapper::createMapper(mapper));
}


bool SprDataFeeder::addMultiClassLearner(
				  const SprAbsTrainedMultiClassLearner* c,
				  const char* name,
				  SprCoordinateMapper* mapper) 
{
  if( c != 0 ) {
    // sanity check
    if( mode_ == 1 ) {
      cerr << "Unable to add multi class learner: " 
	   << "DataFeeder is in the regular classifier mode." << endl;
      return false;
    }
    else
      mode_ = 2;

    // check dimensionality
    if( c->dim()!=data_->dim() && mapper==0 ) {
      cerr << "Dimensionality of classifier " << name 
	   << " does not match dimensionality of input data." << endl;
      return false;
    }

    // add multi class learner
    string sname = name;
    multiclass_.push_back(c);
    multiSpecificMappers_.push_back(mapper);
    vector<string> axes;
    c->outputNames(sname.c_str(),axes);
    for( int i=0;i<axes.size();i++ )
      writer_->addAxis(axes[i].c_str());
    writer_->addAxis(sname.c_str());
  }

  // exit
  return true;
}


bool SprDataFeeder::feed(int nout) const
{
  // sanity checks
  assert( data_ != 0 );
  assert( writer_ != 0 );
  if( classifiers_.empty() && multiclass_.empty() ) {
    cout << "Warning: no classifiers specified for SprDataFeeder. " 
	 << "Will save data only." << endl;
  }

  // loop through data
  for( int i=0;i<data_->size();i++ ) {
    // message
    if( nout>0 && (i%nout)==0 ) {
      cout << "Feeder storing point " << i 
	   << " out of " << data_->size() << endl;
    }

    // get point
    const SprPoint* pTuple = (*data_)[i];
    const SprPoint* pResp = ( mapper_==0 ? pTuple : mapper_->output(pTuple) );

    // add classifiers
    vector<double> f;
    for( int j=0;j<classifiers_.size();j++ ) {
      const SprPoint* pSpecific = 0;
      if( specificMappers_[j] != 0 ) 
	pSpecific = specificMappers_[j]->output(pTuple);
      else
	pSpecific = pResp;
      assert( pSpecific != 0 );
      f.push_back(classifiers_[j]->response(pSpecific));
      if( specificMappers_[j] != 0 ) specificMappers_[j]->clear();
    }

    // add multi class learners
    for( int j=0;j<multiclass_.size();j++ ) {
      map<int,double> resp;
      const SprPoint* pSpecific = 0;
      if( multiSpecificMappers_[j] != 0 ) 
	pSpecific = multiSpecificMappers_[j]->output(pTuple);
      else
	pSpecific = pResp;
      assert( pSpecific != 0 );
      int cls = multiclass_[j]->response(pSpecific,resp);
      for( map<int,double>::const_iterator 
	     iter=resp.begin();iter!=resp.end();iter++ ) 
	f.push_back(iter->second);
      f.push_back(double(cls));
      if( multiSpecificMappers_[j] != 0 ) multiSpecificMappers_[j]->clear();
    }

    // write
    writer_->write(data_->w(i),pTuple,f);

    // clean up
    if( mapper_ != 0 ) mapper_->clear();
  }

  // exit
  bool closed = writer_->close();
  if( closed )
    cout << "Writer successfully closed." << endl;
  else
    cout << "Writer was unable to close successfully." << endl;
  return closed;
}
