//$Id: SprBinaryEncoder.cc,v 1.2 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprBinaryEncoder.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprEmptyFilter.hh"
#include "StatPatternRecognition/SprAbsClassifier.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"

#include <algorithm>

using namespace std;


SprBinaryEncoder::~SprBinaryEncoder()
{
  this->destroy();
}


SprBinaryEncoder::SprBinaryEncoder(SprAbsFilter* data, 
				   SprAbsClassifier* c,
				   const std::vector<int>& classes)
  :
  SprAbsMultiClassLearner(data,classes),
  convertedData_(0),
  trainable_(c),
  trained_(0)
{
  assert( trainable_ != 0 );
}


void SprBinaryEncoder::destroy()
{
  delete convertedData_;
  convertedData_ = 0;
  delete trained_;
  trained_ = 0;
}


bool SprBinaryEncoder::train(int verbose)
{
  // if data has not been set, do it now
  if( convertedData_ == 0 ) {
    if( !this->setData(data_) ) {
      cerr << "Unable to set training data for SprBinaryEncoder." << endl;
      return false;
    }
  }

  // train
  if( !trainable_->train(verbose) ) {
    cerr << "Unable to train binary classifier for SprBinaryEncoder." << endl;
    return false;
  }

  // make trained classifier
  trained_ = trainable_->makeTrained();
  if( trained_ == 0 ) {
    cerr << "Unable to make trained binaruy classifier for SprBinaryEncoder." 
	 << endl;
    return false;
  }

  // exit
  return true;
}


bool SprBinaryEncoder::reset()
{
  delete trained_;
  trained_ = 0;
  if( trainable_->reset() ) {
    cerr << "Unable to reset trainable classifier for SprBinaryEncoder." 
	 << endl;
    return false;
  }
  return true;
}


bool SprBinaryEncoder::setData(SprAbsFilter* data)
{
  // set global data
  data_ = data;

  // clear data
  if( convertedData_ != 0 ) delete convertedData_;

  // convert data into binary format
  convertedData_ = this->convertData();
  if( convertedData_ == 0 ) {
    cerr << "Unable to convert data for SprBinaryEncoder." << endl;
    return false;
  }

  // set data for the trainable classifier
  if( !trainable_->setData(convertedData_) ) {
    cerr << "Unable to set data for binary classifier for SprBinaryEncoder." 
	 << endl;
    return false;
  }

  // set classes for trainable classifier
  SprClass cls0(0), cls1(1);
  if( !trainable_->setClasses(cls0,cls1) ) {
    cerr << "Unable to set classes for binary classifier for SprBinaryEncoder."
	 << endl;
    return false;
  }

  // clenup trained
  delete trained_;
  trained_ = 0;

  // exit
  return true;
}


void SprBinaryEncoder::print(std::ostream& os) const
{
  assert( trained_ != 0 );
  os << "Trained BinaryEncoder " << SprVersion << endl;

  // print classes
  os << "Classes: " << mapper_.size() << endl;
  for( int ic=0;ic<mapper_.size();ic++ )
    os << mapper_[ic] << " ";
  os << endl;

  // print classifier
  trained_->print(os);
}


SprTrainedBinaryEncoder* SprBinaryEncoder::makeTrained()
{
  if( trained_ == 0 ) return 0;
  SprTrainedBinaryEncoder* t = 0;
  if( SprEnforceLowMemory ) {
    t = new SprTrainedBinaryEncoder(mapper_,trained_,true);
    trained_ = 0;
  }
  else {
    t = new SprTrainedBinaryEncoder(mapper_,trained_->clone(),true);
  }
  assert( t != 0 );
  return t;
}


SprEmptyFilter* SprBinaryEncoder::convertData() const
{
  // get variable list
  vector<string> vars;
  data_->vars(vars);
  string sclass = "TrueClass";
  vector<string>::const_iterator found = find(vars.begin(),vars.end(),sclass);
  if( found != vars.end() ) {
    cerr << "Variable " << sclass.c_str() << " is already included " 
	 << "in the input list for SprBinaryEncoder." << endl;
    return false;
  }
  vars.push_back(sclass);

  // make new data
  bool ownPoints = true;
  SprData* newData = new SprData("ConvertedDataWithTrueClass",vars,ownPoints);

  // get the number of classes
  double wSignal = mapper_.size();

  // add points expanding with true class
  vector<double> weights;
  for( int i=0;i<data_->size();i++ ) {
    const SprPoint* p = (*data_)[i];
    if( find(mapper_.begin(),mapper_.end(),p->class_) == mapper_.end() )
      continue;
    double w = data_->w(i);
    for( int ic=0;ic<mapper_.size();ic++ ) {
      int icls = mapper_[ic];
      vector<double> x = p->x_;
      x.push_back(icls);
      int cls = ( icls==p->class_ ? 1 : 0 );
      newData->insert(p->index_,cls,x);
      double factor = ( cls==0 ? 1. : wSignal );
      weights.push_back(w*factor);
    }
  }

  // make new filter
  bool ownData = true;
  return new SprEmptyFilter(newData,weights,ownData);
}
