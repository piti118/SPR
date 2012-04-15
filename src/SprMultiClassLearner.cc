//$Id: SprMultiClassLearner.cc,v 1.10 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprMultiClassLearner.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprAbsClassifier.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprClass.hh"

#include <fstream>
#include <iomanip>

using namespace std;


SprMultiClassLearner::~SprMultiClassLearner()
{
  this->destroy();
}


SprMultiClassLearner::SprMultiClassLearner(SprAbsFilter* data, 
					   SprAbsClassifier* c,
					   const std::vector<int>& classes,
					   const SprMatrix& indicator,
					   MultiClassMode mode)
  :
  SprAbsMultiClassLearner(data,classes),
  mode_(mode),
  indicator_(indicator),
  trainable_(c),
  trained_(),
  weights_(),
  keepWeights_(false),
  includeZeroIndicator_(true)
{
  assert( trainable_ != 0 );
  assert( mode_!=User || indicator_.num_col()!=0 );
  bool status = this->setClasses();
  assert( status );
}


void SprMultiClassLearner::setClassifierWeights(
			       const std::vector<double>& weights)
{
  assert( weights.size() == indicator_.num_col() );
  weights_ = weights;
  keepWeights_ = true;
}


bool SprMultiClassLearner::setClasses()
{
  // check classes
  if( !this->checkClasses() ) {
    cerr << "Incorrect class list supplied to SprMultiClassLearner." << endl;
    return false;
  }

  // check indicator matrix
  if( mode_ == User ) {
    if( indicator_.num_row() != mapper_.size() ) {
      cerr << "Number of rows of the indicator matrix is not equal " 
	   << "to the specified number of classes."<< endl;
      return false;
    }
  }

  // fill out indicator matrix
  if(      mode_ == OneVsAll ) {
    unsigned n = mapper_.size();
    SprMatrix mat(n,n,0);
    indicator_ = mat;
    for( int i=0;i<n;i++ ) { 
      for( int j=0;j<n;j++ ) indicator_[i][j] = -1;
    }
    for( int i=0;i<n;i++ ) indicator_[i][i] = 1;
  }
  else if( mode_ == OneVsOne ) {
    unsigned n = mapper_.size();
    unsigned m = n*(n-1)/2;
    SprMatrix mat(n,m,0);
    indicator_ = mat;
    int jstart = 0;
    int jend = 0;
    for( int i=0;i<n;i++ ) {
      jstart = jend;
      jend += n-1-i;
      for( int j=jstart;j<jend;j++ ) indicator_[i][j] = 1;
      int j = jstart;
      for( int k=i+1;k<n;k++ ) indicator_[k][j++] = -1;
    }
  }

  // show matrix
  this->printIndicatorMatrix(cout);

  // exit
  return true;
}


bool SprMultiClassLearner::train(int verbose)
{
  // reset
  this->destroy();
  trained_.resize(indicator_.num_col());

  // need to compute classifier weights?
  bool computeClassifierWeights = false;
  if( weights_.empty() ) {
    if( !includeZeroIndicator_ )
      computeClassifierWeights = true;
    weights_.resize(indicator_.num_col(),1.);
  }

  // build a map of classes
  map<int,unsigned> mapper;
  for( int i=0;i<mapper_.size();i++ )
    mapper.insert(pair<const int,unsigned>(mapper_[i],i));

  // loop thru columns of the indicator matrix
  for( int j=0;j<indicator_.num_col();j++ ) {
    vector<int> classes0, classes1;
    for( map<int,unsigned>::const_iterator
	   i=mapper.begin();i!=mapper.end();i++ ) {
      int cls = i->first;
      double flag = indicator_[i->second][j];
      if(      flag < -0.5 )
	classes0.push_back(cls);
      else if( flag > 0.5 )
	classes1.push_back(cls);
    }	
    SprClass cls0(classes0);
    SprClass cls1(classes1);

    // message
    double w0 = data_->weightInClass(cls0);
    double w1 = data_->weightInClass(cls1);
    if( computeClassifierWeights ) weights_[j] = w0+w1;
    unsigned n0 = data_->ptsInClass(cls0);
    unsigned n1 = data_->ptsInClass(cls1);
    if( verbose > 0 ) {
      cout << "Training classifier for matrix column " << j << " with "
	   << "     W0=" << w0 << " W1=" << w1
	   << "     N0=" << n0 << " N1=" << n1 << endl;
    }

    // apply a classifier
    if( !trainable_->reset() ) {
      cerr << "Unable to reset classifier for indicator column " << j << endl;
      return false;
    }
    if( !trainable_->setClasses(cls0,cls1) ) {
      cerr << "Unable to reset classes for indicator column " << j << endl;
      return false;
    }
    if( !trainable_->train(verbose) ) {
      cerr << "Unable to train classifier for indicator column " << j << endl;
      return false;
    }

    // make trained classifier
    const SprAbsTrainedClassifier* t = trainable_->makeTrained();
    if( t == 0 ) {
      cerr << "Unable to train classifier for column " << j << endl;
      return false;
    }
    trained_[j] = pair<const SprAbsTrainedClassifier*,bool>(t,true);
  }// end loop over columns

  // exit
  return true;
}


bool SprMultiClassLearner::reset()
{
  if( !trainable_->reset() ) {
    cerr << "Unable to reset trainable classifier." << endl;
    return false;
  }
  this->destroy();
  return true;
}


void SprMultiClassLearner::destroy()
{
  for( int i=0;i<trained_.size();i++ )
    if( trained_[i].second )
      delete trained_[i].first;
  trained_.clear();
  if( !keepWeights_ ) weights_.clear();
}


bool SprMultiClassLearner::setData(SprAbsFilter* data)
{
  assert( data != 0 );
  if( !trainable_->setData(data) ) {
    cerr << "Unable to set data for trainable classifier." << endl;
    return false;
  }
  data_ = data;
  this->destroy();
  return true;
}


void SprMultiClassLearner::print(std::ostream& os) const 
{
  os << "Trained MultiClassLearner " << SprVersion << endl;

  // print matrix
  this->printIndicatorMatrix(os);

  // print weights
  os << "Weights:";
  for( int i=0;i<weights_.size();i++ )
    os << " " << weights_[i];
  os << endl;

  // print classifiers
  for( int i=0;i<trained_.size();i++ ) {
    os << "Multi class learner classifier: " << i << endl;
    trained_[i].first->print(os);
  }
}


SprTrainedMultiClassLearner* SprMultiClassLearner::makeTrained()
{
  // make
  SprTrainedMultiClassLearner* t = 0;
  if( SprEnforceLowMemory ) {
    t = new  SprTrainedMultiClassLearner(indicator_,mapper_,trained_);
    for( int i=0;i<trained_.size();i++ ) trained_[i].second = false;
  }
  else {
    vector<pair<const SprAbsTrainedClassifier*,bool> > 
      trained(trained_.size());
    for( int i=0;i<trained_.size();i++ ) {
      trained[i] 
	= pair<const SprAbsTrainedClassifier*,bool>(trained_[i].first->clone(),
						    true);
    }
    t = new SprTrainedMultiClassLearner(indicator_,mapper_,trained);
  }
  assert( t != 0 );

  // weights
  t->setClassifierWeights(weights_);

  // treatment zero elements of the indicator matrix
  if( includeZeroIndicator_ )
    t->includeZeroClasses();
  else
    t->excludeZeroClasses();

  // set missing values
  if( !t->setDefaultMissing(validRange_,defaultMissing_) ) {
    cerr << "Unable to set medians for trained MultiClass learner." << endl;
    delete t;
    return 0;
  }

  // vars
  vector<string> vars;
  data_->vars(vars);
  t->setVars(vars);

  // exit
  return t;
}


void SprMultiClassLearner::printIndicatorMatrix(std::ostream& os) const
{
  os << "Indicator matrix:" << endl;
  os << setw(20) << "Classes/Classifiers" << " : " 
     << mapper_.size() << " " << indicator_.num_col() << endl;
  os << "=========================================================" << endl;
  for( int i=0;i<indicator_.num_row();i++ ) {
    os << setw(20) << mapper_[i] << " : ";
    for( int j=0;j<indicator_.num_col();j++ ) 
      os << setw(2) << indicator_[i][j] << " ";
    os << endl;
  }
  os << "=========================================================" << endl;
}


void SprMultiClassLearner::setTrained(const SprMatrix& indicator, 
				      const std::vector<int>& classes,
				      const std::vector<std::pair<
				      const SprAbsTrainedClassifier*,bool> >& 
				      trained,
				      const std::vector<double>& weights) 
{
  indicator_ = indicator;
  mode_ = User;
  mapper_ = classes;
  trained_ = trained;
  assert( mapper_.size() == indicator_.num_row() );
  assert( trained_.size() == indicator_.num_col() );
  assert( !mapper_.empty() );
  assert( !trained_.empty() );
  weights_ = weights;
  assert( weights_.size() == trained_.size() );
}
