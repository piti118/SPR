//$Id: SprTrainedBinaryEncoder.cc,v 1.1 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprTrainedBinaryEncoder.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"

#include <cassert>
#include <utility>
#include <algorithm>
#include <functional>

using namespace std;


struct STBECmpPairSecond
  : public binary_function<pair<const int,double>,
			   pair<const int,double>,bool> {
  bool operator()(const pair<const int,double>& l, 
		  const pair<const int,double>& r) const {
    return (l.second < r.second);
  }
};


SprTrainedBinaryEncoder::SprTrainedBinaryEncoder(
			  const std::vector<int>& mapper,
			  const SprAbsTrainedClassifier* classifier,
			  bool ownClassifier)
  :
  SprAbsTrainedMultiClassLearner(mapper),
  classifier_(classifier),
  ownClassifier_(ownClassifier)
{
  assert( classifier_ != 0 );
  bool normalized = classifier_->normalized();
  assert( normalized );
}


SprTrainedBinaryEncoder::SprTrainedBinaryEncoder(
			       const SprTrainedBinaryEncoder& other)
  :
  SprAbsTrainedMultiClassLearner(other),
  classifier_(other.classifier_->clone()),
  ownClassifier_(true)
{
  assert( classifier_ != 0 );
  bool normalized = classifier_->normalized();
  assert( normalized );
}


void SprTrainedBinaryEncoder::destroy()
{
  if( ownClassifier_ ) {
    delete classifier_;
    classifier_ = 0;
  }
}


int SprTrainedBinaryEncoder::response_one(const std::vector<double>& input,
					  std::map<int,double>& output) const
{
  // init
  output.clear();

  // loop thru classes
  for( int ic=0;ic<mapper_.size();ic++ ) {
    vector<double> v = input;
    v.push_back(mapper_[ic]);
    output.insert(pair<const int,double>(mapper_[ic],
					 1.-classifier_->response(v)));
  }

  // find minimal loss
  map<int,double>::const_iterator iter 
    = min_element(output.begin(),output.end(),STBECmpPairSecond());
  return iter->first;
}


void SprTrainedBinaryEncoder::print(std::ostream& os) const 
{
  os << "Trained BinaryEncoder " << SprVersion << endl;

  // print classes
  os << "Classes: " << mapper_.size() << endl;
  for( int ic=0;ic<mapper_.size();ic++ )
    os << mapper_[ic] << " ";
  os << endl;

  // print classifier
  classifier_->print(os);
}

