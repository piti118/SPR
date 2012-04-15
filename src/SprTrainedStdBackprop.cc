//$Id: SprTrainedStdBackprop.cc,v 1.7 2008-01-18 22:56:59 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprTrainedStdBackprop.hh"
#include "StatPatternRecognition/SprTransformation.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <cmath>
#include <iomanip>
#include <cassert>
#include <algorithm>

using namespace std;


SprTrainedStdBackprop::SprTrainedStdBackprop()
  : 
  SprAbsTrainedClassifier()
  , nNodes_(0)
  , nLinks_(0)
  , structure_()
  , nodeType_()
  , nodeActFun_()
  , nodeNInputLinks_()
  , nodeFirstInputLink_()
  , linkSource_()
  , nodeBias_()
  , linkWeight_()
{
  this->setCut(SprUtils::lowerBound(0.5));
}


SprTrainedStdBackprop::SprTrainedStdBackprop(
			const char* structure,
                        const std::vector<SprNNDefs::NodeType>& nodeType,
			const std::vector<SprNNDefs::ActFun>& nodeActFun,
			const std::vector<int>& nodeNInputLinks,
			const std::vector<int>& nodeFirstInputLink,
			const std::vector<int>& linkSource,
			const std::vector<double>& nodeBias,
			const std::vector<double>& linkWeight)
  :
  SprAbsTrainedClassifier(),
  nNodes_(0),
  nLinks_(0),
  structure_(structure),
  nodeType_(nodeType),
  nodeActFun_(nodeActFun),
  nodeNInputLinks_(nodeNInputLinks),
  nodeFirstInputLink_(nodeFirstInputLink),
  linkSource_(linkSource),
  nodeBias_(nodeBias),
  linkWeight_(linkWeight)
{
  nNodes_ = nodeType_.size();
  assert( nNodes_ == nodeActFun_.size() );
  assert( nNodes_ == nodeNInputLinks_.size() );
  assert( nNodes_ == nodeFirstInputLink_.size() );
  assert( nNodes_ == nodeBias_.size() );
  nLinks_ = linkSource_.size();
  assert( nLinks_ == linkWeight_.size() );
  this->setCut(SprUtils::lowerBound(0.5));
}


SprTrainedStdBackprop::SprTrainedStdBackprop(
			  const SprTrainedStdBackprop& other)
  :
  SprAbsTrainedClassifier(other)
  , nNodes_(other.nNodes_)
  , nLinks_(other.nLinks_)
  , structure_(other.structure_)
  , nodeType_(other.nodeType_)
  , nodeActFun_(other.nodeActFun_)
  , nodeNInputLinks_(other.nodeNInputLinks_)
  , nodeFirstInputLink_(other.nodeFirstInputLink_)
  , linkSource_(other.linkSource_)
  , nodeBias_(other.nodeBias_)
  , linkWeight_(other.linkWeight_)
{}


double SprTrainedStdBackprop::activate(double x, SprNNDefs::ActFun f) const 
{
  switch (f) 
    {
    case SprNNDefs::ID :
      return x;
      break;
    case SprNNDefs::LOGISTIC :
      return SprTransformation::logit(x);
      break;
    default :
      cerr << "FATAL ERROR: Unknown activation function " 
	   << f << " in SprTrainedStdBackprop::activate" << endl;
      return 0;
    }
  return 0;
}


void SprTrainedStdBackprop::print(std::ostream& os) const 
{
  os << "Trained StdBackprop with configuration " 
     << structure_.c_str() << " " << SprVersion << endl; 
  os << "Activation functions: Identity=1, Logistic=2" << endl;
  os << "Cut: " << cut_.size();
  for( int i=0;i<cut_.size();i++ )
    os << "      " << cut_[i].first << " " << cut_[i].second;
  os << endl;
  os << "Nodes: " << nNodes_ << endl;
  for( int i=0;i<nNodes_;i++ ) {
    char nodeType;
    switch( nodeType_[i] )
      {
      case SprNNDefs::INPUT :
	nodeType = 'I';
	break;
      case SprNNDefs::HIDDEN :
	nodeType = 'H';
	break;
      case SprNNDefs::OUTPUT :
	nodeType = 'O';
	break;
      }
    int actFun = 0;
    switch( nodeActFun_[i] )
      {
      case SprNNDefs::ID :
	actFun = 1;
	break;
      case SprNNDefs::LOGISTIC :
	actFun = 2;
	break;
      }
    os << setw(6) << i
       << "    Type: "           << nodeType
       << "    ActFunction: "    << actFun
       << "    NInputLinks: "    << setw(6) << nodeNInputLinks_[i]
       << "    FirstInputLink: " << setw(6) << nodeFirstInputLink_[i]
       << "    Bias: "           << nodeBias_[i]
       << endl;
  }
  os << "Links: " << nLinks_ << endl;
  for( int i=0;i<nLinks_;i++ ) {
    os << setw(6) << i
       << "    Source: " << setw(6) << linkSource_[i]
       << "    Weight: " << linkWeight_[i]
       << endl;
  }
}


double SprTrainedStdBackprop::response(const std::vector<double>& v) const
{
  // Initialize and process input nodes
  vector<double> nodeOut(nNodes_,0);
  int d = 0;
  for( int i=0;i<nNodes_;i++ ) {
    if( nodeType_[i] == SprNNDefs::INPUT ) {
      assert( d < v.size() );
      nodeOut[i] = v[d++];
    }
    else
      break;
  }
  assert( d == v.size() );

  // Process hidden and output nodes
  for( int i=0;i<nNodes_;i++ ) {
    double nodeAct = 0;
    if( nodeNInputLinks_[i] > 0 ) {
      for( int j=nodeFirstInputLink_[i];
	   j<nodeFirstInputLink_[i]+nodeNInputLinks_[i];j++ ) {
        nodeAct += nodeOut[linkSource_[j]] * linkWeight_[j];
      }
      nodeOut[i] = this->activate(nodeAct+nodeBias_[i],nodeActFun_[i]);
    }
  }

  // Find output node and return result
  return nodeOut[nNodes_-1];
}


bool SprTrainedStdBackprop::varImportance(std::vector<
					  std::pair<std::string,double> >& 
					  importance) const
{
  // init
  importance.clear();

  // count input nodes
  unsigned dim = 0;
  for( int i=0;i<nNodes_;i++ ) {
    if( nodeType_[i] == SprNNDefs::INPUT ) dim++;
  }
  if( dim == 0 ) {
    cerr << "Unable to find input nodes in the neural net." << endl;
    return false;
  }

  // init structure for holding importance:
  // vector of size dim at each node
  vector<vector<double> > imp(nNodes_,vector<double>(dim,0));
  int d = 0;
  for( int i=0;i<nNodes_;i++ ) {
    if( nodeType_[i] == SprNNDefs::INPUT ) imp[i][d++] = 1;
  }

  // update importance moving thru the network
  for( int i=0;i<nNodes_;i++ ) {
    if( nodeNInputLinks_[i] > 0 ) {
      for( int j=nodeFirstInputLink_[i];
	   j<nodeFirstInputLink_[i]+nodeNInputLinks_[i];j++ ) {
	int k = linkSource_[j];
	double w = fabs(linkWeight_[j]);
	for( int d=0;d<dim;d++ )
	  imp[i][d] += w*imp[k][d];
      }
    }
  }

  // normalize
  vector<double>& imp_norm = imp[nNodes_-1];
  /*
  double sum = 0;
  for( int i=0;i<imp_norm.size();i++ ) sum += imp_norm[i];
  if( sum < SprUtils::eps() ) {
    cerr << "Importance of all input variables is too small." << endl;
    return false;
  }
  for( int i=0;i<imp_norm.size();i++ ) imp_norm[i] /= sum;  
  */

  // compose
  assert( dim == vars_.size() );
  importance.resize(dim);
  for( int d=0;d<dim;d++ )
    importance[d] = pair<string,double>(vars_[d],imp_norm[d]);

  // exit
  return true;
}
