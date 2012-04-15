// $Id: SprInputNormalizer.cc,v 1.3 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh" 
#include "StatPatternRecognition/SprAbsFilter.hh" 
#include "StatPatternRecognition/SprInputNormalizer.hh" 
#include "StatPatternRecognition/SprDataMoments.hh" 
#include "StatPatternRecognition/SprDefs.hh"
#include "math/SprVector.hh"
#include "math/SprSymMatrix.hh"

#include <cassert>
#include <cmath>
#include <algorithm>

using namespace std;


SprInputNormalizer::SprInputNormalizer()
  :
  SprAbsVarTransformer(),
  mean_(),
  sigma_()
{}


SprInputNormalizer::SprInputNormalizer(const std::vector<double>& mean,
				       const std::vector<double>& sigma)
  :
  SprAbsVarTransformer(),
  mean_(mean),
  sigma_(sigma)
{}


SprInputNormalizer::SprInputNormalizer(const SprInputNormalizer& other)
  :
  SprAbsVarTransformer(other),
  mean_(other.mean_),
  sigma_(other.sigma_)
{}  


bool SprInputNormalizer::allVarsIndependent() const
{
  if( oldVars_.size() != newVars_.size() ) {
    cerr << "Variable lists sizes do not match in "
	 << "SprInputNormalizer::allVarsIndependent()" << endl;
    return false;
  }
  for( int i=0;i<oldVars_.size();i++ ) {
    if( oldVars_[i] != newVars_[i] ) {
      cerr << "Variable mismatch in SprInputNormalizer::allVarsIndependent(): "
	   << oldVars_[i] << " " << newVars_[i] << endl;
      return false;
    }
  }
  return true;
}


bool SprInputNormalizer::reduceVars(const std::vector<std::string>& vars)
{
  // sanity check
  if( !this->allVarsIndependent() ) {
    cerr << "SprInputNormalizer cannot reduce variable list. "
	 << "Independence check fails." << endl;
    return false;
  }

  // prepare new vectors
  vector<string> newVars;
  vector<double> mean, sigma;
  int N = oldVars_.size();
  assert( mean_.size() == N );
  assert( sigma_.size() == N );
  for( int i=0;i<N;i++ ) {
    if( find(vars.begin(),vars.end(),oldVars_[i]) != vars.end() ) {
      newVars.push_back(oldVars_[i]);
      mean.push_back(mean_[i]);
      sigma.push_back(sigma_[i]);
    }
  }

  // reassign
  mean_ = mean;
  sigma_ = sigma;
  oldVars_ = newVars;
  newVars_ = newVars;

  // exit
  return true;
}


bool SprInputNormalizer::train(const SprAbsFilter* data, int verbose)
{
  // get vars
  data->vars(oldVars_);
  assert( !oldVars_.empty() );
  newVars_ = oldVars_;
  int dim = oldVars_.size();

  // compute moments
  SprDataMoments moms(data);
  SprSymMatrix cov;
  SprVector mean;
  if( !moms.covariance(cov,mean) ) {
    cerr << "Unable to compute mean and covariance for  SprInputNormalizer."
	 << endl;
    return false;
  }
  assert( cov.num_row() == dim );
  assert( mean.num_row() == dim );
  mean_ = mean.std();
  sigma_.clear();
  sigma_.resize(dim);
  for( int i=0;i<dim;i++ )
    sigma_[i] = ( cov[i][i]>0 ? sqrt(cov[i][i]) : 0 );

  // exit
  return true;
}


void SprInputNormalizer::transform(const std::vector<double>& in, 
				    std::vector<double>& out) const
{
  assert( in.size() == mean_.size() );
  out.clear();
  out.resize(in.size());
  for( int i=0;i<in.size();i++ ) {
    if( sigma_[i] > 0 )
      out[i] = (in[i]-mean_[i])/sigma_[i];
  }
}


void SprInputNormalizer::inverse(const std::vector<double>& in, 
				 std::vector<double>& out) const
{
  assert( in.size() == mean_.size() );
  out.clear();
  out.resize(in.size());
  for( int i=0;i<in.size();i++ ) {
    if( sigma_[i] > 0 )
      out[i] = in[i]*sigma_[i] + mean_[i];
  }
}


void SprInputNormalizer::print(std::ostream& os) const
{
  // save name
  os << "VarTransformer: " << this->name().c_str() 
     << " " << SprVersion.c_str() << endl;

  // init
  int dim = oldVars_.size();
  vector<string> vars(oldVars_);

  // protect againt spaces in var names
  for( int i=0;i<dim;i++ ) {
    if( vars[i].find(' ') != string::npos )
	vars[i].erase(vars[i].find_first_of(' '));
  }

  // store dimensionality
  os << "Dim: " << dim << endl;

  // store transformations
  for( int i=0;i<dim;i++ ) {
    os << i << " " << vars[i].c_str() << "=("
       << vars[i].c_str() << "- " << mean_[i] << " )/ " << sigma_[i] << endl;
  }
}


