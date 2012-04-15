//$Id: SprIntegerBootstrap.cc,v 1.4 2008-02-06 22:01:57 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprIntegerBootstrap.hh"

#include <vector>

using namespace std;


void SprIntegerBootstrap::set(unsigned dim, unsigned nsample)
{
  assert( dim > 0 );
  assert( nsample > 0 );
  dim_ = dim;
  nsample_ = nsample;
}


bool SprIntegerBootstrap::replica(std::vector<unsigned>& v, int npts)
{
  // init
  if( npts <= 0 ) npts = nsample_;

  // reset input vector
  v.clear();

  // make array of uniform random numbers on [0,1]
  double* r = new double [npts];
  generator_.sequence(r,npts);
  unsigned iuse = 0;
  for( int i=0;i<npts;i++ ) {
    iuse = unsigned(r[i] * dim_);
    if( iuse < dim_ ) v.push_back(iuse);
  }
  delete [] r;

  // exit
  return (v.size()==npts);
}


bool SprIntegerBootstrap::replica(std::set<unsigned>& v, int npts)
{
  // init
  if( npts <= 0 ) npts = nsample_;

  // reset input set
  v.clear();

  // make array of uniform random numbers on [0,1]
  double* r = new double [npts];
  generator_.sequence(r,npts);
  unsigned iuse = 0;
  for( int i=0;i<npts;i++ ) {
    iuse = unsigned(r[i] * dim_);
    if( iuse < dim_ ) v.insert(iuse);
  }
  delete [] r;

  // exit
  return !v.empty();
}



