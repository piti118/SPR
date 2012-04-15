//$Id: SprRandomNumber.cc,v 1.4 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprRandomNumber.hh"

#include <sys/time.h>

using namespace std;

int SprRandomNumber::timeSeed(int seed)
{
  if( seed < 0 ) {
    struct timeval tp;
    gettimeofday(&tp, 0);
    return tp.tv_usec;
  }
  return seed;
}

void SprRandomNumber::init(int seed)
{
  theRanluxEngine_.setSeed(SprRandomNumber::timeSeed(seed));
}

void SprRandomNumber::sequence(double* seq, int npts)
{
  // make array of uniform random numbers on [0,1]

  theRanluxEngine_.flatArray(npts, seq);
}
