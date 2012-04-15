//$Id: SprTrainedGEP.cc,v 1.2 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprTrainedGEP.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprUtils.hh"

using namespace std;


SprTrainedGEP::SprTrainedGEP() 
  :  
  SprAbsTrainedClassifier(),
  chromosome_()
{
  this->setCut(SprUtils::lowerBound(0.));
}


SprTrainedGEP::SprTrainedGEP(const SprTrainedGEP& other)
  : 
  SprAbsTrainedClassifier(other),
  chromosome_(other.chromosome_)
{}


SprTrainedGEP::SprTrainedGEP(const SprChromosome& best) 
  :  
  SprAbsTrainedClassifier(),
  chromosome_(best)
{
  this->setCut(SprUtils::lowerBound(0.));
}


double SprTrainedGEP::response(const std::vector<double>& v) const
{
  // evaluate the best chrom. with the vector v
  return chromosome_.Evaluate(v,0);
}


void SprTrainedGEP::print(std::ostream& os) const
{
  os << "Trained GEP " << SprVersion << endl;
  chromosome_.print(os);
}

