//$Id: SprTrainedGEP.hh,v 1.3 2008-05-08 19:57:43 narsky Exp $
// File and Version Information:
//
// Description:
//      Class SprTrainedGEP :
//          Interface for trained classifiers.
//          The purpose of this class is to generate response of 
//          a trained classifier on validation or test data.
//
// Environment:
//      Software developed at Caltech's Center for Advanced Computing Research
//
// Author List:
//      Julian Bunn                     Original author
//
// Copyright Information:
//      Copyright (C) 2008              California Institute of Technology
//
//------------------------------------------------------------------------
 
#ifndef _SprTrainedGEP_HH
#define _SprTrainedGEP_HH

#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprChromosome.hh"

#include <iostream>
#include <string>
#include <vector>


class SprTrainedGEP : public SprAbsTrainedClassifier
{
public:
  virtual ~SprTrainedGEP() {}

  SprTrainedGEP(); 

  SprTrainedGEP(const SprTrainedGEP& other);

  SprTrainedGEP(const SprChromosome& best);


  /*
    Returns classifier name.
  */
  std::string name() const {
    return "GEP";
  }

  /*
    Make a clone.
  */
  SprTrainedGEP* clone() const {
    return new SprTrainedGEP(*this);
  }

  /*
    Classifier response for a data point. 
    Binary split produces binary response only: 
    0 for background and 1 for signal.
  */
  double response(const std::vector<double>& v) const;

  bool generateCode(std::ostream& os) const { 
    return false; 
  } 

  // Change normalization.
  void useStandard() {}
  void useNormalized() {}
  bool normalized() const { return false; }

  /*
    Print out.
  */
  void print(std::ostream& os) const;

  // Local methods.

private:
  SprChromosome chromosome_;
};

#endif
