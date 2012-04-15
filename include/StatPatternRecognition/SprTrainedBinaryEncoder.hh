// File and Version Information:
//      $Id: SprTrainedBinaryEncoder.hh,v 1.1 2008-04-02 23:36:44 narsky Exp $
//
// Description:
//      Class SprTrainedBinaryEncoder :
//          Interface for trained multiclass methods.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2008              California Institute of Technology
//------------------------------------------------------------------------
 
#ifndef _SprTrainedBinaryEncoder_HH
#define _SprTrainedBinaryEncoder_HH

#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"

#include <vector>
#include <map>
#include <iostream>
#include <string>

class SprAbsTrainedClassifier;


class SprTrainedBinaryEncoder : public SprAbsTrainedMultiClassLearner
{
public:
  virtual ~SprTrainedBinaryEncoder() { this->destroy(); }

  SprTrainedBinaryEncoder(const std::vector<int>& mapper,
			  const SprAbsTrainedClassifier* classifier,
			  bool ownClassifier=false);

  SprTrainedBinaryEncoder(const SprTrainedBinaryEncoder& other);

  /*
    Returns classifier name.
  */
  std::string name() const { return "BinaryEncoder"; }

  /*
    Make a clone.
  */
  SprTrainedBinaryEncoder* clone() const {
    return new SprTrainedBinaryEncoder(*this);
  }

  /*
    Generate code.
  */
  bool generateCode(std::ostream& os) const {
    return false;
  }

  /*
    Print out.
  */
  void print(std::ostream& os) const;

protected:
  int response_one(const std::vector<double>& input,
		   std::map<int,double>& output) const;

private:
  void destroy();

  const SprAbsTrainedClassifier* classifier_;
  bool ownClassifier_;
};

#endif
