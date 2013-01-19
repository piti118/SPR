// File and Version Information:
//      $Id: SprTrainedMultiClassLearner.hh,v 1.7 2008-04-02 23:36:44 narsky Exp $
//
// Description:
//      Class SprTrainedMultiClassLearner :
//          Trained classifier implementing the Allwein-Schapire-Singer
//            method.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2005-2008         California Institute of Technology
//------------------------------------------------------------------------
 
#ifndef _SprTrainedMultiClassLearner_HH
#define _SprTrainedMultiClassLearner_HH

#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprMatrix.hh"

#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <string>

class SprAbsTrainedClassifier;


class SprTrainedMultiClassLearner : public SprAbsTrainedMultiClassLearner
{
public:
  typedef double (*SprPerEventLoss)(int,double);
  typedef double (*Spr1DTransformer)(double);

  virtual ~SprTrainedMultiClassLearner() { this->destroy(); }

  SprTrainedMultiClassLearner(const SprMatrix& indicator,
			      const std::vector<int>& mapper,
			      const std::vector<std::pair<
			      const SprAbsTrainedClassifier*,bool> >& 
			      classifiers);

  SprTrainedMultiClassLearner(const SprTrainedMultiClassLearner& other);

  /*
    Returns classifier name.
  */
  std::string name() const { return "MultiClassLearner"; }

  /*
    Make a clone.
  */
  SprTrainedMultiClassLearner* clone() const {
    return new SprTrainedMultiClassLearner(*this);
  }

  /*
    Set weights for classifier normalization. By default all weights are 1.
  */
  void setClassifierWeights(const std::vector<double>& weights);

  /*
    Include non-participating classes in loss evaluation?
  */
  void excludeZeroClasses() { includeZeroIndicator_ = false; }
  void includeZeroClasses() { includeZeroIndicator_ = true; }

  /*
    Set appropriate loss and transformation.
  */
  void setLoss(SprPerEventLoss loss, Spr1DTransformer trans=0) {
    loss_ = loss; trans_ = trans;
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
  void printIndicatorMatrix(std::ostream& os) const;

protected:
  int response_one(const std::vector<double>& input,
		   std::map<int,double>& output) const;

private:
  void destroy();

  SprMatrix indicator_;
  std::vector<std::pair<const SprAbsTrainedClassifier*,bool> > classifiers_;
  std::vector<double> weights_;// classifier weights for normalization
  bool includeZeroIndicator_;
  SprPerEventLoss loss_;
  Spr1DTransformer trans_;
};

#endif
