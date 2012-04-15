// File and Version Information:
//      $Id: SprTrainedBagger.hh,v 1.7 2008-01-03 20:51:58 narsky Exp $
//
// Description:
//      Class SprTrainedBagger :
//         Classifies a new event using a simple majority vote of all 
//         registered subclassifiers.
//
//         Discrete=true forces each subclassifier to return either 0 or 1.
//         Discrete=false makes each subclassifier run response() method
//         which can return a continuous value, depending on what subclassifier
//         is used. Both continuous and discrete values of subclassifiers
//         used to form an algebraic sum which represents the response value
//         of the bagger.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2005              California Institute of Technology
//
//------------------------------------------------------------------------
 
#ifndef _SprTrainedBagger_HH
#define _SprTrainedBagger_HH

#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"

#include <vector>
#include <utility>
#include <string>
#include <iostream>


class SprTrainedBagger : public SprAbsTrainedClassifier
{
public:
  virtual ~SprTrainedBagger() { this->destroy(); }

  SprTrainedBagger(const std::vector<
		   std::pair<const SprAbsTrainedClassifier*,bool> >& 
		   trained, bool discrete=false);

  SprTrainedBagger(const SprTrainedBagger& other);

  SprTrainedBagger* clone() const {
    return new SprTrainedBagger(*this);
  }

  /*
    Classifier name.
  */
  std::string name() const { return "Bagger"; }

  /*
    Classifier response for a data point. 
    Works only for problems with two categories, e.g., signal and background.
  */
  double response(const std::vector<double>& v) const;

  /*
    Generate code.
  */
  bool generateCode(std::ostream& os) const;

  // Change normalization.
  void useStandard() {}
  void useNormalized() {}
  bool normalized() const { return true; }

  // print out
  void print(std::ostream& o) const;

  /*
    Local accessors.
  */
  const SprAbsTrainedClassifier* classifier(int i) const {
    if( i>=0 && i<trained_.size() )
      return trained_[i].first;
    return 0;
  }

  void classifierList(std::vector<const SprAbsTrainedClassifier*>& classifiers)
    const {
    classifiers.clear();
    classifiers.resize(trained_.size());
    for( int i=0;i<trained_.size();i++ )
      classifiers[i] = trained_[i].first;
  }

  unsigned nClassifiers() const { return trained_.size(); }

  /*
    Local modifiers.
  */
  void setContinuous() { discrete_ = false; }
  void setDiscrete() { discrete_ = true; }
  bool discrete() const { return discrete_; }

  // addition
  SprTrainedBagger& operator+=(const SprTrainedBagger& other);

  /*
    useNClassifiers() lets you use a smaller number of classifiers
    than were trained by the Bagger. If the user supplies zero, all
    classifiers are used. All classifiers are used by default.
  */
  void useNClassifiers(unsigned nClassifiers=0) {
    nUsedClassifiers_ = nClassifiers;
  }
  unsigned nUsedClassifiers() const { return nUsedClassifiers_; }

protected:
  friend const SprTrainedBagger operator+(const SprTrainedBagger& l,
					  const SprTrainedBagger& r);

  void destroy();

  std::vector<std::pair<const SprAbsTrainedClassifier*,bool> > trained_;
  bool discrete_;
  unsigned nUsedClassifiers_;
};

#endif

