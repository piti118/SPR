//------------------------------------------------------------------------
// File and Version Information:
//      $Id: SprFomCalculator.hh,v 1.5 2008-04-02 23:36:44 narsky Exp $
//
// Description:
//      Class SprFomCalculator :
//         Computes FOM for specified data and classifier.
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
 
#ifndef _SprFomCalculator_HH
#define _SprFomCalculator_HH

#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprPlotter.hh"
#include "StatPatternRecognition/SprMultiClassPlotter.hh"

#include <map>
#include <vector>

class SprAbsFilter;
class SprAbsTwoClassCriterion;
class SprAverageLoss;
class SprAbsTrainedClassifier;
class SprClass;
class SprAbsTrainedMultiClassLearner;


class SprFomCalculator
{
public:
  virtual ~SprFomCalculator() {}

  SprFomCalculator() {}

  static SprValueWithError fom(
	const std::vector<std::vector<SprPlotter::Response> >& responses,
	const SprAbsTwoClassCriterion* crit, 
	SprAverageLoss* loss, 
	bool integrate, int verbose=0);

  static SprValueWithError fom(
		    const std::vector<const SprAbsFilter*>& data, 
		    const std::vector<const SprAbsTrainedClassifier*>& trained,
		    const SprAbsTwoClassCriterion* crit, 
		    SprAverageLoss* loss, 
		    const SprClass& cls0, const SprClass& cls1,
		    bool integrate, int verbose=0);

  static SprValueWithError loss(
              const std::vector<int>& classes,
              const std::vector<SprMultiClassPlotter::Response>& responses,
	      SprAverageLoss* loss,
	      SprClassificationTable& classificationTable, 
	      std::map<int,double>& weightInClass,
	      int verbose=0);

  static SprValueWithError loss(
	     const std::vector<const SprAbsFilter*>& data, 
	     const std::vector<const SprAbsTrainedMultiClassLearner*>& trained,
	     SprAverageLoss* loss,
	     SprClassificationTable& classificationTable, 
	     std::map<int,double>& weightInClass,
	     int verbose=0);
};

#endif

