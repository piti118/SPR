//------------------------------------------------------------------------
// File and Version Information:
//      $Id: SprCrossValidator.hh,v 1.7 2008-04-02 23:36:44 narsky Exp $
//
// Description:
//      Class SprCrossValidator :
//         Cross-validates data by dividing the data into a specified
//         number of equal-sized pieces, training classifiers on all
//         pieces but one and computing FOM for the leftover piece.
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
 
#ifndef _SprCrossValidator_HH
#define _SprCrossValidator_HH

#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <vector>
#include <map>
#include <cassert>

class SprAbsFilter;
class SprAbsTwoClassCriterion;
class SprAverageLoss;
class SprAbsClassifier;
class SprClass;
class SprAbsMultiClassLearner;


class SprCrossValidator
{
public:
  virtual ~SprCrossValidator();

  SprCrossValidator(const SprAbsFilter* data, unsigned nPieces);

  /*
    Returns a vector of cross-validated FOMs for the vector
    of supplied classifiers. Note that, due to some technicalities, this method
    works under assumption that the content of the supplied data has not
    changed since the time SprCrossValidator was constructed.

    If "integrate" flag is set to true, FOM is integrated over all
    points in the supplied data. If "integrate" is false, a separate
    FOM is computed for each validation piece and then an average of
    these FOMs is estimated and an error is derived from the spread of
    the FOM values around the average. If "integrate" is true, no
    error is estimated.
  */
  bool validate(const SprAbsTwoClassCriterion* crit,
		SprAverageLoss* loss,
		const std::vector<SprAbsClassifier*>& classifiers,
		const SprClass& cls0, const SprClass& cls1,
		std::vector<SprValueWithError>& crossFom,
		bool integrate,
		int verbose=0) const;

  /*
    Returns the loss and classification table for the multiclass learner.
  */
  bool validate(SprAbsMultiClassLearner* mcLearner,
		SprAverageLoss* loss,
		SprValueWithError& realLoss,
		SprClassificationTable& classificationTable,
		std::map<int,double>& weightInClass,
		int verbose) const;

private:
  // methods
  bool divide(unsigned nPieces);

  // data
  const SprAbsFilter* data_;
  std::vector<const SprAbsFilter*> samples_;
};

#endif

