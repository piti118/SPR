// File and Version Information:
//      $Id: SprClassifierEvaluator.hh,v 1.7 2008-04-02 23:36:44 narsky Exp $
//
// Description:
//      Class SprClassifierEvaluator :
//          Evaluates various quantities for trained classifiers.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2007              California Institute of Technology
//
//------------------------------------------------------------------------
 
#ifndef _SprClassifierEvaluator_HH
#define _SprClassifierEvaluator_HH

#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <vector>
#include <utility>
#include <string>
#include <set>

class SprAbsFilter;
class SprAbsClassifier;
class SprAbsTrainedClassifier;
class SprAbsMultiClassLearner;
class SprAbsTrainedMultiClassLearner;
class SprCoordinateMapper;
class SprAbsTwoClassCriterion;
class SprAverageLoss;


struct SprClassifierEvaluator
{
  typedef std::pair<std::string,SprValueWithError> NameAndValue;
  typedef std::pair<std::set<std::string>,SprValueWithError> SetVarsAndValue;
  typedef 
  std::pair<std::vector<std::string>,SprValueWithError> ListVarsAndValue;

  /*
    Computes variable importance by randomly permuting class labels
    in the data and calculating the increase in the quadratic loss
    for two-class classifiers or misid rate for the multi-class learner.
    Returns a vector of loss increases with statistical errors.
  */
  static bool variableImportance(const SprAbsFilter* data,
				 SprAbsTrainedClassifier* trained,
				 SprAbsTrainedMultiClassLearner* mcTrained,
				 SprAverageLoss* loss,
				 SprCoordinateMapper* mapper,
				 unsigned nPerm,
				 std::vector<NameAndValue>& lossIncrease);

  /*
    Computes interaction between two specified subsets of variables.
  */
  static bool twoSubsetInteraction(const SprAbsFilter* data,
				  SprAbsTrainedClassifier* trained,
				  SprAbsTrainedMultiClassLearner* mcTrained,
				  SprCoordinateMapper* mapper,
				  const std::vector<std::string>& subset1,
				  const std::vector<std::string>& subset2,
				  unsigned nPoints,
				  SprValueWithError& interaction,
				  int verbose=0);

  /*
    Computes interaction between two subsets of variables. The syntax
    for defining subsets of variables is the same as the syntax
    for defining classes: Groups of variables are separated by ':' and
    variables within on egroup are separated by ','. If a group is
    specified as '.', this means all variables but those in the other group,
    for example, 'x1:.' means that the method will compute interaction
    between x1 and all variables excluding x1. If an empty string
    is entered, the method computes interaction between each variable
    known to the trained classifier and all other variables.

    Interaction is defined as Correlation(F(S1),F(S2)), where

        F(S1) is the classifier response at a given point integrated
	      over all variables not included in subset S1
  */
  static bool variableInteraction(const SprAbsFilter* data,
				  SprAbsTrainedClassifier* trained,
				  SprAbsTrainedMultiClassLearner* mcTrained,
				  SprCoordinateMapper* mapper,
				  const char* vars,
				  unsigned nPoints,
				  std::vector<NameAndValue>& interaction,
				  int verbose=0);

  /*
    Selects the most interacting variable sequentially: first selects
    the variable that interacts most with all other variables, then
    selects the 2nd variable that interacts most with all other
    variables except the variable selected at step 1, then selects 3rd
    variable that interacts most with all other variables except
    variables 1 and 2 etc.  The outcome is a list of variable subsets
    with corresponding interactions. The variable subsets in this list
    are ordered by size: the 1st subset consists of one variable, the
    2nd subset consists of 2 variables etc, with each subset being the
    most powerful subset among the subsets of the same size based on
    the strength of interaction.
  */
  static bool sortByInteraction(const SprAbsFilter* data,
				SprAbsTrainedClassifier* trained,
				SprAbsTrainedMultiClassLearner* mcTrained,
				SprCoordinateMapper* mapper,
				unsigned nPoints,
				std::vector<ListVarsAndValue>& interaction,
				int verbose=0);

  /*
    Builds models using the supplied trainable classifier on subsets
    of input variables. Starts with the null set, adds n variables
    sequentially one at a time, then removes r variables sequentially
    one at a time, at each step inserting the most efficient or
    removing the least efficient variable, as specfied by the chosen
    loss. If no test data are supplied, the models are averaged by
    using nCross pieces of the training data for
    cross-validation. Otherwise the models are always built on
    trainData and tested on testData. The method returns a vector of
    length equal to the dimensionality of the supplied data. Each
    element of the vector shows the most efficient combination of
    variables and associated loss. The two loss types recommended for
    this method are: quadratic and correct_id.

    See SprCrossValidator.hh for description of the "integrate" flag.
  */
  static bool addNremoveR(const SprAbsFilter* trainData,
			  const SprAbsFilter* testData,
			  SprAbsClassifier* trainable,
			  SprAbsMultiClassLearner* mcTrainable,
			  unsigned N, unsigned R,
			  SprAverageLoss* loss,
			  unsigned nCross,
			  std::vector<SetVarsAndValue>& vars_and_loss,
			  bool integrate,
			  int verbose=0);

private:
  static std::pair<SprValueWithError,bool> lossPerVar(
				    const std::set<std::string>& attempt_model,
				    const SprAbsFilter* trainData,
				    const SprAbsFilter* testData,
				    unsigned nCross,
				    SprAbsClassifier* trainable,
				    SprAbsMultiClassLearner* mcTrainable,
				    SprAverageLoss* aveLoss,
				    const std::vector<SprClass>& classes,
				    bool integrate,
				    int verbose);
};

#endif
