// File and Version Information:
//      $Id: SprTreeNode.hh,v 1.9 2007-08-11 22:08:10 narsky Exp $
//
// Description:
//      Class SprTreeNode :
//          Implements a node of the decision tree.
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
 
#ifndef _SprTreeNode_HH
#define _SprTreeNode_HH

#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprBoxFilter.hh"

#include <vector>
#include <utility>

class SprAbsFilter;
class SprAbsTwoClassCriterion;
class SprIntegerBootstrap;
class SprTrainedNode;


class SprTreeNode
{
public:
  virtual ~SprTreeNode();

  SprTrainedNode* makeTrained() const;

  double fom() const { return fom_; }

  double w0() const { return w0_; }
  double w1() const { return w1_; }

  unsigned n0() const { return n0_; }
  unsigned n1() const { return n1_; }

  void box(SprBox& limits) const {
    limits = limits_;
  }

  SprInterval limits(int d) const;

  int id() const { return id_; }

  int nodeClass() const { return nodeClass_; }

private:
  friend class SprDecisionTree;
  friend class SprTopdownTree;

  SprTreeNode(const SprAbsTwoClassCriterion* crit,
	      const SprAbsFilter* data,
	      bool allLeafsSignal,
	      int nmin,
	      bool discrete,
	      bool canHavePureNodes,
	      bool fastSort,
	      SprIntegerBootstrap* bootstrap=0);

  SprTreeNode(const SprAbsTwoClassCriterion* crit,
	      const SprBoxFilter& data,
	      bool allLeafsSignal,
	      int nmin,
	      bool discrete,
	      bool canHavePureNodes,
	      bool fastSort,
	      const SprClass& cls0,
	      const SprClass& cls1,
	      const SprTreeNode* parent,
	      const SprBox& limits,
	      SprIntegerBootstrap* bootstrap=0);

  bool split(std::vector<SprTreeNode*>& nodesToSplit, 
	     std::vector<std::pair<int,double> >& countTreeSplits,
	     int verbose=0);

  bool setClasses(const SprClass& cls0, const SprClass& cls1);

  bool sort(unsigned d, std::vector<int>& sorted,
	    std::vector<double>& division);

  bool prepareExit(bool status);

  const SprAbsTwoClassCriterion* crit_;
  SprBoxFilter* data_;
  bool allLeafsSignal_;
  int nmin_;// minimal number of events per node
  bool discrete_;// type of node output: discrete (0 or 1) or continuous
  bool canHavePureNodes_;// true if allow nodes w/ only signal or bgrnd events
  bool fastSort_;
  SprClass cls0_;
  SprClass cls1_;
  const SprTreeNode* parent_;
  SprTreeNode* left_;
  SprTreeNode* right_;
  double fom_;
  double deltaFom_;
  double w0_;
  double w1_;
  unsigned n0_;
  unsigned n1_;
  SprBox limits_;
  int id_;
  int nodeClass_;
  int d_;
  double cut_;
  SprIntegerBootstrap* bootstrap_;

  static int counter_;
};

#endif
