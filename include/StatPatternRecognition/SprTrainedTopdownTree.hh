// File and Version Information:
//      $Id: SprTrainedTopdownTree.hh,v 1.4 2007-07-11 19:52:10 narsky Exp $
//
// Description:
//      Class SprTrainedTopdownTree :
//         Implements a trained decision tree. The constructor input
//         is a vector of all (terminal and non-terminal) nodes of the tree. 
//         By convention, the 0th element in the vector has to be the root
//         node with id=0.
//
//         Compared to SprTrainedDecisionTree, this implementation offers
//         a faster response algorithm, however the configuration file
//         is less interpretable. The tree is now represented by the full
//         path from top to bottom. The response time is reduced by a factor
//         of log(N)/N where N is the number of terminal nodes in the tree.
//         For large trees typically produced by the bagging algorithm,
//         this is is a substantial improvement. Because the tree is no longer
//         saved as a list of boxes, it is harder for the user to interpret
//         the classifier configuration file.
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
 
#ifndef _SprTrainedTopdownTree_HH
#define _SprTrainedTopdownTree_HH

#include "StatPatternRecognition/SprTrainedDecisionTree.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprDefs.hh"

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

class SprTrainedNode;


class SprTrainedTopdownTree : public SprTrainedDecisionTree
{
public:
  virtual ~SprTrainedTopdownTree();

  SprTrainedTopdownTree(const std::vector<const SprTrainedNode*>& nodes, 
			bool ownTree=false)
    :
    SprTrainedDecisionTree(std::vector<SprBox>()),
    nodes_(nodes),
    ownTree_(ownTree)
  {
    assert( !nodes_.empty() );
    this->setCut(SprUtils::lowerBound(0.5));
  }

  SprTrainedTopdownTree(const SprTrainedTopdownTree& other)
    :
    SprTrainedDecisionTree(other),
    nodes_(),
    ownTree_(true)
  {
    bool status = this->replicate(other.nodes_);
    assert( status );
  }

  /*
    Returns classifier name.
  */
  std::string name() const {
    return "TopdownTree";
  }

  /*
    Make a clone.
  */
  SprTrainedTopdownTree* clone() const {
    return new SprTrainedTopdownTree(*this);
  }

  /*
    Classifier response for a data point. 
  */
  double response(const std::vector<double>& v) const;

  /*
    Generate code.
  */
  bool generateCode(std::ostream& os) const {
    this->printFunction(os,0);
    return true;
  }

  /*
    Print out.
  */
  void print(std::ostream& os) const;

  // Local methods.

protected:
  bool replicate(const std::vector<const SprTrainedNode*>& nodes);
  void printFunction(std::ostream& os,
		     const SprTrainedNode* currentNode,
		     int indentLevel=2) const;

  std::vector<const SprTrainedNode*> nodes_;
  bool ownTree_;
};

#endif
