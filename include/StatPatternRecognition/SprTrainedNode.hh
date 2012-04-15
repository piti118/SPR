// File and Version Information:
//      $Id: SprTrainedNode.hh,v 1.4 2007-02-05 21:49:45 narsky Exp $
//
// Description:
//      Class SprTrainedNode :
//          Keeps info about a trained node of the decision tree.
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
 
#ifndef _SprTrainedNode_HH
#define _SprTrainedNode_HH

struct SprTrainedNode
{
  virtual ~SprTrainedNode() {}

  SprTrainedNode() 
    : id_(-1), score_(0.5), d_(-1), cut_(0),
      toDau1_(0), toDau2_(0), toParent_(0),
      n0_(0),n1_(0),w0_(0.0),w1_(0.0),
      deltaFom_(0.0)
  {}

  SprTrainedNode(const SprTrainedNode& other)
    : 
    id_(other.id_),
    score_(other.score_),
    d_(other.d_),
    cut_(other.cut_),
    toDau1_(0),
    toDau2_(0),
    n0_(other.n0_),n1_(other.n1_),w0_(other.w0_),w1_(other.w1_),
    deltaFom_(0.0),
    toParent_(0)
  {}

  int id_;// id of this node
  double score_;// score function for this node; must vary from 0 to 1
  int d_;// dimension on which the node is split
  double cut_;// the split
  int n0_,n1_;//raw number of how many of class 1 or zero
  double w0_,w1_;//sum of weights of class 1 or zero
  const SprTrainedNode* toDau1_;
  const SprTrainedNode* toDau2_;
  const SprTrainedNode* toParent_;
  double deltaFom_;
};

#endif
