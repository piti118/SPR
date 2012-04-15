// File and Version Information:
//      $Id: SprTrainedRBF.hh,v 1.7 2008-01-03 20:51:58 narsky Exp $
//
// Description:
//      Class SprTrainedRBF :
//          Returns response of a trained Radial Basis Function
//          SNNS network.
/*
  For description of radial basis function methods see
     Haykin "Neural Networks, a Comprehensive Foundation"
*/
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
#ifndef SprTrainedRBF_HH
#define SprTrainedRBF_HH

#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprNNDefs.hh"

#include <vector>
#include <iostream>
#include <string>


class SprTrainedRBF : public SprAbsTrainedClassifier
{
public:
  enum ActRBF { Gauss=1, MultiQ, ThinPlate };

  class Link;
  class Node;

  struct Link {
    ~Link() {}
    Link() : source_(0), target_(0), weight_(0) {}
    Link(const Link& other) 
      :
      source_(other.source_),
      target_(other.target_),
      weight_(other.weight_) 
    {}

    Node* source_;
    Node* target_;
    double weight_;
  };

  struct Node {
    ~Node() {}

    Node() : 
      index_(0), type_(), actFun_(), actRBF_(), outFun_(), act_(0), bias_(0), 
      incoming_(), outgoing_() {}

    Node(const Node& other) :
      index_(other.index_),
      type_(other.type_),
      actFun_(other.actFun_),
      actRBF_(other.actRBF_),
      outFun_(other.outFun_),
      act_(other.act_),
      bias_(other.bias_),
      incoming_(other.incoming_),
      outgoing_(other.outgoing_)
    {}

    unsigned index_;
    SprNNDefs::NodeType type_;
    SprNNDefs::ActFun actFun_;
    ActRBF actRBF_;
    SprNNDefs::OutFun outFun_;
    double act_;
    double bias_;
    std::vector<Link*> incoming_;
    std::vector<Link*> outgoing_;
  };

  virtual ~SprTrainedRBF() { this->destroy(); }

  SprTrainedRBF()
    : SprAbsTrainedClassifier(),
      initialized_(false), nodes_(), links_()
  {}

  SprTrainedRBF(const SprTrainedRBF& other) 
    : 
    SprAbsTrainedClassifier(other), 
    initialized_(other.initialized_), 
    nodes_(),
    links_()
  {
    this->correspondence(other);
  }

  SprTrainedRBF* clone() const {
    return new SprTrainedRBF(*this);
  }

  bool readNet(const char* netfile);

  void printNet(std::ostream& os) const;
  void print(std::ostream& os) const { this->printNet(os); }

  std::string name() const { return "RBF"; }

  double response(const std::vector<double>& v) const;

  bool generateCode(std::ostream& os) const {
    return false;
  }

  void useStandard() {}
  void useNormalized() {}
  bool normalized() const { return true; }

private:
  void destroy();
  void correspondence(const SprTrainedRBF& other);
  double rbf(double r2, double p, ActRBF act) const;// RBF function
  double act(double x, double p, SprNNDefs::ActFun act) const;

  bool initialized_;
  std::vector<Node*> nodes_;
  std::vector<Link*> links_;
};

#endif
