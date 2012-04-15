//$Id: SprTrainedStdBackprop.hh,v 1.7 2008-01-18 22:56:59 narsky Exp $
//--------------------------------------------------------------------------
//
//      Class SprTrainedStdBackprop:
//        Output of a trained SNNS network with standard backpropagation
//        algorithm and either identity or logistic (aka sigmoid) 
//        activation function.
//
// Copied from Harvard folks and modified for StatPatternRecognition 
//    compliance    - Ilya Narsky, 2005
//
// SNNSNetwork.hh
//
// 29 June 2002, Masahiro Morii, Harvard Univeristy
//
// Reads SNNS .net files and calculates output from inputs.
// Based on InfoSimpleSNNSNetwork by J. Beringer.
//
//------------------------------------------------------------------------
#ifndef SprTrainedStdBackprop_HH
#define SprTrainedStdBackprop_HH

#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprNNDefs.hh"

#include <vector>
#include <iostream>
#include <string>
#include <utility>


class SprTrainedStdBackprop : public SprAbsTrainedClassifier
{
public:
  virtual ~SprTrainedStdBackprop() {}

  SprTrainedStdBackprop(); 

  SprTrainedStdBackprop(const char* structure,
			const std::vector<SprNNDefs::NodeType>& nodeType,
			const std::vector<SprNNDefs::ActFun>& nodeActFun,
			const std::vector<int>& nodeNInputLinks,
			const std::vector<int>& nodeFirstInputLink,
			const std::vector<int>& linkSource,
			const std::vector<double>& nodeBias,
			const std::vector<double>& linkWeight);

  SprTrainedStdBackprop(const SprTrainedStdBackprop& other);

  SprTrainedStdBackprop* clone() const {
    return new SprTrainedStdBackprop(*this);
  }

  std::string name() const { return "StdBackprop"; }

  double response(const std::vector<double>& v) const;

  bool generateCode(std::ostream& os) const {
    return false;
  }

  void useStandard() {}
  void useNormalized() {}
  bool normalized() const { return true; }

  void print(std::ostream& os) const;

  bool varImportance(std::vector<std::pair<std::string,double> >& 
		     importance) const;

private:
  double activate(double x, SprNNDefs::ActFun f) const;

  int     nNodes_; // Total number of nodes (input+hidden+output)
  int     nLinks_; // Total number of connections between nodes

  std::string structure_;
  std::vector<SprNNDefs::NodeType>   nodeType_;
  std::vector<SprNNDefs::ActFun>     nodeActFun_;
  std::vector<int>                   nodeNInputLinks_;
  std::vector<int>                   nodeFirstInputLink_;
  std::vector<int>                   linkSource_;
  std::vector<double>                nodeBias_;
  std::vector<double>                linkWeight_;
};

#endif
