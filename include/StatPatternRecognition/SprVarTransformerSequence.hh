// File and Version Information:
//      $Id: SprVarTransformerSequence.hh,v 1.2 2008-02-08 23:22:09 narsky Exp $
//
// Description:
//      Class SprVarTransformerSequence :
//          Transformer sequence.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2007              California Institute of Technology
//------------------------------------------------------------------------
 
#ifndef _SprVarTransformerSequence_HH
#define _SprVarTransformerSequence_HH

#include "StatPatternRecognition/SprAbsVarTransformer.hh"

#include <vector>
#include <utility>
#include <iostream>

class SprAbsFilter;


class SprVarTransformerSequence : public SprAbsVarTransformer
{
public:
  virtual ~SprVarTransformerSequence();

  SprVarTransformerSequence();
  SprVarTransformerSequence(const std::vector<
			    std::pair<SprAbsVarTransformer*,bool> >& 
			    transformers);
  SprVarTransformerSequence(const SprVarTransformerSequence& other);

  // Clone
  SprVarTransformerSequence* clone() const {
    return new SprVarTransformerSequence(*this);
  }

  /*
    VarTransformer name.
  */
  std::string name() const { return "TransformerSequence"; }

  /*
    Computes transformation using input data. 
    Returns true on success, false otherwise.
  */
  bool train(const SprAbsFilter* data, int verbose=0);

  // Applies transformation.
  void transform(const std::vector<double>& in,
		 std::vector<double>& out) const;

  // Applies inverse transformation.
  void inverse(const std::vector<double>& in,
	       std::vector<double>& out) const;

  // Status of the transformer - if returns true, ready to transform.
  bool ready() const;

  // Returns true if all transformations for individual variables
  // are independent of each other.
  bool allVarsIndependent() const;

  // Eliminates transformation variables which are not in the supplied list.
  bool reduceVars(const std::vector<std::string>& vars);

  // Output
  void print(std::ostream& os) const;

private:
  bool initVars();

  std::vector<std::pair<SprAbsVarTransformer*,bool> > transformers_;
};

#endif
