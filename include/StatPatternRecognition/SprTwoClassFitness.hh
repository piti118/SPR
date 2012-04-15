//$Id: SprTwoClassFitness.hh,v 1.2 2008-04-02 23:36:44 narsky Exp $
// File and Version Information:
//
// Description:
//      Class SprTwoClassFitness :
//        Returns the Fitness = Specificity*Sensitivity 
//		  where Sensitivity = wcor1/(wcor1+wmis1);
//        and Specificity = wcor0/(wcor0+wmis0);
//
// Environment:
//
// Author List:
//      Julian Bunn                     Original author
//
// Copyright Information:
//      Copyright (C) 2008              California Institute of Technology
//
//------------------------------------------------------------------------
 
#ifndef _SprTwoClassFitness_HH
#define _SprTwoClassFitness_HH

#include "StatPatternRecognition/SprAbsTwoClassCriterion.hh"
#include "StatPatternRecognition/SprUtils.hh"

#include <cmath>
#include <iostream>


class SprTwoClassFitness : public SprAbsTwoClassCriterion
{
public:
  virtual ~SprTwoClassFitness() {}

  SprTwoClassFitness() : SprAbsTwoClassCriterion() {}

  double fom(double wcor0, double wmis0, double wcor1, double wmis1) const {
    if( (wcor1+wmis1) < SprUtils::eps() ) return this->min();
    if( (wcor0+wmis0) < SprUtils::eps() ) return this->min();
    
    double sensitivity = wcor1/(wcor1+wmis1);
    double specificity = wcor0/(wcor0+wmis0);
    return sensitivity*specificity;
  }

  bool symmetric() const { return true; }

  double min() const { return -1; }
  double max() const { return  1; }

  double dfom_dwmis0(double wcor0, double wmis0, 
		     double wcor1, double wmis1) const {
    std::cerr << "Derivative for fitness not implemented." << std::endl; 
    return 0;
  }

  double dfom_dwcor1(double wcor0, double wmis0, 
		     double wcor1, double wmis1) const {
    std::cerr << "Derivative for fitness not implemented." << std::endl; 
    return 0;
  }
};

#endif
