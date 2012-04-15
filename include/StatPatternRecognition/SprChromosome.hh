//$Id: SprChromosome.hh,v 1.2 2008-05-08 19:57:43 narsky Exp $
// File and Version Information:
//
//
// Description:

//
// Environment:
//      Software developed at Caltech's Center for Advanced Computing Research
//
// Author List:
//      Julian Bunn                     Original author
//      Ilya Narsky                     massive cleanup
//
// Copyright Information:
//      Copyright (C) 2008              California Institute of Technology
//------------------------------------------------------------------------
 
#ifndef _SprChromosome_HH
#define _SprChromosome_HH

#include <string>
#include <iostream>
#include <vector>

class SprRandomNumber;


class SprGene 
{
public:
  virtual ~SprGene() {}

  SprGene() 
    : head_(), tail_() {}

  SprGene(const SprGene& other)
    :
    head_(other.head_),
    tail_(other.tail_)
  {}

  SprGene& operator=(const SprGene& other) {
    head_ = other.head_;
    tail_ = other.tail_;
    return *this;
  }

  void setHead(const std::vector<int>& head) { head_ = head; }
  void setTail(const std::vector<int>& tail) { tail_ = tail; }

  std::vector<int> getHead() { return head_; }
  std::vector<int> getTail() { return tail_; }

  void generateTree(std::vector<std::vector<int> >& tree) const;
  
  void print(std::ostream& os) const;

private:
  std::vector<int> head_;
  std::vector<int> tail_;
};


class SprChromosome 
{
public:
  unsigned num_genes_;
  char link_function_;
  std::vector<SprGene> gene_;
  std::vector<double> constants_;
  double fitness_;
  unsigned length_;
  double shift_;
		
  virtual ~SprChromosome();

  SprChromosome();
  SprChromosome(SprRandomNumber* rndm);
  SprChromosome(SprRandomNumber* rndm,
		unsigned num_genes,
		char link_function);
  
  double Evaluate(const std::vector<double> &data, int verbose) const;
  void generateCode();
  
  bool mutate(const std::vector<int>& funcs, 
	      unsigned dim, 
	      unsigned constants_per_chromosome);
  bool RIS();
  bool IS();
  bool OnePoint(SprChromosome &other);
  bool TwoPoint(SprChromosome &other);
  bool WholeGene(SprChromosome &other);
  bool ConstantMutation(double range);
  bool ConstantSwap(SprChromosome &other);

  void setLength(unsigned num_genes, unsigned num_constants, unsigned length) {
    num_genes_ = num_genes;
    gene_.resize(num_genes);
    constants_.resize(num_constants);
    length_ = length;
  }

  void print(std::ostream& os) const;

private:
  SprRandomNumber* rndm_;
  bool ownRndm_;
};

#endif
