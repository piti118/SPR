//$Id: SprGEP.cc,v 1.3 2008-05-09 21:25:26 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprGEP.hh"
#include "StatPatternRecognition/SprChromosome.hh"
#include "StatPatternRecognition/SprTrainedGEP.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprAbsTwoClassCriterion.hh"
#include "StatPatternRecognition/SprAverageLoss.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprTwoClassIDFraction.hh"

#include <stdio.h>
#include <cassert>
#include <algorithm>
#include <functional>
#include <cmath>
#include <sstream>

using namespace std;

const int SprGEP::Variable = 500;
const int SprGEP::Constant = 400; 
const int SprGEP::Max_Depth = 40;


struct SGEPCmpPairDDFirst 
  : public binary_function<pair<double,double>,pair<double,double>,bool> {
  bool operator()(const pair<double,double>& l, const pair<double,double>& r)
    const {
    return (l.first < r.first);
  }
};


struct SGEPCmpPairDDFirstNumber
  : public binary_function<pair<double,double>,double,bool> {
  bool operator()(const pair<double,double>& l, double r) const {
    return (l.first < r);
  }
};


SprGEP::SprGEP(SprAbsFilter* data, 
	       const SprAbsTwoClassCriterion* crit,
	       unsigned gene_head_length,
	       unsigned genes_per_chromosome,
	       unsigned constants_per_chromosome,
	       unsigned population_size,
	       int Max_Epoch,
	       const string& functions,
	       int seed)
  : 
  SprAbsClassifier(data),
  rndm_(seed),
  crit_(crit),
  gene_head_length_(gene_head_length),
  genes_per_chromosome_(genes_per_chromosome),
  constants_per_chromosome_(constants_per_chromosome),
  population_size_(population_size),
  Max_Epoch_(Max_Epoch),
  function_arg_max_(0),
  chromosomes_(population_size,SprChromosome(&rndm_)),
  mutation_rate_(0.05), // 0.05
  RIS_rate_(0.1), // 0.1
  IS_rate_(0.1), // 0.1
  One_Point_rate_(0.3), // 0.3
  Two_Point_rate_(0.3), // 0.3
  Whole_Gene_rate_(0.1), // 0.1
  Constant_mutation_(0.1), // 0.1 0.7
  Constant_swap_(0.1), // 0.1
  Constant_range_(10.0), // 10.0
  fractionOfHeadsForFuncts_(0.8), // 0.8
  best_(&rndm_),
  cls0_(0), 
  cls1_(1), 
  valData_(0),
  valPrint_(0)
{
  assert (gene_head_length > 1);
  assert (genes_per_chromosome > 0);
  assert (population_size > 1);

  cout << "Begin GEP" << endl 
       << " Genes per Chromosome         = " 
       << genes_per_chromosome_ << endl 
       << " Gene head length             = " <<	gene_head_length_ << endl 
       << " Constants per Chromosome     = " 
       << constants_per_chromosome_ << endl 
       << " Chromosome population size   = " << population_size_ << endl 
       << " Extra functions requested    = " << functions << endl 
       << " Maximum Epochs/Generations   = " << Max_Epoch_ << endl;

  FunctionSet = "+-*/AQELGO";

  this->setClasses();
  this->createFunctionSet(functions);
  this->createPopulation();
}

void SprGEP::setRates(double mutation_rate,
		      double RIS_rate,
		      double IS_rate,
		      double One_Point_rate,
		      double Two_Point_rate,
		      double Whole_Gene_rate,
		      double Constant_mutation,
		      double Constant_swap,
		      double Constant_range,
		      double fractionOfHeadsForFuncts) 
{
  mutation_rate_ = mutation_rate;
  RIS_rate_ = RIS_rate;
  IS_rate_ = IS_rate;
  One_Point_rate_ = One_Point_rate;
  Two_Point_rate_ = Two_Point_rate;
  Whole_Gene_rate_ = Whole_Gene_rate;
  Constant_mutation_ = Constant_mutation;
  Constant_swap_ = Constant_swap;
  Constant_range_ = Constant_range;
  fractionOfHeadsForFuncts_ = fractionOfHeadsForFuncts;

  cout <<
    " Mutation Rate                = " << mutation_rate_ << endl <<
    " Root Insertion Rate          = " << RIS_rate_ << endl <<
    " Insertion Rate               = " << IS_rate_ << endl <<
    " One Point Transposition Rate = " << One_Point_rate_ << endl <<
    " Two Point Transposition Rate = " << Two_Point_rate_ << endl <<
    " Whole Gene Transpose Rate    = " << Whole_Gene_rate_ << endl <<
    " Constant Mutation Rate       = " << Constant_mutation_ << endl <<
    " Constant Swap Rate           = " << Constant_swap_ << endl <<
    " Constant Generation Range    = " << Constant_range_ << endl <<
    " Fraction of Heads for Functs = " << fractionOfHeadsForFuncts_ << endl;
}


void SprGEP::createFunctionSet(const std::string& functions) 
{
  // The functions +-*/ are always included
  FunctionsInUse.clear();
  FunctionsInUse.push_back(Plus);
  FunctionsInUse.push_back(Minus);
  FunctionsInUse.push_back(Times);
  FunctionsInUse.push_back(Divide);
  
  for(int i=0;i<functions.size();i++) {
    if( FunctionSet.find(functions[i]) != string::npos ) {
      int f = functionFromChar(functions[i]);
      if( find(FunctionsInUse.begin(),FunctionsInUse.end(),f) 
	  == FunctionsInUse.end()) {
	FunctionsInUse.push_back(f);
      }
    } 
    else {
      cerr << "Unknown function \"" << functions[i] 
	   << "\" requested." << endl;
    }
  }
  cout << "Using functions " << FunctionsInUse << endl;
  function_arg_max_ = 2;
}


bool SprGEP::setValidation(const SprAbsFilter* valData, unsigned valPrint) 
{
  valData_ = valData;
  valPrint_ = valPrint;
  return true;
}


int SprGEP::functionFromChar(char f) 
{
  switch(f) 
    {
    case '+':
      return SprGEP::Plus;
    case '-':
      return SprGEP::Minus;
    case '*':
      return SprGEP::Times;
    case '/':
      return SprGEP::Divide;
    case 'Q':
      return SprGEP::Sqrt;
    case 'E':
      return SprGEP::Exp;
    case 'L':
      return SprGEP::LessThan;
    case 'G':
      return SprGEP::GreaterThan;
    case 'A':
      return SprGEP::Abs;
    case 'O':
      return SprGEP::Log;
    default:
      return -1;
    }
  return -1;
}


char SprGEP::charFromFunction(int f) 
{
  switch(f) 
    {
    case SprGEP::Plus:
      return '+';
    case SprGEP::Minus:
      return '-';
    case SprGEP::Times:
      return '*';
    case SprGEP::Divide:
      return '/';
    case SprGEP::Sqrt:
      return 'Q';
    case SprGEP::Exp:
      return 'E';
    case SprGEP::LessThan:
      return 'L';
    case SprGEP::GreaterThan:
      return 'G';
    case SprGEP::Abs:
      return 'A';	
    case SprGEP::Log:
      return 'O';	
    default:
      return '!';
    }
}


void SprGEP::createPopulation() 
{
  cout << "Creating Chromosome population of size " 
       << population_size_ << endl;
  chromosomes_.clear();
  chromosomes_.resize(population_size_,SprChromosome(&rndm_));

  unsigned size = data_->size();
  unsigned dim = data_->dim();

  unsigned head_length = gene_head_length_;
  unsigned tail_length = gene_head_length_ * (function_arg_max_ - 1) + 1;
  unsigned ch_length = gene_head_length_ + tail_length;

  cout << "Each gene will have a tail of length " << tail_length 
       << " and a total length of " << ch_length << endl;

  // Create the chromosomes, each of 1 or more genes
  int nchromosome = 0;
  for( vector<SprChromosome>::iterator it=chromosomes_.begin();
       it!=chromosomes_.end();it++ ) {
    SprChromosome &thisChromosome = *it;
    thisChromosome.setLength(genes_per_chromosome_,
			     constants_per_chromosome_,
			     genes_per_chromosome_*ch_length);

    for( int ch=0;ch<genes_per_chromosome_;ch++ ) {
      vector<int> head;
      vector<int> tail;

      for( int ih=0;ih<head_length;ih++ ) {
	// head can contain functions or terminals (variables or constants)
	// first position in head must be a function
	if( ih == 0 ) {
	  int ipos = int(floor(rndm_.one()*FunctionsInUse.size()));
	  head.push_back(FunctionsInUse[ipos]);
	} 
	else {
	  // how many head positions will contain functions
	  if( rndm_.one() < fractionOfHeadsForFuncts_ ) {
	    int ipos = int(floor(rndm_.one()*FunctionsInUse.size()));
	    head.push_back(FunctionsInUse[ipos]);
	  } 
	  else {
	    int ipos 
	      = int(floor(rndm_.one()*(dim + constants_per_chromosome_)));
	    if( ipos >= dim ) {
	      head.push_back(Constant);
	    } 
	    else {
	      int var = Variable + ipos;
	      head.push_back(var);
	    }
	  }
	}
      }
      
      for(int j=0;j<tail_length;j++) {
	int idim = int(floor(rndm_.one()*(dim + constants_per_chromosome_)));
	if( idim >= dim ) {
	  tail.push_back(Constant);
	} 
	else {
	  int var = Variable + idim;
	  tail.push_back(var);
	}
      }
      
      thisChromosome.gene_[ch].setHead(head);
      thisChromosome.gene_[ch].setTail(tail);	
      
    }
    // Generate random constants for this chromosome if required
    if( constants_per_chromosome_ > 0 ) {
      thisChromosome.constants_.clear();
      for(int ic=0;ic<constants_per_chromosome_;ic++) {
	double val = Constant_range_*rndm_.one();
	thisChromosome.constants_.push_back(val);
      }
      //cout << "  Constants generated: " << thisChromosome.constants_ << endl;
    }
    thisChromosome.print(cout);
    nchromosome++;
  }
}


void SprGEP::setClasses() 
{
  vector<SprClass> classes;
  data_->classes(classes);
  int size = classes.size();
  if( size > 0 ) cls0_ = classes[0];
  if( size > 1 ) cls1_ = classes[1];
  cout << "Classes for SprGEP are set to " 
       << cls0_ << " " << cls1_ << endl;
}

SprTrainedGEP* SprGEP::makeTrained()
{
  return new SprTrainedGEP(best_);
}


bool SprGEP::reset()
{
  return true;
}


bool SprGEP::setData(SprAbsFilter* data)
{
  assert( data != 0 );
  data_ = data;
  return this->reset();
}


bool SprGEP::printValidation(unsigned cycle) const
{
  // no print-out for zero training cycle
 
  if( cycle == 0 ) return true;

  double sigsig = 0;
  double sigbac = 0;
  double bacsig = 0;
  double bacbac = 0;

  double wtot = 0;
  double loss = 0;

  // loop through validation data
  
  int vsize = valData_->size();
  
  for( int i=0;i<vsize;i++ ) {
    const SprPoint* p = (*valData_)[i];
    double w = valData_->w(i);
    int cls = 0;
    if(      p->class_ == cls0_ ) 
      cls = -1;
    else if( p->class_ == cls1_ ) 
      cls = +1;
    else
      continue;

    double value = best_.Evaluate(p->x_,0);
    wtot += w;

    if( crit_ == 0 ) {
      loss += w*exp(-cls*value);
    }
    else {
      if( cls > 0 ) {
	if( value > 0. )
	  sigsig += w;
	else 
	  sigbac += w;
      } 
      else {
	if( value <= 0. )
	  bacbac += w;
	else 
	  bacsig += w;
      }
    }
  }// end loop thru data

  double fom = SprUtils::min();
  if( crit_ == 0 )
    fom = -loss/wtot;
  else
    fom = crit_->fom(bacbac,bacsig,sigsig,sigbac);
	
  cout << "Validation FOM=" << fom << " at Epoch " << cycle
       << " (SS=" << sigsig << ",BB=" << bacbac 
       << ",SB=" << sigbac << ",BS=" << bacsig <<")" << endl;

  return true;
}


bool SprGEP::train(int verbose) {


  unsigned size = data_->size();
  unsigned dim = data_->dim();
  
  unsigned Epoch = 0; // counts the chromosome generations
  
  double best_fitness = SprUtils::min();
  int lastbest = 0;

  // if no criterion is specified, use correctly classified fraction
  const SprAbsTwoClassCriterion* crit = crit_;
  bool ownCrit = false;
  bool computeExpLoss = false;
  if( crit == 0 ) {
    crit = new SprTwoClassIDFraction();
    ownCrit = true;
    computeExpLoss = true;
  }
  

  // main loop
  while( Epoch < Max_Epoch_ ) {
    
    Epoch++;

    // For each chromosome in the population we calculate its fitness
    // by calculating its specificity and sensitivity over the whole
    // data set
    
    best_fitness = SprUtils::min();
    int best_chromosome = 0;
    
    double sumfitness = 0;
    
    for(int ig=0;ig<population_size_;ig++) {

      pair<double,double> f_and_s
	= this->fitnessWithShift(crit,chromosomes_[ig],computeExpLoss,verbose);

      double fitness = f_and_s.first;
      double shift = f_and_s.second;

      chromosomes_[ig].fitness_ = fitness;
      chromosomes_[ig].shift_  += shift;

      if( fitness > best_fitness ) {
	best_fitness = fitness;
	best_chromosome = ig;
      }

      sumfitness += fitness;
    }// end ig=0;ig<population_size_;ig++

    double avgfitness = sumfitness / population_size_;

    best_ = chromosomes_[best_chromosome];

    if( lastbest != best_chromosome ) {
      lastbest = best_chromosome;	
      if( verbose > 0 ) {
	cout << endl << "Epoch " << Epoch 
	     << " Best chromosome=" << best_chromosome 
	     << " FOM=" << best_fitness << endl;
	best_.print(cout);
      }
    } 
    else if( verbose > 1 ) {
      cout << "Epoch " << Epoch << " Best/Average fitness " << best_fitness 
	   << "/" << avgfitness << endl;
    } 
    else {
      //if(Epoch % 10) cout << "." ;
    }

    if( valPrint_ != 0 && Epoch % valPrint_ == 0 ) {
      if( !this->printValidation(Epoch) ) {
	cerr << "Error printing validation data" << endl;
	if( ownCrit ) delete crit;
	return false;
      }
    }

    // Chromosome Mutation
    // Chromosome Transposition
    // a) by Root Insertion Sequence
    // b) by Insertion Sequence (not at the Root)
    // c) by Gene Transposition (not treated yet, and irrelevant for "+" linking function)
    // Chromosome Recombination
    // a) One Point
    // b) Two Point
    // c) Whole Gene
    // Constant mutation
    // Constant swap
		
    for( int ig=0;ig<population_size_;ig++ ) {
       // the best chromosome is always left alone
      if( ig == best_chromosome ) continue;

      SprChromosome &chromosome = chromosomes_[ig];

      // Mutation
      if( rndm_.one() < mutation_rate_ ) {
	if( verbose > 2 ) { 
	  cout << "Before mutation " << endl; 
	  chromosome.print(cout);
	}
	bool mut = chromosome.mutate(FunctionsInUse, 
				     dim, constants_per_chromosome_);
	if(verbose > 2) { 
	  cout << "After mutation " << endl; 
	  chromosome.print(cout);}
      }

      // Root insertion
      if( rndm_.one() < RIS_rate_ ) {
	if( verbose > 2 ) { 
	  cout << "Before RIS " << endl; 
	  chromosome.print(cout);
	}
	bool ris = chromosome.RIS();
	if( verbose > 2 ) { 
	  cout << "After RIS " << endl; 
	  chromosome.print(cout);
	}
      }
			
      // Insertion
      if( rndm_.one() < IS_rate_ ) {
	if(verbose > 2) { 
	  cout << "Before IS " << endl; 
	  chromosome.print(cout);
	}
	bool is = chromosome.IS();
	if( verbose > 2 ) { 
	  cout << "After IS " << endl; 
	  chromosome.print(cout);
	}
      }

      // One point recombination
      if( rndm_.one() < One_Point_rate_) {
	// Pick another chromosome from the population, randomly
	// but not the best chromosome (which must survive unchanged)
	int iother = best_chromosome;
	while( iother == best_chromosome || iother == ig ) 
	  iother = int(floor(rndm_.one()*population_size_));
	SprChromosome &otherChromosome = chromosomes_[iother];
	if( verbose > 2 ) { 
	  cout << "Before One Point " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);
	}
	bool io = chromosome.OnePoint(otherChromosome);
	if( verbose > 2 ) { 
	  cout << "After One Point " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);
	}
      }

      // Two point recombination
      if( rndm_.one() < Two_Point_rate_ ) {
	// Pick another chromosome from the population, randomly
	// but not the best chromosome (which must survive unchanged)
	int iother = best_chromosome;
	while( iother == best_chromosome || iother == ig ) 
	  iother = int(floor(rndm_.one()*population_size_)); 
	SprChromosome &otherChromosome = chromosomes_[iother];
	if( verbose > 2 ) { 
	  cout << "Before Two Point " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);
	}
	bool it = chromosome.TwoPoint(otherChromosome);
	if( verbose > 2 ) { 
	  cout << "After Two Point " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);
	}
      }

      // Whole gene recombination
      if( rndm_.one() < Whole_Gene_rate_ ) {
	// Pick another chromosome from the population, randomly
	// but not the best chromosome (which must survive unchanged)
	int iother = best_chromosome;
	while( iother == best_chromosome || iother == ig ) 
	  iother = int(floor(rndm_.one()*population_size_)); 
	SprChromosome &otherChromosome = chromosomes_[iother];
	if( verbose > 2 ) { 
	  cout << "Before Whole Gene " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);
	}
	bool wg = chromosome.WholeGene(otherChromosome);
	if(verbose > 2) { 
	  cout << "After Whole Gene " << endl; 
	  chromosome.print(cout); 
	  otherChromosome.print(cout);}
      }

      // Constant mutation
      if( constants_per_chromosome_ > 0 ) {
	if( rndm_.one() < Constant_mutation_ ) {
	  if( verbose > 2 ) { 
	    cout << "Before Constant Mutation " << endl; 
	    chromosome.print(cout);
	  }
	  bool cm = chromosome.ConstantMutation(Constant_range_);
	  if( verbose > 2 ) { 
	    cout << "After Constant Mutation " << endl; 
	    chromosome.print(cout);}
	}
	if( rndm_.one() < Constant_swap_ ) {
	  // Pick another chromosome from the population, randomly
	  // but not the best chromosome (which must survive unchanged)
	  int iother = best_chromosome;
	  while( iother == best_chromosome || iother == ig ) 
	    iother = int(floor(rndm_.one()*population_size_)); 
	  SprChromosome &otherChromosome = chromosomes_[iother];
	  if( verbose > 2 ) { 
	    cout << "Before Constant Swap " << endl; 
	    chromosome.print(cout); 
	    otherChromosome.print(cout);
	  }
	  bool cs = chromosome.ConstantSwap(otherChromosome);
	  if( verbose > 2 ) { 
	    cout << "After Constant Swap " << endl; 
	    chromosome.print(cout); 
	    otherChromosome.print(cout);
	  }
	}
      }
    }
  }// while(Epoch < Max_Epoch_)

  cout << endl << "Finished training. Epoch " << Epoch 
       << " Best chromosome FOM " << best_fitness << endl;
  best_.print(cout);
  if( valPrint_ != 0 ) {
    if( !this->printValidation(Epoch) ) {
      cerr << "Error printing validation data" << endl;
      if( ownCrit ) delete crit;
      return false;
    }
  }

  cout << " C++ code that implements the trained chromosome: " << endl;
  best_.generateCode();

  if( ownCrit ) delete crit;

  return true;
}



char SprGEP::charFromInt(int i) 
{
  char c = '!';
  if(      i < Constant ) {
    c = charFromFunction(i);
  } 
  else if( i == Constant ) {
    c = '?';
  } 
  else if( i >= Variable ) {
    int ic = (i-Variable) % 26; // we map them all to a-z
    c = 'a' + ic;
  }
  return c;
}


void SprGEP::print(std::ostream& os) const
{
  os << "Trained GEP " << SprVersion << endl;
  best_.print(os);
}


std::pair<double,double> SprGEP::fitnessWithShift(
				 const SprAbsTwoClassCriterion* crit,
				 const SprChromosome& chromosome,
				 bool computeExpLoss,
				 int verbose) const
{
  assert( crit != 0 );

  // init cut
  double fitness = SprUtils::min();
  double zcut = 0;

  // store classifier output
  vector<pair<double,double> > r0, r1;
  vector<double> r(data_->size());
  for( int i=0;i<data_->size();i++ ) {
    const SprPoint* p = (*data_)[i];
    double rsp = chromosome.Evaluate(p->x_,0);
    r[i] = rsp;
    if(      p->class_ == cls0_ )
      r0.push_back(pair<double,double>(rsp,data_->w(i)));
    else if( p->class_ == cls1_ ) 
      r1.push_back(pair<double,double>(rsp,data_->w(i)));
  }
    
  // init weights
  double wcor0(0), wmis1(0);
  double wmis0 = 0;
  for( int i=0;i<r0.size();i++ )
    wmis0 += r0[i].second;
  double wcor1 = 0;
  for( int i=0;i<r1.size();i++ )
    wcor1 += r1[i].second;
  assert( wmis0>0 && wcor1>0 );
  double wtot = wmis0+wcor1;
  if( verbose > 2 ) 
    cout << "Optimizing cut for W0=" << wmis0 << " W1=" << wcor1 << endl;
    
  // sort, init min and max, and init iterators
  stable_sort(r.begin(),r.end(),less<double>());
  stable_sort(r0.begin(),r0.end(),SGEPCmpPairDDFirst());
  stable_sort(r1.begin(),r1.end(),SGEPCmpPairDDFirst());
  vector<pair<double,double> >::iterator i0start(r0.begin()), 
    i0split(r0.begin());
  vector<pair<double,double> >::iterator i1start(r1.begin()),
    i1split(r1.begin());
  if( verbose > 2 ) {
    if( !r0.empty() ) {
      cout << "Classifier range for 0: " 
	   << r0[0].first << " " << r0[r0.size()-1].first << endl;
    }
    if( !r1.empty() ) {
      cout << "Classifier range for 1: " 
	   << r1[0].first << " " << r1[r1.size()-1].first << endl;
    }
  }
    
  // fill out divisions
  vector<double> division;
  division.push_back(SprUtils::min());
  double xprev = r[0];
  for( int k=1;k<r.size();k++ ) {
    double xcurr = r[k];
    if( (xcurr-xprev) > SprUtils::eps() ) {
      division.push_back(0.5*(xcurr+xprev));
      xprev = xcurr;
    }
  }
  division.push_back(SprUtils::max());

  // find optimal point
  int ndiv = division.size();
  vector<double> f(ndiv);
  for( int k=0;k<ndiv;k++ ) {
    double z = division[k];
    i0split = find_if(i0start,r0.end(),
		      not1(bind2nd(SGEPCmpPairDDFirstNumber(),z)));
    i1split = find_if(i1start,r1.end(),
		      not1(bind2nd(SGEPCmpPairDDFirstNumber(),z)));
    for( vector<pair<double,double> >::iterator iter=i0start;
	 iter!=i0split;iter++ ) {
      wcor0 += iter->second;
      wmis0 -= iter->second;
    }
    for( vector<pair<double,double> >::iterator iter=i1start;
	 iter!=i1split;iter++ ) {
      wmis1 += iter->second;
      wcor1 -= iter->second;
    }
    i0start = i0split;
    i1start = i1split;
    f[k] = crit->fom(wcor0,wmis0,wcor1,wmis1);
  }

  // find point giving min misclassification fraction
  vector<double>::const_iterator imax = max_element(f.begin(),f.end());
  int k = imax - f.begin();
  fitness = *imax;
  zcut = division[k];

  // update exponential loss
  if( computeExpLoss ) {
    double wtot = 0;
    double aveloss = 0;
    for( int i=0;i<r0.size();i++ ) {
      aveloss += exp(r0[i].first);
      wtot += r0[i].second;
    }
    for( int i=0;i<r1.size();i++ ) {
      aveloss += exp(-r1[i].first);
      wtot += r1[i].second;
    }
    assert( wtot > 0 );
    aveloss /= wtot;
    fitness = -aveloss;
  }

  // exit
  return pair<double,double>(fitness,zcut);
}
