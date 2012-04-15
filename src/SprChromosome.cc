//$Id: SprChromosome.cc,v 1.2 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprGEP.hh"
#include "StatPatternRecognition/SprChromosome.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprRandomNumber.hh"

#include <stdio.h>
#include <cassert>
#include <utility>
#include <algorithm>
#include <functional>
#include <cmath>
#include <sstream>

using namespace std;


SprChromosome::~SprChromosome() 
{
  if( ownRndm_ ) {
    delete rndm_;
    rndm_ = 0;
  }
}


SprChromosome::SprChromosome() 
  :
  num_genes_(1),
  link_function_('+'),
  gene_(1),
  constants_(),
  fitness_(0),
  length_(0),
  shift_(0),
  rndm_(new SprRandomNumber()),
  ownRndm_(true)
{}

SprChromosome::SprChromosome(SprRandomNumber* rndm) 
  :
  num_genes_(1),
  link_function_('+'),
  gene_(1),
  constants_(),
  fitness_(0),
  length_(0),
  rndm_(rndm),
  ownRndm_(false)
{}


SprChromosome::SprChromosome(SprRandomNumber* rndm,
			     unsigned num_genes,
			     char link_function)
  :
  num_genes_(num_genes),
  link_function_(link_function),
  gene_(1),
  constants_(),
  fitness_(0),
  length_(0),
  shift_(0),
  rndm_(rndm),
  ownRndm_(false)
{
  assert( num_genes >= 1 );
  gene_.resize(num_genes);
}


void SprChromosome::generateCode() 
{
  ostringstream outs;
  int iconstant = 0;

  cout << "bool Signal(double x[]) {" << endl;

  if( link_function_ == '/' || link_function_ == '*' ) {
    cout << "     double Chromosome_value = 1;" << endl;
  } 
  else {
    cout << "     double Chromosome_value = 0;" << endl;
  }

  for( int ch=0; ch<num_genes_; ch++ ) {

    SprGene &gene = gene_[ch];
    
    vector<vector<int> > tree;
    gene.generateTree(tree);
    
    int level = tree.size()-1;
    
    vector<string> vrow, vrowbelow;

    string value;
    double val;

    for( int i=level;i>=0;i-- ) {
      int iposbelow = 0;
      for( int ic=tree[i].size()-1;ic>=0;ic-- ) {
	int c = tree[i][ic];
	if(      c == SprGEP::Constant ) {
	  // pick the next constant from the chromosome
	  val = constants_.at(iconstant++);
	  if( iconstant >= constants_.size() ) iconstant = 0;
	  outs << val;
	  vrow.push_back(outs.str());
	  outs.str("");
	} 
	else if( c >= SprGEP::Variable ) {
	  int index = c - SprGEP::Variable;
	  outs << index;
	  value = "x[" + outs.str() + "]";
	  outs.str("");
	  vrow.push_back(value);
	} 
	else {
	  switch(c) 
	    {
	    case SprGEP::Sqrt:
	      value = "sqrt(" + vrowbelow.at(iposbelow) + ")";
	      iposbelow++;
	      vrow.push_back(value);
	      break;
	    case SprGEP::Exp:
	      value = "exp(" + vrowbelow.at(iposbelow) + ")";
	      iposbelow++;
	      vrow.push_back(value);
	      break;
	    case SprGEP::Abs:
	      value = "abs(" + vrowbelow.at(iposbelow) + ")";
	      iposbelow++;
	      vrow.push_back(value);
	      break;
	    case SprGEP::Log:
	      value = "log(" + vrowbelow.at(iposbelow) + ")";
	      iposbelow++;
	      vrow.push_back(value);
	      break;
	    case SprGEP::Times:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ "*" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    case SprGEP::Divide:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ "/" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    case SprGEP::Plus:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ "+" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    case SprGEP::Minus:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ "-" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    case SprGEP::LessThan:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ "<" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    case SprGEP::GreaterThan:
	      value = "(" + vrowbelow.at(iposbelow+1) 
		+ ">" + vrowbelow.at(iposbelow) + ")";
	      iposbelow += 2;
	      vrow.push_back(value);
	      break;
	    default:
	      cerr << "Error parsing value tree" << endl;
	    }
	}
      }
      vrowbelow = vrow;
      vrow.clear();
    }
    assert (vrowbelow.size() == 1);
    
    cout << "     Chromosome_value " << outs.str() << link_function_ 
	 << "= " << vrowbelow[0] << ";" << endl;

  }

  cout << "     Chromosome_value += " << -shift_ << ";" << endl;

  cout << "     return Chromosome_value > 0;" << endl;
  cout << "}" << endl;
}
	

double SprChromosome::Evaluate(const std::vector<double>& data, 
			       int verbose) const
{
  // Parse the gene(s) and evaluate each with the given data
  vector<double> gene_value;
  
  int iconstant = 0;
  
  for(int ch=0; ch<num_genes_; ch++) {
    
    const SprGene &gene = gene_[ch];
    
    vector<vector<int> > tree;
    gene.generateTree(tree);
    
    int level = tree.size();
    
    if( verbose > 1 ) {
      cout << "Evaluating Gene #" <<  ch << ":";
      gene.print(cout);
      cout << endl;
      //if(constants_.size() > 0) cout << "  Constants " << constants_ << endl;
    }
    
    if( verbose > 1 ) {
      cout << "Finished parsing Gene" << endl;
      for( int i=0;i<level;i++ ) {
	cout << "Level " << i << ": ";
	for( int it=0;it<tree[i].size();it++ ) 
	  cout << SprGEP::charFromInt(tree[i][it]);
	cout << endl;
      }
    }

    level--;

    vector<double> vrow;
    vector<double> vbelow;
    double val;
    for( int i=level;i>=0;i-- ) {
      int iposbelow = 0;
      for( int ic=tree[i].size()-1;ic>=0;ic-- ) {
	int c = tree[i][ic];
	if(     c == SprGEP::Constant ) {
	  // pick the next constant from the chromosome
	  val = constants_.at(iconstant++);
	  if( iconstant >= constants_.size() ) iconstant = 0;
	  vrow.push_back(val);	
	} 
	else if( c >= SprGEP::Variable ) {
	  int index = c - SprGEP::Variable;
	  val = data.at(index);
	  vrow.push_back(val);
	} 
	else {
	  switch(c) {
	  case SprGEP::Sqrt:
	    // need to protect against sqrt neg number
	    if( vbelow.at(iposbelow) < 0 ) {
	      val = 0;
	    } 
	    else {
	      val = sqrt(vbelow.at(iposbelow));
	    }
	    iposbelow++;
	    vrow.push_back(val);
	    break;
	  case SprGEP::Abs:
	    // absolute value
	    val = abs(vbelow.at(iposbelow));
	    iposbelow++;
	    vrow.push_back(val);
	    break;
	  case SprGEP::Log:
	    // log base e
	    val = log(vbelow.at(iposbelow));
	    iposbelow++;
	    vrow.push_back(val);
	    break;
	  case SprGEP::Exp:
	    // exponential
	    val = exp(vbelow.at(iposbelow));
	    iposbelow++;
	    vrow.push_back(val);
	    break;
	  case SprGEP::LessThan:
	    // Less than
	    val = 0;
	    if(vbelow.at(iposbelow+1) < vbelow.at(iposbelow)) val = 1;
	    iposbelow += 2;
	    vrow.push_back(val);
	    break;
	  case SprGEP::GreaterThan:
	    // Greater than
	    val = 0;
	    if(vbelow.at(iposbelow+1) > vbelow.at(iposbelow)) val = 1;
	    iposbelow += 2;
	    vrow.push_back(val);
	    break;
	  case SprGEP::Times:
	    val = vbelow.at(iposbelow+1) * vbelow.at(iposbelow);
	    iposbelow += 2;
	    vrow.push_back(val);
	    break;
	  case SprGEP::Divide:
	    // need to protect against divide by zero
	    if( vbelow.at(iposbelow) == 0 ) {
	      val = SprUtils::max();
	    } 
	    else {
	      val = vbelow.at(iposbelow+1) / vbelow.at(iposbelow);
	    }
	    iposbelow += 2;
	    vrow.push_back(val);
	    break;
	  case SprGEP::Plus:
	    val = vbelow.at(iposbelow+1) + vbelow.at(iposbelow);
	    iposbelow += 2;
	    vrow.push_back(val);
	    break;
	  case SprGEP::Minus:
	    val = vbelow.at(iposbelow+1) - vbelow.at(iposbelow);
	    iposbelow += 2;
	    vrow.push_back(val);
	    break;
	  default:
	    cerr << "Error parsing value tree" << endl;
	    return 0;
	  }
	}
      }
      //if(verbose) cout << "Expression tree row " 
      // << i << " (reversed) evaluated to " << vrow << endl;
      vbelow = vrow;
      vrow.clear();
    }
    assert (vbelow.size() == 1);
    gene_value.push_back(vbelow.at(0));
    if(verbose > 1) cout << "Evaluates to: " << vbelow.at(0) << endl << endl;
  }
  
  // Calculate the return value of the chromosome. 
  // This depends on the linking function and
  // the number of genes. If there is only one gene, then the linking function
  // is not used, and the chromosome has the same value as the gene.
  
  double result = gene_value.at(0);

  if(num_genes_ == 1) return result;
	
  for(int ch=1;ch<num_genes_;ch++) {
    if(      link_function_ == '*' ) {
      result *= gene_value.at(ch);
    } 
    else if( link_function_ == '/' ) {
      result /= gene_value.at(ch);
    } 
    else if( link_function_ == '+' ) {
      result += gene_value.at(ch);
    } 
    else if( link_function_ == '-' ) {
      result -= gene_value.at(ch);
    }
  }

  // shift
  result -= shift_;

  if(verbose > 1) cout << "Chromosome evaluates to: " << result << endl;
  
  return result;
}


void SprGene::print(std::ostream& os) const 
{
  vector<int> sequence(head_);
  int lhead = sequence.size();
  sequence.insert(sequence.end(),tail_.begin(),tail_.end());

  for(int is=0;is<sequence.size();is++) {
    if(is == lhead) os << "|";
    int val = sequence.at(is);
    os << SprGEP::charFromInt(val);
  }
}


void SprChromosome::print(std::ostream& os) const 
{
  for(int ch=0;ch<num_genes_;ch++) {
    gene_.at(ch).print(os);
    os << " ";
  }
  os << " L= " << link_function_;
  if( !constants_.empty() ) os << " C=";
  for( int i=0;i<constants_.size();i++ )
    os << " " << constants_[i];
  os << " S= " << shift_;
  os << endl;
}


bool SprChromosome::ConstantMutation(double range) 
{
  // pick a constant at random
  int ic = int(floor(rndm_->one()*constants_.size()));
  // replace it with a new one
  constants_[ic] = rndm_->one() * range;
  return true;
}

bool SprChromosome::ConstantSwap(SprChromosome &other) 
{
  // choose random constants from each chromosome and swap them
  int ic1 = int(floor(rndm_->one()*constants_.size()));
  int ic2 = int(floor(rndm_->one()*constants_.size()));
  swap(constants_[ic1],other.constants_[ic2]);
  /*
    double val = constants_[ic1];
    constants_[ic1] = other.constants_[ic2];
    other.constants_[ic2] = val;
  */
  return true;
}

bool SprChromosome::WholeGene(SprChromosome &other) 
{
  // Choose a random gene from both chromosomes, and swap them
  int ig1 = int(floor(rndm_->one()*num_genes_));
  int ig2 = int(floor(rndm_->one()*num_genes_));
  swap(gene_[ig1],other.gene_[ig2]);
  /*
    vector<int> head1 = gene_[ig1].getHead();
    vector<int> tail1 = gene_[ig1].getTail();
    gene_[ig1].setHead(other.gene_[ig2].getHead());
    gene_[ig1].setTail(other.gene_[ig2].getTail());
    other.gene_[ig2].setHead(head1);
    other.gene_[ig2].setTail(tail1);
  */
  return true;
}


bool SprChromosome::OnePoint(SprChromosome &other) 
{
  // We choose a random point in the chromosome, and exchange material
  // The assumptions is that chromosomes are the same size as are their genes
  vector<int> chromosome1 = gene_[0].getHead();
  int lhead = chromosome1.size();
  vector<int> tail1 = gene_[0].getTail();
  int ltail = tail1.size();
  chromosome1.insert(chromosome1.end(),tail1.begin(),tail1.end());
  vector<int> chromosome2 = other.gene_[0].getHead();
  vector<int> tail2 = other.gene_[0].getTail();
  chromosome2.insert(chromosome2.end(),tail2.begin(),tail2.end());
  
  for(int i=1;i<num_genes_;i++) {
    vector<int> head1 = gene_[i].getHead();
    tail1 = gene_[i].getTail();
    vector<int> head2 = other.gene_[i].getHead();
    tail2 = other.gene_[i].getTail();
    chromosome1.insert(chromosome1.end(),head1.begin(),head1.end());
    chromosome1.insert(chromosome1.end(),tail1.begin(),tail1.end());
    chromosome2.insert(chromosome2.end(),head2.begin(),head2.end());
    chromosome2.insert(chromosome2.end(),tail2.begin(),tail2.end());
  }
  assert(chromosome1.size() == chromosome2.size());
  int len = chromosome1.size();
  int pos = int(floor(rndm_->one()*len));
  
  // Swap the contents from position pos to the end of the chromosomes
  swap_ranges(chromosome1.begin()+pos, 
	      chromosome1.end(),
	      chromosome2.begin()+pos);
  
  pos = 0;
  int ig = 0;
  while( pos < len && ig < num_genes_ ) {
    vector<int> new1(chromosome1.begin()+pos,chromosome1.begin()+pos+lhead);
    vector<int> new2(chromosome2.begin()+pos,chromosome2.begin()+pos+lhead);
    gene_[ig].setHead(new1);
    other.gene_[ig].setHead(new2);
    pos += lhead;
    vector<int> new3(chromosome1.begin()+pos,chromosome1.begin()+pos+ltail);
    vector<int> new4(chromosome2.begin()+pos,chromosome2.begin()+pos+ltail);
    gene_[ig].setTail(new3);
    other.gene_[ig].setTail(new4);
    pos += ltail;
    ig++;
  }
  
  return true;
}

bool SprChromosome::TwoPoint(SprChromosome &other) 
{
  // We choose two random points in the chromosome, and exchange material
  // The assumptions is that chromosomes are the same size as are their genes
  vector<int> chromosome1 = gene_[0].getHead();
  int lhead = chromosome1.size();
  vector<int> tail1 = gene_[0].getTail();
  int ltail = tail1.size();
  chromosome1.insert(chromosome1.end(),tail1.begin(),tail1.end());
  vector<int> chromosome2 = other.gene_[0].getHead();
  vector<int> tail2 = other.gene_[0].getTail();
  chromosome2.insert(chromosome2.end(),tail2.begin(),tail2.end());
  
  for(int i=1;i<num_genes_;i++) {
    vector<int> head1 = gene_[i].getHead();
    tail1 = gene_[i].getTail();
    vector<int> head2 = other.gene_[i].getHead();
    tail2 = other.gene_[i].getTail();
    chromosome1.insert(chromosome1.end(),head1.begin(),head1.end());
    chromosome1.insert(chromosome1.end(),tail1.begin(),tail1.end());
    chromosome2.insert(chromosome2.end(),head2.begin(),head2.end());
    chromosome2.insert(chromosome2.end(),tail2.begin(),tail2.end());
  }
  assert(chromosome1.size() == chromosome2.size());
  int len = chromosome1.size();
  
  int pos1 = int(floor(rndm_->one()*(len-1)));
  int pos2 = pos1 + int(floor(rndm_->one()*(len-pos1)));
  
  assert(pos1 <= pos2);
  int nchars = pos2-pos1+1;
  
  // Swap the contents from position pos1 to position pos2
  swap_ranges(chromosome1.begin()+pos1, chromosome1.begin()+pos2,
	      chromosome2.begin()+pos1);
  
  int pos = 0;
  int ig = 0;
  while(pos < len && ig < num_genes_) {
    vector<int> new1(chromosome1.begin()+pos,chromosome1.begin()+pos+lhead);
    vector<int> new2(chromosome2.begin()+pos,chromosome2.begin()+pos+lhead);
    gene_[ig].setHead(new1);
    other.gene_[ig].setHead(new2);
    pos += lhead;
    vector<int> new3(chromosome1.begin()+pos,chromosome1.begin()+pos+ltail);
    vector<int> new4(chromosome2.begin()+pos,chromosome2.begin()+pos+ltail);
    gene_[ig].setTail(new3);
    other.gene_[ig].setTail(new4);
    pos += ltail;
    ig++;
  }
  
  return true;
}

bool SprChromosome::mutate(const std::vector<int>& funcs, 
			   unsigned dim, unsigned nconsts) 
{

  // select a random gene
  int ic = int(floor(rndm_->one()*num_genes_));
  vector<int> head = gene_[ic].getHead(); 
  vector<int> tail = gene_[ic].getTail();
  int lhead = head.size();
  int ltail = tail.size();
  // select a random position in the gene
  int posinc = int(floor(rndm_->one()*(lhead+ltail)));
  // in the head we can change to a function or a terminal
  // in all cases, we force a real mutation, 
  // i.e. the terminal/function must change
  if( posinc < lhead ) {
    // note that the first position in the head must be a function
    int headval = head.at(posinc);
    if( posinc == 0 || rndm_->one() < 0.5 ) {
      // mutate to a (different) function
      int ifunc = funcs.at(int(floor(rndm_->one()*funcs.size())));
      while(ifunc == headval) 
	ifunc = funcs.at(int(floor(rndm_->one()*funcs.size())));
      head[posinc] = ifunc;
    } 
    else {
      int i = int(floor(rndm_->one()*(dim+nconsts)));
      int cterm = SprGEP::Constant;
      if(i < dim) {
	cterm = SprGEP::Variable + i;	
	while(cterm == headval) 
	  cterm = SprGEP::Variable + int(floor(rndm_->one()*dim));
      }
      head[posinc] = cterm;
    }
  } 
  else {
    // in the tail we can only mutate to a terminal or constant
    int tailval = tail[posinc-lhead];
    int i = int(floor(rndm_->one()*(dim+nconsts)));
    int cterm = SprGEP::Constant;
    if(i < dim) {
      cterm = SprGEP::Variable + int(floor(rndm_->one()*dim));
      while(cterm == tailval) 
	cterm = SprGEP::Variable + int(floor(rndm_->one()*dim));
    }
    tail[posinc-lhead] = cterm;
  }
  gene_[ic].setHead(head);
  gene_[ic].setTail(tail);
  return true;
}


bool SprChromosome::RIS() 
{
  // select a random gene
  int ic = int(floor(rndm_->one()*num_genes_));
  vector<int> head = gene_[ic].getHead(); 
  int lhead = head.size();
  // select a random position in the head, but not the start
  int poshead = int(floor(rndm_->one()*lhead));
  while( poshead == 0 ) poshead = int(floor(rndm_->one()*lhead));
  // find first downstream function
  for( int i=poshead;i<head.size();i++ ) {
    if(head[i] < SprGEP::Constant) { // in the function space
      int len = lhead-i;
      vector<int> ris(head.begin()+i,head.end());
      // move this function headed sequence to the head of the gene
      head.insert(head.begin(),ris.begin(),ris.end());
      // discard the parts of the head that were shifted out
      head.erase(head.begin()+lhead,head.end());
      gene_[ic].setHead(head);
      return true;
    }
  }
  return false;
}

bool SprChromosome::IS() 
{
  // select a random gene
  int ic = int(floor(rndm_->one()*num_genes_));
  vector<int> gene(gene_[ic].getHead());
  int lhead = gene.size();
  vector<int> tail(gene_[ic].getTail());
  gene.insert(gene.end(),tail.begin(),tail.end());
  int len = gene.size();
  // select a random position in the gene
  int pos = int(floor(rndm_->one()*len));
  // select a random length for the transposon
  int lseq = 1 + int(floor(rndm_->one()*(len-pos-1)));
  vector<int> transposon(gene.begin()+pos,gene.begin()+pos+lseq);
  // select a random position in the chromosome
  // first select a random gene
  ic = int(floor(rndm_->one()*num_genes_));
  vector<int> head = gene_[ic].getHead();
  pos = int(floor(rndm_->one()*lhead));
  // ensure the position is not at the first position in the head
  while( pos == 0 ) pos = int(floor(rndm_->one()*lhead));
  head.insert(head.begin()+pos,transposon.begin(),transposon.end());
  head.erase(head.begin()+lhead,head.end());
  gene_[ic].setHead(head);
  return true;
}


void SprGene::generateTree(std::vector<std::vector<int> >& tree) const 
{
  vector<int> geneVector = head_;
  geneVector.insert(geneVector.end(),tail_.begin(),tail_.end());

  tree.clear();
  tree.resize(SprGEP::Max_Depth);

  vector<int> nargsreq(SprGEP::Max_Depth,0);
  int ipos = 0;
  int level = 0;
  nargsreq[level] = 1;
  while( ipos < geneVector.size() ) {
    int c = geneVector[ipos];
    tree[level].push_back(c);
    nargsreq[level]--;
    if(c == SprGEP::Constant || c >= SprGEP::Variable) { 
    } 
    else {
      switch(c) 
	{
	case SprGEP::Sqrt:
	case SprGEP::Exp:
	case SprGEP::Abs:
	case SprGEP::Log:
	  nargsreq[level+1] += 1;
	  break;
	case SprGEP::Plus:
	case SprGEP::Minus:
	case SprGEP::Times:
	case SprGEP::Divide:
	case SprGEP::LessThan:
	case SprGEP::GreaterThan:
	  nargsreq[level+1] += 2;
	  break;
	default:
	  cerr << "Unrecognised function in gene \"" << c << "\"" << endl;
	  return;
	}
    }
    if(nargsreq[level] == 0) {
      if(nargsreq[level+1] == 0) break;
      level++;
    }
    ipos++;
    if(level+1 >= SprGEP::Max_Depth) {
      cerr << "Max expression tree depth " << SprGEP::Max_Depth 
	   << " reached for C=";
      this->print(cout);
      cout << endl;
      return;
    }
  }
}
