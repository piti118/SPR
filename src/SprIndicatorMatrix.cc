//$Id: SprIndicatorMatrix.cc,v 1.2 2008-05-09 21:25:26 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprIndicatorMatrix.hh"

#include <stdlib.h>
#include <cassert>
#include <string>
#include <cstring>
#include <cmath>

using namespace std;

const double SprIndicatorMatrix::tolerance_ = 0.5;
const int SprIndicatorMatrix::maxExhaustiveRows_ = 30;
SprRandomNumber SprIndicatorMatrix::rndm_;

SprIndicatorMatrix::SprIndicatorMatrix(int p,int q)
  : 
  nrow_(p), 
  ncol_(q), 
  size_(p*q), 
  m_(new int[size_])
{}

SprIndicatorMatrix::SprIndicatorMatrix(const SprMatrix& matrix)
  : 
  nrow_(matrix.num_row()),
  ncol_(matrix.num_col()),
  size_(matrix.num_row()*matrix.num_col()),
  m_(new int[size_])
{  
  for(int i=0;i<nrow_;i++){
    for(int j=0;j<ncol_;j++){
      double doubleCode = matrix[i][j];
      int code = ( doubleCode>tolerance_ ? 1 : 
		   (( doubleCode<(-1*tolerance_)) ? -1 : 0 ) );
      this->set(i,j,code);
    }
  }
}

SprIndicatorMatrix::SprIndicatorMatrix(const SprIndicatorMatrix& matrix)
  :
  nrow_(matrix.nrow_),
  ncol_(matrix.ncol_),
  size_(matrix.size_),
  m_(0)
{
  m_ = new int[size_];
  memcpy(m_, matrix.m_, size_*sizeof(int));
}

SprIndicatorMatrix* SprIndicatorMatrix::randomDense(
	int numrow, int numcolumn, 
	double psignal, double pbackground, 
	int numtrial,
	MatrixMeasure measure)
{
  return SprIndicatorMatrix::randomSparse(numrow, 
					  numcolumn, 
					  psignal, 
					  pbackground, 
					  0 , 
					  numtrial, 
					  measure);
}

SprIndicatorMatrix* SprIndicatorMatrix::randomSparse(
	int numrow, int numcolumn, 
	double psignal, double pbackground, double pignore, 
	int numtrial,
	MatrixMeasure measure)
{
  assert(psignal>=0);
  assert(pbackground>=0);
  assert(pignore>=0);
  assert(psignal!=0 || pbackground!=0 || pignore!=0);
  
  //normalize the probablity
  double normalizationConstant = psignal+pbackground+pignore;
  double npsignal = psignal/normalizationConstant;
  //double npbackground = pbackground/normalizationConstant;
  double npignore = pignore/normalizationConstant;
  
  //cutoff for (flat) randomly generated number between 0-1 
  //  to be called into each category
  double ignoreCutOff=npignore;
  double signalCutOff=ignoreCutOff+npsignal; 
  //double backgroundCutOff = 1;//ignore rounding error
  
  //measure function
  SprIndicatorMatrix temp(numrow,numcolumn);
  SprIndicatorMatrix* currentBest=0;
  double currentBestMeasure=0;
  bool first=true;
  for(int n=0; n<numtrial; n++){//loop over numtrails
    //generate matrix
    for(int i=0;i<numrow;i++){
      for(int j=0;j<numcolumn;j++){
	double thisNumber = rndm_.one();
	//getting code for given random number using cutoffs
	int thisCode 
	  = ( thisNumber>ignoreCutOff ? ((thisNumber>signalCutOff)?-1:1) : 0 );
	temp.set(i,j,thisCode);	
      }
    }
    //check the code
    
    if( temp.checkMatrix() ){
      //calculate hamming distance and update the curretBest if necessary
      double thisMeasure = temp.evaluate(measure);
      if(thisMeasure>currentBestMeasure || first){
	first = false;
	if(currentBest!=0){delete currentBest; currentBest=0;}
	currentBest = new SprIndicatorMatrix(temp);
	currentBestMeasure = thisMeasure;
      }
    }
  }
  
  return currentBest;
}

SprIndicatorMatrix* SprIndicatorMatrix::randomHillClimbing(
        int row, int column, int numstep, double pKeepBadChange,
	MatrixMeasure measure)
{
  // init
  SprIndicatorMatrix current(row,column);
  bool valid = false;

  //fill in current matrix 
  for( int stepcount=0;stepcount<numstep;stepcount++ ){
    for(int i=0;i<row;i++){
      for(int j=0;j<column;j++){
	current.set(i,j,( rndm_.one()<0.5 ? -1 : 1 )); 
      }
    }
    valid = current.checkMatrix();
    if( valid ) break;
  }
  if( !valid ) return 0;

  // minimize
  unsigned validSteps(0), goodSteps(0);
  double bestMeasure = current.evaluate(measure);
  for( int stepcount=0;stepcount<numstep;stepcount++ ) {
    bool breakEqualCols = ( measure==MINROW || rndm_.one()<0.5 );

    vector<pair<int,int> > pairs;
    if( breakEqualCols )
      pairs = current.closestRowPairs();
    else
      pairs = current.closestColPairs();
    int pairindex = int(floor(rndm_.one()*pairs.size()));
    assert( pairindex < pairs.size() );
    pair<int,int> mypair = pairs[pairindex];
    int i1 = mypair.first;
    int i2 = mypair.second;
    int swap1 = ( rndm_.one()<0.5 ? i1 : i2 );

    vector<int> swaps;
    if( breakEqualCols )
      swaps = current.equalColsForRows(i1,i2);
    else
      swaps = current.equalRowsForCols(i1,i2);
    int swapindex = int(floor(rndm_.one()*swaps.size()));
    assert( swapindex < swaps.size() );
    int swap2 = swaps[swapindex];

    //do the swap
    int swaprow = swap1;
    int swapcol = swap2;
    if( !breakEqualCols ) {
      swaprow = swap2;
      swapcol = swap1;
    }
    current.swap(swaprow,swapcol);

    //check if it's valid
    if(!current.checkMatrix()){
      //not valid undo the change and continue
      //TODO: caching
      current.swap(swaprow,swapcol);		
      continue;
    }

    //if it doesn't improve the measure then undo and continue
    validSteps++;
    double thisMeasure = current.evaluate(measure);
    if(thisMeasure > bestMeasure){
      //keep the change
      bestMeasure = thisMeasure;
      goodSteps++;
    }else{
      //undo the change
      if( rndm_.one()>pKeepBadChange || stepcount==(numstep-1) )
	current.swap(swaprow,swapcol);		
    }
  }// end stepcount loop

  // how many valid steps
  cout << "Random hill climbing converged after "
       << validSteps << " valid steps; of these " 
       << goodSteps << " steps improved the separation measure." << endl;

  // exit  
  return new SprIndicatorMatrix(current);
}

void SprIndicatorMatrix::swap(int swaprow, int swapcol)
{
  this->set(swaprow,swapcol,-1*this->get(swaprow,swapcol));
}

vector<int> SprIndicatorMatrix::equalColsForRows(int row1, int row2) const
{
  int numcol = this->num_col();
  vector<int> result;
  for(int i=0;i<numcol;i++){
    if(this->get(row1,i)==this->get(row2,i)){
      result.push_back(i);
    }
  }
  return result;
}

vector<int> SprIndicatorMatrix::equalRowsForCols(int col1, int col2) const
{
  int numrow = this->num_row();
  vector<int> result;
  for(int i=0;i<numrow;i++){
    if(this->get(i,col1)==this->get(i,col2)){
      result.push_back(i);
    }
  }
  return result;
}

//return vector of two element vector for row#'s of cloest pairs
vector<pair<int,int> > SprIndicatorMatrix::closestRowPairs()
{
  int numrow = this->num_row();
  int numcol = this->num_col();
  int currentMin = numcol;
  
  vector<pair<int,int> > result;
  for(int i=0;i<numrow;i++){
    for(int j=i+1;j<numrow;j++){
      int hammingDistance = this->rowHammingDistance(i,j);
      if(hammingDistance < currentMin){
	currentMin = hammingDistance;
	result.clear();
	result.push_back(pair<int,int>(i,j));
      }
      else if(hammingDistance==currentMin){
	result.push_back(pair<int,int>(i,j));
      } 
    }
  }
  return result;
}

//return vector of two element vector for row#'s of cloest pairs
vector<pair<int,int> > SprIndicatorMatrix::closestColPairs()
{
  int numrow = this->num_row();
  int numcol = this->num_col();
  int currentMin = numcol;
  
  vector<pair<int,int> > result;
  for(int i=0;i<numcol;i++){
    for(int j=i+1;j<numcol;j++){
      int hammingDistance = this->columnHammingDistance(i,j);
      if(hammingDistance < currentMin){
	currentMin = hammingDistance;
	result.clear();
	result.push_back(pair<int,int>(i,j));
      }
      else if(hammingDistance==currentMin){
	result.push_back(pair<int,int>(i,j));
      } 
    }
  }
  return result;
}

//matrix of One Vs One append to One Vs All
//this kind of matrix will pick only easy binary problems
SprIndicatorMatrix* SprIndicatorMatrix::ovoova(int row)
{
  int numcol = row*(row-1)/2+row;
  SprIndicatorMatrix* result = new SprIndicatorMatrix(row,numcol);
  for(int i=0;i<row;i++){
    for(int j=0;j<row;j++){//fill in one vs all
      result->set(i,j,i==j?1:-1);
    }
    for(int j=row;j<numcol;j++){//fill in 0 for the one vs one
      result->set(i,j,0);
    }
  }
  int l=row;
  for(int j=0;j<row;j++){//fill in one vs one
    for(int k=j+1;k<row;k++){
      result->set(j,l,1);
      result->set(k,l,-1);
      l++;
    }
  }  
  return result;
}

//return an exhaustive code matrix
SprIndicatorMatrix* SprIndicatorMatrix::exhaustive(int row)
{
  if( row > maxExhaustiveRows_ ) { 
    cerr << "too many rows. remember number of columns grows like 2^(n-1)-1)"
	 << endl; 
    return 0; 
  }
  int numcol = (1<<(row-1))-1;//2^(n-1)-1
  SprIndicatorMatrix* result = new SprIndicatorMatrix(row,numcol);
  for(int i=0;i<row;i++){
    for(int j=0;j<numcol;j++){
      //the line below means
      // (col / 2^row) mod 2 then (result*2)-1 to convert from 0,1 to -1,1
      int code = (((j+1)/(1<<i))%2<<1)-1;//believe me that's all you need
      result->set(i,j,code);
    }
  }
  return result;
}

bool SprIndicatorMatrix::checkMatrix() const
{
  int numrow = this->num_row();
  int numcol = this->num_col();
  //check for duplicated row 
  for(int firstRow=0; firstRow< numrow; firstRow++){
    for(int secondRow=firstRow+1; secondRow<numrow; secondRow++){
      bool dupRow = true;
      for(int j=0; j<numcol; j++){
	int x=this->get(firstRow,j);
	int y=this->get(secondRow,j);
	if(x!=0 || y!=0){	//ignore the ignore bit
	  dupRow &= (x==y); //duplication
	}
	if(!dupRow){
	  break; //these two rows are fine
	}
      }
      if(dupRow){
	return false;
      }
    }		
  }
  
  //now check for duplicated/complement column	
  //check for duplicated row 
  for(int firstCol=0; firstCol< numcol; firstCol++){
    for(int secondCol=firstCol+1; secondCol<numcol; secondCol++){
      bool dupCol = true;
      bool compCol = true;
      for(int i=0; i<numrow; i++){
	
	int x=this->get(i,firstCol);
	int y=this->get(i,secondCol);
	if(x!=0 || y!=0){	//ignore the ignore bit
	  dupCol &= (x==y); //duplication
	  compCol &= (x+y==0); //complement (x==~y)
	}
	if(!dupCol && !compCol){
	  break; //these two columns are fine
	}
      }
      if(dupCol || compCol){
	return false;
      }
    }		
  }
  return true;
}

//calculate hamming distance for given pair of row
inline int SprIndicatorMatrix::rowHammingDistance(int row1, int row2) const
{
  int result = 0;
  int numcol = this->num_col();
  for(int j=0;j<numcol;j++){
    int value1 = this->get(row1,j);
    int value2 = this->get(row2,j);
    if(value1==0 || value2==0){continue;}
    result += (value1==value2)?0:1;
  }
  return result;
}

//calculate hamming distance for given pair of columns
inline int SprIndicatorMatrix::columnHammingDistance(int col1, int col2) const
{
  int result = 0;
  int resultComplement = 0;//HD for treating col1 and col2 as complement
  int numrow = this->num_row();
  for(int i=0;i<numrow;i++){
    int value1 = this->get(i,col1);
    int value2 = this->get(i,col2);
    if(value1==0 || value2==0){continue;}
    if(value1==value2){
      result++;
    }
    else{
      resultComplement++;
    }
  }
  return ( result<resultComplement ? result : resultComplement );
}

//find minimum pairwise row hamming distance for the matrix
int SprIndicatorMatrix::minRowHammingDistance() const
{
  int numcol=this->num_col();
  int numrow=this->num_row();
  int currentmin=numcol;
  for(int row1=0;row1<numrow;row1++){
    for(int row2=row1+1;row2<numrow;row2++){
      int hammingDistance = this->rowHammingDistance(row1,row2);
      if(hammingDistance < currentmin){
        currentmin = hammingDistance;
      }
    }
  }
  return currentmin;
}

//find minimum column hamming distance for the matrix
int SprIndicatorMatrix::minColumnHammingDistance() const
{
  int numcol=this->num_col();
  int numrow=this->num_row();
  int currentmin=numrow;
  for(int col1=0;col1<numcol-1;col1++){
    for(int col2=col1+1;col2<numcol;col2++){
      int hammingDistance = this->columnHammingDistance(col1,col2);
      if(hammingDistance < currentmin){
        currentmin = hammingDistance;
      }
    }
  }
  return currentmin;
}

double SprIndicatorMatrix::hammingMeasure() const 
{
  return (this->minRowHammingDistance()+this->minColumnHammingDistance());
}

double SprIndicatorMatrix::minRowHammingMeasure() const
{
  return this->minRowHammingDistance();
}

double SprIndicatorMatrix::minColumnHammingMeasure() const
{
  return this->minColumnHammingDistance();
}

double SprIndicatorMatrix::diversityMeasure() const
{
  return (this->rowDiversity()+this->columnDiversity());
}


//find sum pairwise row hamming distance for the matrix
//the result is normalized by number of row pairs (row*(row-1)/2)
double SprIndicatorMatrix::rowDiversity() const
{
  //int numcol=this->num_col();
  int numrow=this->num_row();
  int currentsum=0;
  for(int row1=0;row1<numrow-1;row1++){
    for(int row2=row1+1;row2<numrow;row2++){
      int hammingDistance = this->rowHammingDistance(row1,row2);
      currentsum += hammingDistance;
    }
  }
  return double(currentsum)/double(numrow*(numrow-1)/2);
}

//find sum pairwise column hamming distance for the matrix
//the result is normalized by number of col pairs (col*(col-1)/2)
double SprIndicatorMatrix::columnDiversity() const
{
  int numcol=this->num_col();
  //int numrow=this->num_row();
  int currentsum=0;
  for(int col1=0;col1<numcol-1;col1++){
    for(int col2=col1+1;col2<numcol;col2++){
      int hammingDistance = this->columnHammingDistance(col1,col2);
      currentsum += hammingDistance;
    }
  }
  return double(currentsum)/double(numcol*(numcol-1)/2);
}

SprMatrix* SprIndicatorMatrix::toSprMatrix() const
{
  int numrow = this->num_row();
  int numcol = this->num_col();
  SprMatrix* toReturn = new SprMatrix(numrow, numcol);
  for(int i=0;i<numrow;i++){
    for(int j=0;j<numcol;j++){
      (*toReturn)[i][j]=this->get(i,j);
    }
  }
  return toReturn;
}


double SprIndicatorMatrix::evaluate(MatrixMeasure measure) const
{
  switch(measure)
    {
    case MINROW:
      return this->minRowHammingMeasure();
      break;
    case HAMMING:
      return this->hammingMeasure();
      break;
    case DIVERSITY:
      return this->diversityMeasure();
      break;
    }
  return 0;
}


void SprIndicatorMatrix::print(std::ostream& os) const 
{
  for(int i=0;i<nrow_;i++){
    for(int j=0;j<ncol_;j++){
      const char* toprint 
	= ( m_[this->toVectorPosition(i,j)]==0 ? "  0" 
	: ( m_[this->toVectorPosition(i,j)]==-1 ? " -1" : "  1" ) );
      os << toprint;
    }
    os << endl;
  }
}
