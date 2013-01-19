//$Id: SprIndicatorMatrix.hh,v 1.2 2008-05-09 21:25:26 narsky Exp $

/*************************
 *
 * Author: Piti Ongmongkolkul
 * Adapted from SprMatrix
 * Reference: Reducing Multiclass to Binary: 
 *  A Unifying Approach for Margin Classifiers
 *            Erin L. Allwein, Robert E Schapire, Yoram Singer: 
 *              Journal of Learning Machine (2000)
 *
 */

#ifndef _SprIndicatorMatrix_HH
#define _SprIndicatorMatrix_HH

#include "StatPatternRecognition/SprRandomNumber.hh"
#include "StatPatternRecognition/SprMatrix.hh"

#include <vector>
#include <utility>
#include <iostream>

class SprIndicatorMatrix {

public:
  ~SprIndicatorMatrix(){ if(m_!=0) delete[] m_; }

  //constructor
  SprIndicatorMatrix(int p, int q);

  //copy constructor
  SprIndicatorMatrix(const SprIndicatorMatrix& matrix);

  //convert from SprMatrix to SprIndicatorMatrix
  SprIndicatorMatrix(const SprMatrix& matrix);

  // enum for measure use in random algorithm
  // MINROW will use minimum row hamming distance
  // HAMMING will use sum of min row hamming distance and min col hamming dist
  // DIVERSITY will use sum of average row hamming distance 
  //    and sum of average col hamming distance
  enum MatrixMeasure {MINROW, HAMMING, DIVERSITY};

  // randomly fill matrix with -1(background) and 1(signal) 
  //   with the given probabilities
  // numtrial times and check for duplicated rows/columns 
  //   and no complement column (?)
  // if none of the generated matrices passes the check 
  //    then this function returns 0
  // user is responsible for deleting this object
  static SprIndicatorMatrix* randomDense(
  	int row, int column, 
  	double psignal=0.5, double pbackground=0.5, 
  	int numtrial=1000, 
  	MatrixMeasure measure=MINROW);
  
  static SprIndicatorMatrix* randomSparse(
  	int row, int column, 
  	double psignal=0.25, double pbackground=0.25, double pignore=0.5, 
  	int numtrial=1000, 
  	MatrixMeasure measure=MINROW);

  static SprIndicatorMatrix* randomHillClimbing(
	int row, int column, int numstep=1000, double pKeepBadChange=0.0,
  	MatrixMeasure measure=MINROW);
  
  
  //return exhaustive code matrix (guaranteed to be optimal)
  //warning: the number of column is 2^(n-1)-1 
  //    so don't try anything bigger than 10 
  //    unless you have plenty of computing power
  //user is responsible for deleting this object
  static SprIndicatorMatrix* exhaustive(int row);
  
  //matrix of One Vs One append to One Vs All
  //this kind of matrix will pick only easy binary problems
  static SprIndicatorMatrix* ovoova(int row);
  
  //check for duplicated row/column and complement column
  //return true when no such thing is found
  bool checkMatrix() const;
  
  //calculate hamming distance for given rows/column
  int rowHammingDistance(int row1, int row2) const;
  int columnHammingDistance(int col1, int col2) const;
	
  //return minimum pairwise hamming distance 
  //  of the given codeMatrix (O(numrow^2*numcolumn))
  int minRowHammingDistance() const;
  int minColumnHammingDistance() const;
  
  //calculate diversity
  //see Using diversity measures for generating error-correcting output codes 
  //  in classifier ensembles Ludmila I. Kuncheva
  double rowDiversity() const;
  double columnDiversity() const;
  
  // Vector of indices for equal elements
  std::vector<int> equalColsForRows(int row1, int row2) const;
  std::vector<int> equalRowsForCols(int col1, int col2) const;
  
  //return vector of pair position( vector of size 2 of int) where 
  //row/col hamming distance is lowest
  std::vector<std::pair<int,int> > closestRowPairs();
  std::vector<std::pair<int,int> > closestColPairs();
  
  //measures for selecting best matrix
  double evaluate(MatrixMeasure measure) const;
  double diversityMeasure() const;
  double hammingMeasure() const;
  double minRowHammingMeasure() const;
  double minColumnHammingMeasure() const;
  
  inline int num_row() const;
  inline int num_col() const;

  //accessor
  inline int get(int row,int col) const;
  inline void set(int row, int col, int value);
  void swap(int row, int col);//switch 1 to -1 (leave zero)
  void print(std::ostream& os) const;
  //  std::string toString() const;
  SprMatrix* toSprMatrix() const;

private:
  static const double tolerance_;
  static const int maxExhaustiveRows_;
  static SprRandomNumber rndm_;

  int nrow_; 
  int ncol_; 
  int size_;
  int* m_;

  inline int toVectorPosition(int row, int col) const;
};

#endif

inline int SprIndicatorMatrix::num_row() const { return nrow_;}
inline int SprIndicatorMatrix::num_col() const  { return ncol_;}

inline int SprIndicatorMatrix::toVectorPosition(int row, int col) const { 
  return ncol_*row+col;
}

inline int SprIndicatorMatrix::get(int row, int col) const { 
  return m_[this->toVectorPosition(row,col)];
}

inline void SprIndicatorMatrix::set(int row, int col,int value) { 
  m_[this->toVectorPosition(row,col)]=value;
}
