//$Id: SprIndicatorMatrixApp.cc,v 1.1 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprIndicatorMatrix.hh"
#include "StatPatternRecognition/SprConfig.hh"

#include <iostream>
#include <string>
#include <unistd.h>

using namespace std;

//all program constant
namespace SprIndicatorMatrixApp {
  const string algorithmKey="algorithm";
  const string defaultAlgorithm="random";
  
  const string numClassKey="class";
  const string numClassifierKey="column";
  
  const string probSignalKey="probsignal";
  const string probBackgroundKey="probbackground";
  const string probIgnoreKey="probignore";
  const double defaultProbSignal=1.0;
  const double defaultProbBackground=1.0;
  const double defaultProbIgnore=0.0;
  
  const string numTrialKey="numtrial";
  const int defaultNumTrial = 1000;
  
  const string randomAlgoritm="random";
  const string exhaustiveAlgorithm="exhaustive";
  const string ovoovaAlgorithm="ovoova";
  const string randomHillAlgorithm="randomhill";
  
  const string measureKey="measure";
  const string minRowHammingMeasure="minrow";
  const string hammingMeasure="hamming";
  const string diversityMeasure="diversity";
  const string defaultMeasure="minrow";

  const string probKeepBadChangeKey="probKeepBadChange";
  const double defaultProbKeepBadChange=0.0;
}


SprIndicatorMatrix::MatrixMeasure measureString2measureEnum(string& s)
{
  namespace appconst=SprIndicatorMatrixApp;
  if(s==appconst::minRowHammingMeasure){
    return SprIndicatorMatrix::MINROW;
  }
  else if(s==appconst::hammingMeasure){
    return SprIndicatorMatrix::HAMMING;
  }
  else if(s==appconst::diversityMeasure){
    return SprIndicatorMatrix::DIVERSITY;
  }
  //make sure this coincide with the default one
  return SprIndicatorMatrix::MINROW;
}

void help(const char* prog) 
{
  cout << "Usage:  " << prog 
       << " configfile (see data/samplematrixindicator.cfg)"
       << endl;
}

int main(int argc, char ** argv)
{
  // check command line
  if( argc < 2 ) {
    help(argv[0]);
    return 1;
  }

  // decode
  char c;
  namespace appconst=SprIndicatorMatrixApp;
  while( (c = getopt(argc,argv,"h")) != EOF ) {
    switch( c )
      {
      case 'h' :
	help(argv[0]);
	return 0;
      	break;
      }
  }

  // get config file
  string configfile = argv[argc-1];
  if(configfile.empty()){
    cerr << "Please supply the config filename" << endl;
    return 1;
  }
  
  //now read the configfile
  SprConfig config(configfile);
  
  if(!config.keyExists(appconst::numClassKey)){
    cerr << "number of class is not specified "
	 << "in the config file (ex: class=10)" << endl;
    return 1;
  }
  
    //check the configuration key
  string algorithm 
    = config.getStringValue(appconst::algorithmKey,appconst::randomAlgoritm);
  int numclass = config.getIntValue(appconst::numClassKey);
  int numclassifier = config.getIntValue(appconst::numClassifierKey);
  double probsignal 
    = config.getDoubleValue(appconst::probSignalKey,
			    appconst::defaultProbSignal);
  double probbackground 
    = config.getDoubleValue(appconst::probBackgroundKey,
			    appconst::defaultProbBackground);
  double probignore 
    = config.getDoubleValue(appconst::probIgnoreKey,
			    appconst::defaultProbIgnore);
  int numTrial 
    = config.getIntValue(appconst::numTrialKey,appconst::defaultNumTrial);
  string measurevalue 
    = config.getStringValue(appconst::measureKey,appconst::defaultMeasure);
  
  if(algorithm == appconst::randomAlgoritm){
    //check if number of column is specified
    if(!config.keyExists(appconst::numClassifierKey)){
      cerr << "number of classifiers is not specified "
	   << "for random algorithm. (ex: column=10)" << endl;
      return 1;
    }
    SprIndicatorMatrix::MatrixMeasure measure 
      = measureString2measureEnum(measurevalue);
    SprIndicatorMatrix* matrix 
      = SprIndicatorMatrix::randomSparse( numclass, 
					  numclassifier, 
					  probsignal, 
					  probbackground, 
					  probignore, 
					  numTrial,
					  measure);
    
    if(matrix==0){//can't generate matrix
      cerr << "Sorry, none of the generated matrices is valid"<<endl;
      return 1;
    }
    else{
      cout << "# Generated from SprIndicatorMatrixApp using random algorithm" 
	   << endl;
      cout << "# Optimizing on " << measurevalue << endl;
      cout << "# Minimal pairwise row Hamming distance: " 
	   << (matrix->minRowHammingMeasure()) << endl;
      cout << "# Minimal pairwise column Hamming distance: " 
	   << (matrix->minColumnHammingMeasure()) << endl;
      cout << "# RowDiversity: " << matrix->rowDiversity() << endl;
      cout << "# ColumnDiversity: " << matrix->columnDiversity() << endl;
      cout << "# MatrixDiversity: " << matrix->diversityMeasure() << endl;
      cout << "# HammingMeasure: "<< matrix->hammingMeasure() << endl;
      cout << "# ProbSignal: " << probsignal << endl;
      cout << "# ProbBackground: " << probbackground << endl;
      cout << "# ProbIgnore: " << probignore << endl;
      cout << "# numTrial: " << numTrial << endl;
      cout << numclass << " " << numclassifier << endl;
      matrix->print(cout);
      delete matrix;
    } 
  }
  else if(algorithm==appconst::randomHillAlgorithm){
    double pKeepBadChange
      = config.getDoubleValue(appconst::probKeepBadChangeKey,
			      appconst::defaultProbKeepBadChange);
    if(!config.keyExists(appconst::numClassifierKey)){
      cerr << "number of classifier is not specified "
	   << "for random algorithm. (ex: column=10)" << endl;
      return 1;
    }
    SprIndicatorMatrix::MatrixMeasure measure 
      = measureString2measureEnum(measurevalue);
    SprIndicatorMatrix* matrix 
      = SprIndicatorMatrix::randomHillClimbing(numclass,
					       numclassifier,
					       numTrial,
					       pKeepBadChange,
					       measure);
    if(matrix==0){
      cerr<< "Cannot generate matrix" << endl;
      return 1;
    }
    else{
      cout << "# Generated from SprIndicatorMatrixApp "
	   << " using random hill climbing algorithm" << endl;
      cout << "# Optimizing on " << measurevalue << endl;
      cout << "# Minimal pairwise row Hamming distance: " 
	   << (matrix->minRowHammingMeasure()) << endl;
      cout << "# Minimal pairwise column Hamming distance: " 
	   << (matrix->minColumnHammingMeasure()) << endl;
      cout << "# RowDiversity: " << matrix->rowDiversity() << endl;
      cout << "# ColumnDiversity: " << matrix->columnDiversity() << endl;
      cout << "# MatrixDiversity: " << matrix->diversityMeasure() << endl;
      cout << "# HammingMeasure: "<< matrix->hammingMeasure() << endl;
      cout << "# ProbSignal: " << probsignal << endl;
      cout << "# ProbBackground: " << probbackground << endl;
      cout << "# ProbIgnore: " << probignore << endl;
      cout << "# numTrial: " << numTrial << endl;
      cout << numclass << " " << numclassifier << endl;
      matrix->print(cout);
      delete matrix;
    }
  }
  else if(algorithm == appconst::exhaustiveAlgorithm){
    SprIndicatorMatrix* matrix = SprIndicatorMatrix::exhaustive(numclass);
    if(matrix==0){
      cerr << "Cannot generate exhaustive matrix "
	   <<" (try something lower than 30)" << endl;
      return 1;
    }
    else{
      cout << "# Generated from SprIndicatorMatrixApp "
	   << "using exhaustive algorithm" << endl;
      cout << "# Minimal pairwise row Hamming distance: " 
	   << (matrix->minRowHammingMeasure()) << endl;
      cout << "# Minimal pairwise column Hamming distance: " 
	   << (matrix->minColumnHammingMeasure()) << endl;
      cout << "# RowDiversity: " << matrix->rowDiversity() << endl;
      cout << "# ColumnDiversity: " << matrix->columnDiversity() << endl;
      cout << "# MatrixDiversity: " << matrix->diversityMeasure() << endl;
      cout << "# HammingMeasure: "<< matrix->hammingMeasure() << endl;
      cout << numclass << " " << matrix->num_col() << endl;
      matrix->print(cout);
      delete matrix;
    }
  }
  else if(algorithm==appconst::ovoovaAlgorithm){
    SprIndicatorMatrix* matrix = SprIndicatorMatrix::ovoova(numclass);
    if(matrix==0){
      cerr<< "Cannot generate ovoova matrix" << endl;
      return 1;
    }
    else{
      cout << "# Generated from SprIndicatorMatrixApp using ovoova algorithm" 
	   << endl;
      cout << "# Minimal pairwise row Hamming distance: " 
	   << (matrix->minRowHammingMeasure()) << endl;
      cout << "# Minimal pairwise column Hamming distance: " 
	   << (matrix->minColumnHammingMeasure()) << endl;
      cout << "# RowDiversity: " << matrix->rowDiversity() << endl;
      cout << "# ColumnDiversity: " << matrix->columnDiversity() << endl;
      cout << "# MatrixDiversity: " << matrix->diversityMeasure() << endl;
      cout << "# HammingMeasure: "<< matrix->hammingMeasure() << endl;
      cout << numclass << " " << matrix->num_col() << endl;
      matrix->print(cout);
      delete matrix;
    }
  }
  
  else{
    cerr << "Unknown algorithm \"" << algorithm << "\"."<< endl;
    return 1;
    //unknown algorithm
  }
  return 0;  
}

