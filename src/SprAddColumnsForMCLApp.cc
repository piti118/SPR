//$Id: SprAddColumnsForMCLApp.cc,v 1.1 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprMultiClassReader.hh"
#include "StatPatternRecognition/SprTrainedMultiClassLearner.hh"

#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <memory>

using namespace std;


void help(const char* prog) 
{
  cout << "Usage:  " << prog 
       << " input_file_with_a_list_of_binary_classifiers" 
       << " output_file_with_trained_multiclass_learner" 
       << endl;
  cout << "\t Options: " << endl;
  cout << "\t-h --- help                                             " << endl;
  cout << "\t-v verbose level (0=silent default,1,2)                 " << endl;
}


int main(int argc, char ** argv)
{
  // check command line
  if( argc < 3 ) {
    help(argv[0]);
    return 1;
  }

  // init
  string inputFile, outputFile;
  int verbose = 0;

  // decode command line
  int c;
  extern char* optarg;
  //  extern int optind;
  while( (c = getopt(argc,argv,"hv:")) != EOF ) {
    switch( c )
      {
      case 'h' :
	help(argv[0]);
	return 1;
      case 'v' :
	verbose = (optarg==0 ? 0 : atoi(optarg));
	break;
      }
  }

  // get file names
  inputFile = argv[argc-2];
  outputFile = argv[argc-1];
  if( inputFile.empty() ) {
    cerr << "Input file not specified." << endl;
    return 1;
  }
  if( outputFile.empty() ) {
    cerr << "Output file not specified." << endl;
    return 1;
  }

  // read trained MC learner
  auto_ptr<const SprTrainedMultiClassLearner> 
    trained(SprMultiClassReader::readBinaryList(inputFile.c_str(),verbose));
  if( trained.get() == 0 ) {
    cerr << "Unable to read trained multiclass learner from file " 
	 << inputFile.c_str() << endl;
    return 2;
  }
  cout << "Trained binary classifiers have been read." << endl;

  // save trained MC learner
  cout << "Storing multiclass learner." << endl;
  if( !trained->store(outputFile.c_str()) ) {
    cerr << "Unable to store trained learner into file "
	 << outputFile.c_str() << endl;
    return 3;
  }
  cout << "Multiclass learner has been saved to " 
       << outputFile.c_str() << endl;

  // exit
  return 0;
}
