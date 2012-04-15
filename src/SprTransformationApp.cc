//$Id: SprTransformationApp.cc,v 1.4 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprData.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprAbsReader.hh"
#include "StatPatternRecognition/SprAbsWriter.hh"
#include "StatPatternRecognition/SprDataFeeder.hh"
#include "StatPatternRecognition/SprRWFactory.hh"
#include "StatPatternRecognition/SprStringParser.hh"
#include "StatPatternRecognition/SprAbsVarTransformer.hh"
#include "StatPatternRecognition/SprPCATransformer.hh"
#include "StatPatternRecognition/SprInputNormalizer.hh"
#include "StatPatternRecognition/SprReplaceMissing.hh"
#include "StatPatternRecognition/SprVarTransformerSequence.hh"
#include "StatPatternRecognition/SprTransformerFilter.hh"

#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <memory>

using namespace std;


void help(const char* prog) 
{
  cout << "Usage:  " << prog 
       << " sequence_of_filters_to_be_trained"
       << " training_data_file output_file_for_storing_transformation"  
       << endl;
  cout << "\t Filters must be given in quotes, separated by commas,"
       << " e.g., \'PCA,Normalize\' would imply that first  "   << endl;
  cout << "\t a PCA transformation is applied and then all inputs"
       << " are normalized."                                    << endl;
  cout << "\t List of available transformations:                " << endl;
  cout << "\t PCA - Principal Component Analysis                " << endl;
  cout << "\t Normalize - X -> (X-<X>)/Sx                       " << endl;
  cout << "\t Missing - Find replacements for missing values    " << endl;
  // 
  cout << "\t Options: " << endl;
  cout << "\t-h --- help                                        " << endl;
  cout << "\t-o output Tuple file                               " << endl;
  cout << "\t-a input ascii file mode (see SprSimpleReader.hh)  " << endl;
  cout << "\t-A save output data in ascii instead of Root       " << endl;
  cout << "\t-y list of input classes (see SprAbsFilter.hh)     " << endl;
  cout << "\t-m mode for Missing transformation:                " << endl;
  cout << "\t\t 1 - Median (default)                            " << endl;
  cout << "\t\t 2 - Average                                     " << endl;
  cout << "\t-L lower cutoff for Missing transformation.        " << endl;
  cout << "\t\t Every value below this number is treated as missing." << endl;
  cout << "\t-v verbose level (0=silent default,1,2)            " << endl;
  cout << "\t-V include only these input variables              " << endl;
  cout << "\t-z exclude input variables from the list           " << endl;
  cout << "\t\t Variables must be listed in quotes and separated by commas." 
       << endl;
}


int main(int argc, char ** argv)
{
  // check command line
  if( argc < 2 ) {
    help(argv[0]);
    return 1;
  }

  // init
  string tupleFile;
  int readMode = 0;
  SprRWFactory::DataType writeMode = SprRWFactory::Root;
  int verbose = 0;
  string outFile;
  string includeList, excludeList;
  string inputClassesString;
  string stringVarsDoNotFeed;
  bool setLowCutoff = false;
  double lowCutoff = 0;
  int missingMode = 1;
  bool classBlind = true;

  // decode command line
  int c;
  extern char* optarg;
  //  extern int optind;
  while( (c = getopt(argc,argv,"ho:a:Ay:m:bL:v:K:DV:z:Z:")) != EOF ) {
    switch( c )
      {
      case 'h' :
	help(argv[0]);
	return 1;
      case 'o' :
	tupleFile = optarg;
	break;
      case 'a' :
	readMode = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'A' :
	writeMode = SprRWFactory::Ascii;
	break;
      case 'y' :
	inputClassesString = optarg;
	break;
      case 'm' :
	missingMode = (optarg==0 ? 1 : atoi(optarg));
	break;
      case 'b' :
	classBlind = false;
	break;
      case 'L' :
	if( optarg != 0 ) {
	  setLowCutoff = true;
	  lowCutoff = atof(optarg);
	}
	break;
      case 'v' :
	verbose = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'V' :
	includeList = optarg;
	break;
      case 'z' :
	excludeList = optarg;
	break;
      case 'Z' :
	stringVarsDoNotFeed = optarg;
	break;
      }
  }

  // 2 arguments after all options.
  string filterString = argv[argc-3];
  string trFile = argv[argc-2];
  string transformationFile = argv[argc-1];
  if( trFile.empty() ) {
    cerr << "No training file is specified." << endl;
    return 1;
  }
  if( transformationFile.empty() ) {
    cerr << "No file for storing the trained transformation is specified." 
	 << endl;
    return 1;
  }
  if( filterString.empty() ) {
    cerr << "No transformations specified." << endl;
    return 1;
  }

  // make reader
  SprRWFactory::DataType inputType 
    = ( readMode==0 ? SprRWFactory::Root : SprRWFactory::Ascii );
  auto_ptr<SprAbsReader> reader(SprRWFactory::makeReader(inputType,readMode));

  // include variables
  set<string> includeSet;
  if( !includeList.empty() ) {
    vector<vector<string> > includeVars;
    SprStringParser::parseToStrings(includeList.c_str(),includeVars);
    assert( !includeVars.empty() );
    for( int i=0;i<includeVars[0].size();i++ ) 
      includeSet.insert(includeVars[0][i]);
    if( !reader->chooseVars(includeSet) ) {
      cerr << "Unable to include variables in training set." << endl;
      return 2;
    }
    else {
      cout << "Following variables have been included in optimization: ";
      for( set<string>::const_iterator 
	     i=includeSet.begin();i!=includeSet.end();i++ )
	cout << "\"" << *i << "\"" << " ";
      cout << endl;
    }
  }

  // exclude variables
  set<string> excludeSet;
  if( !excludeList.empty() ) {
    vector<vector<string> > excludeVars;
    SprStringParser::parseToStrings(excludeList.c_str(),excludeVars);
    assert( !excludeVars.empty() );
    for( int i=0;i<excludeVars[0].size();i++ ) 
      excludeSet.insert(excludeVars[0][i]);
    if( !reader->chooseAllBut(excludeSet) ) {
      cerr << "Unable to exclude variables from training set." << endl;
      return 2;
    }
    else {
      cout << "Following variables have been excluded from optimization: ";
      for( set<string>::const_iterator 
	     i=excludeSet.begin();i!=excludeSet.end();i++ )
	cout << "\"" << *i << "\"" << " ";
      cout << endl;
    }
  }

  // read training data from file
  auto_ptr<SprAbsFilter> filter(reader->read(trFile.c_str()));
  if( filter.get() == 0 ) {
    cerr << "Unable to read data from file " << trFile.c_str() << endl;
    return 2;
  }
  vector<string> vars;
  filter->vars(vars);
  cout << "Read data from file " << trFile.c_str() 
       << " for variables";
  for( int i=0;i<vars.size();i++ ) 
    cout << " \"" << vars[i].c_str() << "\"";
  cout << endl;
  cout << "Total number of points read: " << filter->size() << endl;

  // filter training data by class
  vector<SprClass> inputClasses;
  if( !filter->filterByClass(inputClassesString.c_str()) ) {
    cerr << "Cannot choose input classes for string " 
	 << inputClassesString << endl;
    return 2;
  }
  filter->classes(inputClasses);
  assert( inputClasses.size() > 1 );
  cout << "Training data filtered by class." << endl;
  for( int i=0;i<inputClasses.size();i++ ) {
    cout << "Points in class " << inputClasses[i] << ":   " 
	 << filter->ptsInClass(inputClasses[i]) << endl;
  }

  // decode transformations
  vector<vector<string> > transformerNames;
  SprStringParser::parseToStrings(filterString.c_str(),transformerNames);
  if( transformerNames.empty() || transformerNames[0].empty() ) {
    cerr << "Unable to decode list of transformations." << endl;
    return 1;
  }

  // book transformations
  vector<pair<SprAbsVarTransformer*,bool> > transformers;
  for( int i=0;i<transformerNames[0].size();i++ ) {
    if(      transformerNames[0][i] == "PCA" ) {
      transformers.push_back(pair<SprAbsVarTransformer*,
			     bool>(new SprPCATransformer(),true));
    }
    else if( transformerNames[0][i] == "Normalize" ) {
      transformers.push_back(pair<SprAbsVarTransformer*,
			     bool>(new SprInputNormalizer(),true));
    }
    else if( transformerNames[0][i] == "Missing" ) {
      if( !setLowCutoff) {
	cerr << "Low cutoff has not been set to identify missing values." 
	     << endl;
	return 2;
      }
      SprReplaceMissing::Mode mMode = ( missingMode==2 ? 
					SprReplaceMissing::Average : 
					SprReplaceMissing::Median );
      SprReplaceMissing* t 
	= new SprReplaceMissing(mMode,
				SprUtils::lowerBound(lowCutoff),
				classBlind);
      transformers.push_back(pair<SprAbsVarTransformer*,bool>(t,true));
    }
    else {
      cerr << "Unknown transformation requested: " 
	   << transformerNames[0][i].c_str() << endl;
      return 2;
    }
  }
  assert( !transformers.empty() );

  // make transformer sequence
  SprVarTransformerSequence seq(transformers);

  // train
  if( !seq.train(filter.get(),verbose) ) {
    cerr << "Unable to train transformation sequence." << endl;
    return 3;
  }
  cout << "Transformation sequence has been trained." << endl;

  // store
  if( !seq.store(transformationFile.c_str()) ) {
    cerr << "Unable to store transformation coefficients." << endl;
    return 4;
  }
  cout << "Stored the trained transformation" << endl;

  // make histogram if requested
  if( tupleFile.empty() ) return 0;

  // make a writer
  auto_ptr<SprAbsWriter> tuple(SprRWFactory::makeWriter(writeMode,"training"));
  if( !tuple->init(tupleFile.c_str()) ) {
    cerr << "Unable to open output file " << tupleFile.c_str() << endl;
    return 5;
  }

  // transform data
  SprTransformerFilter trans(filter.get());
  bool replaceOriginalData = true;
  if( !trans.transform(&seq,replaceOriginalData) ) {
    cerr << "Unable to transform input data." << endl;
    return 6;
  }

  // feed
  SprDataFeeder feeder(&trans,tuple.get());
  if( !feeder.feed(1000) ) {
    cerr << "Cannot feed data into file " << tupleFile.c_str() << endl;
    return 9;
  }

  // exit
  return 0;
}
