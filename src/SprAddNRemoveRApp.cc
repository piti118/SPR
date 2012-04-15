//$Id: SprAddNRemoveRApp.cc,v 1.5 2008-04-02 23:36:45 narsky Exp $
//
// An executable for variable selction using "add n remove r" strategy.
// See notes in README (Variable Selection).
//

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprAbsClassifier.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprEmptyFilter.hh"
#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprAbsReader.hh"
#include "StatPatternRecognition/SprRWFactory.hh"
#include "StatPatternRecognition/SprStringParser.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprClassifierReader.hh"
#include "StatPatternRecognition/SprClassifierEvaluator.hh"
#include "StatPatternRecognition/SprAbsVarTransformer.hh"
#include "StatPatternRecognition/SprVarTransformerReader.hh"
#include "StatPatternRecognition/SprTransformerFilter.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprAbsTwoClassCriterion.hh"
#include "StatPatternRecognition/SprIntegerBootstrap.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprLoss.hh"
#include "StatPatternRecognition/SprAverageLoss.hh"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <memory>
#include <string>
#include <cassert>

using namespace std;


void prepareExit(vector<SprAbsTwoClassCriterion*>& criteria,
		 vector<SprAbsClassifier*>& classifiers,
		 vector<SprIntegerBootstrap*>& bstraps) 
{
  for( int i=0;i<criteria.size();i++ ) delete criteria[i];
  for( int i=0;i<classifiers.size();i++ ) delete classifiers[i];
  for( int i=0;i<bstraps.size();i++ ) delete bstraps[i];
}


void help(const char* prog) 
{
  cout << "Usage:  " << prog
       << " file_of_classifier_parameters training_data_file" << endl;
  cout << "\t Options: " << endl;
  cout << "\t-h --- help                                        " << endl;
  cout << "\t-y list of input classes (see SprAbsFilter.hh)     " << endl;
  cout << "\t-Q apply variable transformation saved in file     " << endl;
  cout << "\t-a input ascii file mode (see SprSimpleReader.hh)  " << endl;
  cout << "\t-n add n variables at a time (default=1)           "<< endl;
  cout << "\t-r remove r variables at a time (default=0)        "<< endl;
  cout << "\t-g per-event loss                                  " << endl;
  cout << "\t\t 1 - misidentified fraction (default)            " << endl;
  cout << "\t\t 2 - quadratic loss (y-f(x))^2                   " << endl;
  cout << "\t-t test data (must be in same format as training data)"<< endl;
  cout << "\t-x use that many pieces for cross-validation "       << endl;
  cout << "\t\t if no test data is supplied.                    " << endl;
  cout << "\t-i integrate classifier responses over CV data (see README)"
       << endl;
  cout << "\t-v verbose level (0=silent default,1,2)            " << endl;
  cout << "\t-w scale all signal weights by this factor         " << endl;
  cout << "\t-V include only these input variables              " << endl;
  cout << "\t-z exclude input variables from the list           " << endl;
}


int main(int argc, char ** argv)
{
  // check command line
  if( argc < 3 ) {
    help(argv[0]);
    return 1;
  }

  // init
  int readMode = 0;
  int verbose = 0;
  bool scaleWeights = false;
  double sW = 1.;
  string includeList, excludeList;
  string inputClassesString;
  string transformerFile;
  unsigned N = 1;
  unsigned R = 0;
  int iLoss = 0;
  string valFile;
  unsigned nCross = 0;
  bool integrate = false;

  // decode command line
  int c;
  extern char* optarg;
  extern int optind;
  while( (c = getopt(argc,argv,"hy:Q:a:n:r:g:t:x:iv:w:V:z:")) != EOF ) {
    switch( c )
      {
      case 'h' :
	help(argv[0]);
	return 1;
      case 'y' :
	inputClassesString = optarg;
	break;
      case 'Q' :
        transformerFile = optarg;
        break;
      case 'a' :
	readMode = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'n' :
	N = (optarg==0 ? 1 : atoi(optarg));
	break;
      case 'r' :
	R = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'g' :
        iLoss = (optarg==0 ? 1 : atoi(optarg));
        break;
      case 't' :
	valFile = optarg;
	break;
      case 'x' :
	nCross = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'i' :
	integrate = true;
	break;
      case 'v' :
	verbose = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'w' :
	if( optarg != 0 ) {
	  scaleWeights = true;
	  sW = atof(optarg);
	}
	break;
      case 'V' :
	includeList = optarg;
	break;
      case 'z' :
	excludeList = optarg;
	break;
      }
  }

  // Must have 3 arguments on the command line
  string configFile     = argv[argc-2];
  string dataFile       = argv[argc-1];
  if( configFile.empty() ) {
    cerr << "No classifier configuration file is specified." << endl;
    return 1;
  }
  if( dataFile.empty() ) {
    cerr << "No input data file is specified." << endl;
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
  auto_ptr<SprAbsFilter> filter(reader->read(dataFile.c_str()));
  if( filter.get() == 0 ) {
    cerr << "Unable to read data from file " << dataFile.c_str() << endl;
    return 2;
  }
  vector<string> vars;
  filter->vars(vars);
  cout << "Read data from file " << dataFile.c_str() << " for variables";
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

  // scale weights
  if( scaleWeights ) {
    cout << "Signal weights are multiplied by " << sW << endl;
    filter->scaleWeights(inputClasses[1],sW);
  }

  // read validation data from file
  auto_ptr<SprAbsFilter> valFilter;
  if( !valFile.empty() ) {
    auto_ptr<SprAbsReader> 
      valReader(SprRWFactory::makeReader(inputType,readMode));
    if( !includeSet.empty() ) {
      if( !valReader->chooseVars(includeSet) ) {
	cerr << "Unable to include variables in validation set." << endl;
	return 2;
      }
    }
    if( !excludeSet.empty() ) {
      if( !valReader->chooseAllBut(excludeSet) ) {
	cerr << "Unable to exclude variables from validation set." << endl;
	return 2;
      }
    }
    valFilter.reset(valReader->read(valFile.c_str()));
    if( valFilter.get() == 0 ) {
      cerr << "Unable to read data from file " << valFile.c_str() << endl;
      return 2;
    }
    vector<string> valVars;
    valFilter->vars(valVars);
    cout << "Read validation data from file " << valFile.c_str() 
	 << " for variables";
    for( int i=0;i<valVars.size();i++ ) 
      cout << " \"" << valVars[i].c_str() << "\"";
    cout << endl;
    cout << "Total number of points read: " << valFilter->size() << endl;
    cout << "Points in class 0: " << valFilter->ptsInClass(inputClasses[0])
	 << " 1: " << valFilter->ptsInClass(inputClasses[1]) << endl;
  }

  // filter validation data by class
  if( valFilter.get() != 0 ) {
    if( !valFilter->filterByClass(inputClassesString.c_str()) ) {
      cerr << "Cannot choose input classes for string " 
	   << inputClassesString << endl;
      return 2;
    }
    valFilter->classes(inputClasses);
    cout << "Validation data filtered by class." << endl;
    for( int i=0;i<inputClasses.size();i++ ) {
      cout << "Points in class " << inputClasses[i] << ":   " 
	   << valFilter->ptsInClass(inputClasses[i]) << endl;
    }
  }

  // scale weights
  if( scaleWeights && valFilter.get()!=0 )
    valFilter->scaleWeights(inputClasses[1],sW);

  // apply transformation of variables to training and test data
  auto_ptr<SprAbsFilter> garbage_train, garbage_valid;
  if( !transformerFile.empty() ) {
    const SprAbsVarTransformer* t 
      = SprVarTransformerReader::read(transformerFile.c_str());
    if( t == 0 ) {
      cerr << "Unable to read VarTransformer from file "
           << transformerFile.c_str() << endl;
      return 2;
    }
    SprTransformerFilter* t_train = new SprTransformerFilter(filter.get());
    SprTransformerFilter* t_valid = 0;
    if( valFilter.get() != 0 )
      t_valid = new SprTransformerFilter(valFilter.get());
    bool replaceOriginalData = true;
    if( !t_train->transform(t,replaceOriginalData) ) {
      cerr << "Unable to apply VarTransformer to training data." << endl;
      return 2;
    }
    if( t_valid!=0 && !t_valid->transform(t,replaceOriginalData) ) {
      cerr << "Unable to apply VarTransformer to validation data." << endl;
      return 2;
    }
    cout << "Variable transformation from file "
         << transformerFile.c_str() << " has been applied to "
         << "training and validation data." << endl;
    garbage_train.reset(filter.release());
    garbage_valid.reset(valFilter.release());
    filter.reset(t_train);
    valFilter.reset(t_valid);
  }

  // make per-event loss
  auto_ptr<SprAverageLoss> loss;
  switch( iLoss )
    {
    case 1 :
      loss.reset(new SprAverageLoss(&SprLoss::correct_id));
      cout << "Per-event loss set to "
           << "Misidentified fraction of events " << endl;
      break;
    case 2 :
      loss.reset(new SprAverageLoss(&SprLoss::quadratic));
      cout << "Per-event loss set to "
           << "Quadratic loss (y-f(x))^2 " << endl;
      break;
    default :
      loss.reset(new SprAverageLoss(&SprLoss::correct_id));
      cout << "No per-event loss is chosen. " 
	   << "Will use misid fraction of events by default." << endl;
    }

  // open file with classifier configs
  ifstream file(configFile.c_str());
  if( !file ) {
    cerr << "Unable to open file " << configFile.c_str() << endl;
    return 3;
  }

  // prepare vectors of objects
  vector<SprAbsTwoClassCriterion*> criteria;
  vector<SprAbsClassifier*> destroyC;// classifiers to be deleted
  vector<SprIntegerBootstrap*> bstraps;
  vector<SprCCPair> useC;// classifiers and cuts to be used

  // read classifier params
  unsigned nLine = 0;
  bool discreteTree = false;
  bool mixedNodesTree = false;
  bool fastSort = true;
  bool readOneEntry = true;
  if( !SprClassifierReader::readTrainableConfig(file,nLine,filter.get(),
						discreteTree,mixedNodesTree,
						fastSort,criteria,
						bstraps,destroyC,useC,
						readOneEntry) ) {
    cerr << "Unable to read weak classifier configurations from file " 
	 << configFile.c_str() << endl;
    prepareExit(criteria,destroyC,bstraps);
    return 4;
  }
  cout << "Finished reading " << useC.size() << " classifiers from file "
       << configFile.c_str() << endl;

  // get trainable classifier
  assert( useC.size() == 1 );
  SprAbsClassifier* trainable = useC[0].first;
  
  // build models
  std::vector<SprClassifierEvaluator::SetVarsAndValue> vars_and_loss;
  if( !SprClassifierEvaluator::addNremoveR(filter.get(),
					   valFilter.get(),
					   trainable,0,
					   N,R,loss.get(),nCross,
					   vars_and_loss,
					   integrate,
					   verbose) ) {
    cerr << "Unable to build add N remove R models." << endl;
    prepareExit(criteria,destroyC,bstraps);
    return 5;
  }

  //
  // print out models
  //
  cout << "=========================================================" << endl;
  cout << "N_of_variables      Loss" << endl;
  for( int i=0;i<vars_and_loss.size();i++ ) {
    cout << setw(5) << i+1 << "                " 
	 << setw(10) << vars_and_loss[i].second.first 
	 << " +- " << setw(10) << vars_and_loss[i].second.second << endl;
  }
  cout << "Variables:" << endl;
  for( int i=0;i<vars_and_loss.size();i++ ) {
    cout << setw(5) << i+1 << "        ";
    for( set<string>::const_iterator 
	   j=vars_and_loss[i].first.begin();j!=vars_and_loss[i].first.end();
	 j++ ) {
      cout << *j << "   ";
    }
    cout << endl;
  }
  cout << "=========================================================" << endl;

  // exit
  prepareExit(criteria,destroyC,bstraps);
  return 0;
}
