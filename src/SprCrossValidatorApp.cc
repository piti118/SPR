//$Id: SprCrossValidatorApp.cc,v 1.4 2008-04-02 23:36:45 narsky Exp $
//
// An executable for cross-validation by an arbitrary classifier.
// See notes in README (Cross-validation).
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
#include "StatPatternRecognition/SprCrossValidator.hh"
#include "StatPatternRecognition/SprAbsVarTransformer.hh"
#include "StatPatternRecognition/SprVarTransformerReader.hh"
#include "StatPatternRecognition/SprTransformerFilter.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprAbsTwoClassCriterion.hh"
#include "StatPatternRecognition/SprTwoClassSignalSignif.hh"
#include "StatPatternRecognition/SprTwoClassIDFraction.hh"
#include "StatPatternRecognition/SprTwoClassTaggerEff.hh"
#include "StatPatternRecognition/SprTwoClassPurity.hh"
#include "StatPatternRecognition/SprTwoClassGiniIndex.hh"
#include "StatPatternRecognition/SprTwoClassCrossEntropy.hh"
#include "StatPatternRecognition/SprTwoClassUniformPriorUL90.hh"
#include "StatPatternRecognition/SprTwoClassBKDiscovery.hh"
#include "StatPatternRecognition/SprTwoClassPunzi.hh"
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
  cout << "\t-C criterion to be monitored on test data (def=1)  " << endl;
  cout << "\t\t 1 = correctly classified fraction               " << endl;
  cout << "\t\t 2 = signal significance s/sqrt(s+b)             " << endl;
  cout << "\t\t 3 = purity s/(s+b)                              " << endl;
  cout << "\t\t 4 = tagger efficiency Q                         " << endl;
  cout << "\t\t 5 = Gini index (default)                        " << endl;
  cout << "\t\t 6 = cross-entropy                               " << endl;
  cout << "\t\t 7 = 90% Bayesian upper limit with uniform prior " << endl;
  cout << "\t\t 8 = discovery potential 2*(sqrt(s+b)-sqrt(b))   " << endl;
  cout << "\t\t 9 = Punzi's sensitivity s/(0.5*nSigma+sqrt(b))  " << endl;
  cout << "\t\t -P background normalization factor for Punzi FOM" << endl;
  cout << "\t-g per-event loss for (cross-)validation           " << endl;
  cout << "\t\t 1 - quadratic loss (y-f(x))^2                   " << endl;
  cout << "\t\t 2 - exponential loss exp(-y*f(x))               " << endl;
  cout << "\t\t 3 - misid fraction                              " << endl;
  cout << "\t-x use that many pieces for cross-validation (def=2)"<< endl;
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
  double bW = 1.;
  string includeList, excludeList;
  string inputClassesString;
  string transformerFile;
  int iLoss = 0;
  int iCrit = 1;
  unsigned nCross = 2;
  bool integrate = false;

  // decode command line
  int c;
  extern char* optarg;
  extern int optind;
  while( (c = getopt(argc,argv,"hy:Q:a:C:P:g:x:iv:w:V:z:")) != EOF ) {
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
      case 'C' :
        iCrit = (optarg==0 ? 1 : atoi(optarg));
        break;
      case 'P' :
        bW = (optarg==0 ? 1 : atof(optarg));
        break;
      case 'g' :
        iLoss = (optarg==0 ? 1 : atoi(optarg));
        break;
      case 'x' :
	nCross = (optarg==0 ? 2 : atoi(optarg));
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

  // sanity check
  if( nCross < 2 ) {
    cerr << "No pieces for cross-validation are specified. Exiting." << endl;
    return false;
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

  // apply transformation of variables to training data
  auto_ptr<SprAbsFilter> garbage_train;
  if( !transformerFile.empty() ) {
    const SprAbsVarTransformer* t 
      = SprVarTransformerReader::read(transformerFile.c_str());
    if( t == 0 ) {
      cerr << "Unable to read VarTransformer from file "
           << transformerFile.c_str() << endl;
      return 2;
    }
    SprTransformerFilter* t_train = new SprTransformerFilter(filter.get());
    bool replaceOriginalData = true;
    if( !t_train->transform(t,replaceOriginalData) ) {
      cerr << "Unable to apply VarTransformer to training data." << endl;
      return 2;
    }
    cout << "Variable transformation from file "
         << transformerFile.c_str() << " has been applied to "
         << "training data." << endl;
    garbage_train.reset(filter.release());
    filter.reset(t_train);
  }

  // make optimization criterion
  auto_ptr<SprAbsTwoClassCriterion> crit;
  switch( iCrit )
    {
    case 1 :
      crit.reset(new SprTwoClassIDFraction);
      cout << "Monitoring criterion set to "
           << "Fraction of correctly classified events " << endl;
      break;
    case 2 :
      crit.reset(new SprTwoClassSignalSignif);
      cout << "Monitoring criterion set to "
           << "Signal significance S/sqrt(S+B) " << endl;
      break;
    case 3 :
      crit.reset(new SprTwoClassPurity);
      cout << "Monitoring criterion set to "
           << "Purity S/(S+B) " << endl;
      break;
    case 4 :
      crit.reset(new SprTwoClassTaggerEff);
      cout << "Monitoring criterion set to "
           << "Tagging efficiency Q = e*(1-2w)^2 " << endl;
      break;
    case 5 :
      crit.reset(new SprTwoClassGiniIndex);
      cout << "Monitoring criterion set to "
           << "Gini index  -1+p^2+q^2 " << endl;
      break;
    case 6 :
      crit.reset(new SprTwoClassCrossEntropy);
      cout << "Monitoring criterion set to "
           << "Cross-entropy p*log(p)+q*log(q) " << endl;
      break;
    case 7 :
      crit.reset(new SprTwoClassUniformPriorUL90);
      cout << "Monitoring criterion set to "
           << "Inverse of 90% Bayesian upper limit with uniform prior" << endl;
      break;
    case 8 :
      crit.reset(new SprTwoClassBKDiscovery);
      cout << "Monitoring criterion set to "
           << "Discovery potential 2*(sqrt(S+B)-sqrt(B))" << endl;
      break;
    case 9 :
      crit.reset(new SprTwoClassPunzi(bW));
      cout << "Monitoring criterion set to "
           << "Punzi's sensitivity S/(0.5*nSigma+sqrt(B))" << endl;
      break;
    default :
      cerr << "Unable to make initialization criterion." << endl;
      return 3;
    }

  // make per-event loss
  auto_ptr<SprAverageLoss> loss;
  switch( iLoss )
    {
    case 1 :
      loss.reset(new SprAverageLoss(&SprLoss::quadratic));
      cout << "Per-event loss set to "
           << "Quadratic loss (y-f(x))^2 " << endl;
      break;
    case 2 :
      loss.reset(new SprAverageLoss(&SprLoss::purity_ratio));
      cout << "Per-event loss set to "
           << "Exponential loss exp(-y*f(x)) " << endl;
      break;
    case 3 :
      loss.reset(new SprAverageLoss(&SprLoss::correct_id));
      cout << "Per-event loss set to "
	   << "Misid rate int(y==f(x)) " << endl;
      break;
    default :
      cout << "No per-event loss is chosen. Will use the default." << endl;
      break;
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
  bool readOneEntry = false;
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
  assert( !useC.empty() );

  // get trainable classifier
  vector<SprAbsClassifier*> classifiers;
  for( int i=0;i<useC.size();i++ )
    classifiers.push_back(useC[i].first);
  
  // cross-validate
  vector<SprValueWithError> cvFom;
  SprCrossValidator cv(filter.get(),nCross);
  if( !cv.validate(crit.get(),loss.get(),classifiers,
		   inputClasses[0],inputClasses[1],
		   cvFom,integrate,verbose) ) {
    cerr << "Unable to cross-validate." << endl;
    prepareExit(criteria,destroyC,bstraps);
    return 5;
  }
  assert( classifiers.size() == cvFom.size() );

  // print out CV FOMs
  cout << "Cross-validated FOMs:" << endl;
  for( int i=0;i<cvFom.size();i++ ) {
    cout << "Classifier " << setw(5) << i+1
	 << setw(20) << classifiers[i]->name().c_str()
	 << "      FOM=" << setw(10) << cvFom[i].first
	 << " +- " << setw(10) << cvFom[i].second << endl;
  }

  // exit
  prepareExit(criteria,destroyC,bstraps);
  return 0;
}
