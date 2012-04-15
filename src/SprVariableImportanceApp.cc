//$Id: SprVariableImportanceApp.cc,v 1.8 2008-04-02 23:36:45 narsky Exp $
//
// An executable to estimate the relative importance of variables.
// See notes in README (Variable Selection).
//

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprEmptyFilter.hh"
#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprAbsReader.hh"
#include "StatPatternRecognition/SprRWFactory.hh"
#include "StatPatternRecognition/SprStringParser.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprClassifierReader.hh"
#include "StatPatternRecognition/SprMultiClassReader.hh"
#include "StatPatternRecognition/SprCoordinateMapper.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprClassifierEvaluator.hh"
#include "StatPatternRecognition/SprAbsVarTransformer.hh"
#include "StatPatternRecognition/SprVarTransformerReader.hh"
#include "StatPatternRecognition/SprTransformerFilter.hh"
#include "StatPatternRecognition/SprAverageLoss.hh"
#include "StatPatternRecognition/SprLoss.hh"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <set>
#include <vector>
#include <memory>
#include <string>
#include <cassert>
#include <algorithm>
#include <utility>
#include <iomanip>

using namespace std;


void help(const char* prog) 
{
  cout << "Usage:  " << prog << " classifier_config_file"
       << " input_data_file" << endl;
  cout << "\t Options: " << endl;
  cout << "\t-h --- help                                        " << endl;
  cout << "\t-y list of input classes (see SprAbsFilter.hh)     " << endl;
  cout << "\t-Q apply variable transformation saved in file     " << endl;
  cout << "\t-a input ascii file mode (see SprSimpleReader.hh)  " << endl;
  cout << "\t-k keep the specified fraction in input data       " << endl;
  cout << "\t-K keep (1-this_fraction) in input data            " << endl;
  cout << "\t\t For consistency with other executables,         " << endl;
  cout << "\t\t this option will use \"test\" data to estimate "
       << "variable importance." << endl;
  cout << "\t-g per-event loss to be monitored                  " << endl;
  cout << "\t\t 1 - misidentified fraction (default)            " << endl;
  cout << "\t\t 2 - quadratic loss (y-f(x))^2                   " << endl;
  cout << "\t-C use multiclass learner                          " << endl;
  cout << "\t-n number of class permutations per variable       " << endl;
  cout << "\t\t default=0 - variable importance not estimated   " << endl;
  cout << "\t-S subset of variables used to compute interactions" << endl;
  cout << "\t\t See README for a full description of the syntax." << endl;
  cout << "\t-s sort variables by interaction                   " << endl;
  cout << "\t-N number of points used for data integration      " << endl;
  cout << "\t\t -N is only used for computation of interactions "
       << "between variables. "                                   << endl;
  cout << "\t\t The greater the N, the more accurate the estimate." << endl;
  cout << "\t\t The smaller the N, the faster you get results."   << endl;
  cout << "\t\t N cannot exceed data size."                       << endl;
  cout << "\t\t By default all available points in data are used."<< endl;
  cout << "\t-v verbose level (0=silent default,1,2)            " << endl;
  cout << "\t-w scale all signal weights by this factor         " << endl;
  cout << "\t-V include only these input variables              " << endl;
  cout << "\t-z exclude input variables from the list           " << endl;
  cout << "\t-M map variable lists from trained classifiers onto" << endl;
  cout << "\t\t variables available in input data."               << endl;
  cout << "\t\t Variables must be listed in quotes and separated by commas." 
       << endl;
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
  bool mapTrainedVars = false;
  int nPerm = 0;
  bool useMCLearner = false;
  string transformerFile;
  unsigned nPoints = 0;
  bool computeInteraction = false;
  bool sortByInteraction = false;
  string varList;
  bool split = false;
  double splitFactor = 0;
  bool useTrainingData = false;
  int iLoss = 0;

  // decode command line
  int c;
  extern char* optarg;
  extern int optind;
  while( (c = getopt(argc,argv,"hy:Q:a:k:K:g:Cn:S:sN:v:w:V:z:M")) != EOF ) {
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
      case 'k' :
	split = true;
	splitFactor = (optarg==0 ? 0 : atof(optarg));
	useTrainingData = true;
	break;
      case 'K' :
	split = true;
	splitFactor = (optarg==0 ? 0 : atof(optarg));
	useTrainingData = false;
	break;
      case 'g' :
        iLoss = (optarg==0 ? 1 : atoi(optarg));
        break;
      case 'C' :
	useMCLearner = true;
	break;
      case 'n' :
	nPerm = (optarg==0 ? 1 : atoi(optarg));
	break;
      case 'S' :
	computeInteraction = true;
	varList = (optarg==0 ? "" : optarg);
	break;
      case 's' :
	sortByInteraction = true;
	break;
      case 'N' :
	nPoints = (optarg==0 ? 1 : atoi(optarg));
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
      case 'M' :
	mapTrainedVars = true;
	break;
      }
  }

  // sanity check
  if( nPerm<=0 && !computeInteraction && !sortByInteraction ) {
    cerr << "Neither variable importance nor interaction are to be " 
	 << "estimated. Exiting."<< endl;
    return 1;
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

  // read input data from file
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

  // split data
  auto_ptr<SprAbsFilter> garbage_split;
  if( split ) {
    cout << "Splitting input data with factor " << splitFactor << endl;
    vector<double> weights;
    bool splitRandomize = false;
    SprData* splitted = filter->split(splitFactor,weights,splitRandomize);
    if( splitted == 0 ) {
      cerr << "Unable to split input data." << endl;
      return 2;
    }
    if( !useTrainingData ) {
      garbage_split.reset(filter.release());
      bool ownData = true;
      filter.reset(new SprEmptyFilter(splitted,weights,ownData));
    }
    cout << "Input data re-filtered:" << endl;
    for( int i=0;i<inputClasses.size();i++ ) {
      cout << "Points in class " << inputClasses[i] << ":   " 
	   << filter->ptsInClass(inputClasses[i]) << endl;
    }
  }

  // apply transformation of variables to training and test data
  auto_ptr<SprAbsFilter> garbage_trans;
  if( !transformerFile.empty() ) {
    const SprAbsVarTransformer* t 
      = SprVarTransformerReader::read(transformerFile.c_str());
    if( t == 0 ) {
      cerr << "Unable to read VarTransformer from file "
           << transformerFile.c_str() << endl;
      return 2;
    }
    SprTransformerFilter* tf = new SprTransformerFilter(filter.get());
    bool replaceOriginalData = true;
    if( !tf->transform(t,replaceOriginalData) ) {
      cerr << "Unable to apply VarTransformer to training data." << endl;
      return 2;
    }
    cout << "Variable transformation from file "
         << transformerFile.c_str() << " has been applied to "
         << "training data." << endl;
    garbage_trans.reset(filter.release());
    filter.reset(tf);
    filter->vars(vars);
  }

  // read classifier configuration
  auto_ptr<SprAbsTrainedClassifier> trained;
  auto_ptr<SprAbsTrainedMultiClassLearner> mcTrained;
  if( useMCLearner ) {
    mcTrained.reset(SprMultiClassReader::readTrained(configFile.c_str()));
    if( mcTrained.get() == 0 ) {
      cerr << "Failed to read saved multi class learner from file "
           << configFile.c_str() << endl;
      return 3;
    }
    cout << "Read classifier " << mcTrained->name().c_str()
	 << " with dimensionality " << mcTrained->dim() << endl;
  }
  else {
    trained.reset(SprClassifierReader::readTrained(configFile.c_str(),
						   verbose));
    if( trained.get() == 0 ) {
      cerr << "Unable to read classifier configuration from file "
	   << configFile.c_str() << endl;
      return 3;
    }
    cout << "Read classifier " << trained->name().c_str()
	 << " with dimensionality " << trained->dim() << endl;
  }

  // get a list of trained variables
  vector<string> trainedVars;
  if( trained.get() != 0 )
    trained->vars(trainedVars);
  else
    mcTrained->vars(trainedVars);
  if( verbose > 0 ) {
    cout << "Variables:      " << endl;
    for( int j=0;j<trainedVars.size();j++ ) 
      cout << trainedVars[j].c_str() << " ";
    cout << endl;
  }

  // map trained-classifier variables onto data variables
  auto_ptr<SprCoordinateMapper> mapper;
  if( mapTrainedVars || 
      (trained.get()!=0 && trained->name()=="Combiner") ) {
    mapper.reset(SprCoordinateMapper::createMapper(trainedVars,vars));
    if( mapper.get() == 0 ) {
      cerr << "Unable to map trained classifier vars onto data vars." << endl;
      return 4;
    }
  }

  // set loss
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

  // call evaluator
  vector<SprClassifierEvaluator::NameAndValue> lossIncrease, interaction;
  vector<SprClassifierEvaluator::ListVarsAndValue> sortedInteraction;
  if( nPerm>0 && 
      !SprClassifierEvaluator::variableImportance(filter.get(),
						  trained.get(),
						  mcTrained.get(),
						  loss.get(),
						  mapper.get(),
						  nPerm,
						  lossIncrease) ) {
    cerr << "Unable to estimate variable importance." << endl;
    return 5;
  }
  if( computeInteraction &&
      !SprClassifierEvaluator::variableInteraction(filter.get(),
						   trained.get(),
						   mcTrained.get(),
						   mapper.get(),
						   varList.c_str(),
						   nPoints,
						   interaction,
						   verbose) ) {
    cerr << "Unable to estimate variable interactions." << endl;
    return 6;
  }
  if( sortByInteraction &&
      !SprClassifierEvaluator::sortByInteraction(filter.get(),
						 trained.get(),
						 mcTrained.get(),
						 mapper.get(),
						 nPoints,
						 sortedInteraction,
						 verbose) ) {
    cerr << "Unable to sort variables by interaction." << endl;
    return 7;
  }

  //
  // process computed loss
  //
  cout << "===================================================================================" << endl;
  if( !lossIncrease.empty() ) {
    char t [200];
    sprintf(t,"%35s        %15s","Variable","Change in loss");
    cout << t << endl;
    for( int d=0;d<lossIncrease.size();d++ ) {
      char s [200];
      sprintf(s,"%35s      %15.10f +- %15.10f",lossIncrease[d].first.c_str(),
	      lossIncrease[d].second.first,lossIncrease[d].second.second);
      cout << s << endl;
    }
  }
  cout << "===================================================================================" << endl;
  if( !interaction.empty() ) {
    char t [200];
    sprintf(t,"%35s        %15s","Subset","Interaction");
    cout << t << endl;
    for( int d=0;d<interaction.size();d++ ) {
      char s [200];
      sprintf(s,"%35s      %15.10f +- %15.10f",interaction[d].first.c_str(),
	      interaction[d].second.first,interaction[d].second.second);
      cout << s << endl;
    }
  }
  cout << "===================================================================================" << endl;
  if( !sortedInteraction.empty() ) {
    cout << "N_of_variables      Interaction" << endl;
    for( int i=0;i<sortedInteraction.size();i++ ) {
      cout << setw(5) << sortedInteraction[i].first.size() 
	   << "                " 
	   << setw(10) << sortedInteraction[i].second.first 
	   << " +- " << setw(10) << sortedInteraction[i].second.second << endl;
    }
    cout << "Variables:" << endl;
    for( int i=0;i<sortedInteraction.size();i++ ) {
      cout << setw(5) << i+1 << "        ";
      for( int j=0;j<sortedInteraction[i].first.size();j++ ) {
	cout << (sortedInteraction[i].first)[j].c_str() << "   ";
      }
      cout << endl;
    }
  }
  cout << "===================================================================================" << endl;

  // exit
  return 0;
}
