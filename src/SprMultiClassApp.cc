//$Id: SprMultiClassApp.cc,v 1.15 2008-05-08 19:57:43 narsky Exp $
/*
  Note: "-y" option has a different meaning for this executable than
  for other executables in the package. Instead of specifying what
  classes should be treated as background and what classes should be
  treated as signal, the "-y" option simply selects input classes for
  inclusion in the multi-class algorithm. Therefore, entering groups
  of classes separated by semicolons or specifying "." as an input
  class would make no sense. This executable expects a list of classes
  separated by commas.
*/

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprAbsClassifier.hh"
#include "StatPatternRecognition/SprAdaBoost.hh"
#include "StatPatternRecognition/SprBagger.hh"
#include "StatPatternRecognition/SprStdBackprop.hh"
#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprData.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprAbsMultiClassLearner.hh"
#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprMultiClassLearner.hh"
#include "StatPatternRecognition/SprTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprBinaryEncoder.hh"
#include "StatPatternRecognition/SprTrainedBinaryEncoder.hh"
#include "StatPatternRecognition/SprEmptyFilter.hh"
#include "StatPatternRecognition/SprAbsReader.hh"
#include "StatPatternRecognition/SprAbsWriter.hh"
#include "StatPatternRecognition/SprDataFeeder.hh"
#include "StatPatternRecognition/SprRWFactory.hh"
#include "StatPatternRecognition/SprAbsTwoClassCriterion.hh"
#include "StatPatternRecognition/SprIntegerBootstrap.hh"
#include "StatPatternRecognition/SprStringParser.hh"
#include "StatPatternRecognition/SprLoss.hh"
#include "StatPatternRecognition/SprTransformation.hh"
#include "StatPatternRecognition/SprAverageLoss.hh"
#include "StatPatternRecognition/SprMultiClassReader.hh"
#include "StatPatternRecognition/SprClassifierReader.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprMultiClassPlotter.hh"
#include "StatPatternRecognition/SprAbsVarTransformer.hh"
#include "StatPatternRecognition/SprVarTransformerReader.hh"
#include "StatPatternRecognition/SprTransformerFilter.hh"
#include "StatPatternRecognition/SprCrossValidator.hh"
#include "StatPatternRecognition/SprClassifierEvaluator.hh"
#include "StatPatternRecognition/SprCoordinateMapper.hh"
#include "StatPatternRecognition/SprReplaceMissing.hh"

#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <memory>
#include <iomanip>

using namespace std;


void help(const char* prog) 
{
  cout << "Usage:  " << prog << " training_data_file" << endl;
  cout << "\t Options: " << endl;
  cout << "\t-h --- help                                             " << endl;
  cout << "\t-o output Tuple file                                    " << endl;
  cout << "\t-a input ascii file mode (see SprSimpleReader.hh)       " << endl;
  cout << "\t-A save output data in ascii instead of Root            " << endl;
  cout << "\t-y list of input classes                                " << endl;
  cout << "\t\t Classes must be listed in quotes and separated by commas." 
       << endl;
  cout << "\t-w weights to which input classes are normalized (def=unchanged)" 
       << endl;
  cout << "\t-q do not apply class weights to test data " 
       << "(weights are always applied to cross-validated data)"       << endl;
  cout << "\t-Q apply variable transformation saved in file          " << endl;
  cout << "\t-C multiclass algorithm                                 " << endl;
  cout << "\t\t 1 - Allwein-Schapire-Singer (default)                " << endl;
  cout << "\t\t 2 - binary encoder                                   " << endl;
  cout << "\t-W relative weights for classifiers, one weight per column"
       << " of the indicator matrix."                                  << endl;
  cout << "\t\t Weights must be listed in quotes and separated by commas." 
       << endl;
  cout << "\t-e Multi class mode                                     " << endl;
  cout << "\t\t 1 - OneVsAll (default)                               " << endl;
  cout << "\t\t 2 - OneVsOne                                         " << endl;
  cout << "\t\t 3 - user-defined (must use -i option)                " << endl;
  cout << "\t-i input file with user-defined indicator matrix        " << endl;
  cout << "\t-j do not use zero elements of indicator matrix for loss"
       << " computation" << endl;
  cout << "\t-c file with trainable classifier configurations        " << endl;
  cout << "\t-g per-event loss to be used for each row of the " 
       << "indicator matrix"                                           << endl;
  cout << "\t\t 1 - quadratic loss (y-f(x))^2                        " << endl;
  cout << "\t\t 2 - exponential loss exp(-y*f(x))                    " << endl;
  cout << "\t-G per-event loss to be dispalyed as classification error"<< endl;
  cout << "\t\t 1 - misidentified fraction (default)                 " << endl;
  cout << "\t\t 2 - quadratic loss (y-f(x))^2                        " << endl;
  cout << "\t-L replace data values below this cutoff with "
       << "class-specific medians"                                     << endl;
  cout << "\t-v verbose level (0=silent default,1,2)                 " << endl;
  cout << "\t-f store trained multi class learner to file            " << endl;
  cout << "\t-r read multi class learner configuration stored in file" << endl;
  cout << "\t-K keep this fraction in training set and          " << endl;
  cout << "\t\t put the rest into validation set                " << endl;
  cout << "\t-D randomize training set split-up                 " << endl;
  cout << "\t-t read validation/test data from a file           " << endl;
  cout << "\t\t (must be in same format as input data!!!        " << endl;
  cout << "\t-d frequency of print-outs for validation data     " << endl;
  cout << "\t-x use that many pieces for cross-validation       " << endl;
  cout << "\t-N add N variables at a time (default=1) "
       << "using add N remove R algorithm" << endl;
  cout << "\t-R remove R variables at a time (default=0) "
       << "using add N remove R algorithm" << endl;
  cout << "\t-V include only these input variables              " << endl;
  cout << "\t-z exclude input variables from the list           " << endl;
  cout << "\t-Z exclude input variables from the list, "
       << "but put them in the output file " << endl;
  cout << "\t\t Variables must be listed in quotes and separated by commas." 
       << endl;
  cout << "\t-M map variable lists from trained classifiers onto "
       << "variables available in input data."                    << endl; 
}


void prepareExit(vector<SprAbsTwoClassCriterion*>& criteria,
		 vector<SprAbsClassifier*>& classifiers,
		 vector<SprIntegerBootstrap*>& bstraps) 
{
  for( int i=0;i<criteria.size();i++ ) delete criteria[i];
  for( int i=0;i<classifiers.size();i++ ) delete classifiers[i];
  for( int i=0;i<bstraps.size();i++ ) delete bstraps[i];
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
  string resumeFile;
  string configFile;
  string valFile;
  bool setLowCutoff = false;
  double lowCutoff = 0;
  string includeList, excludeList;
  string inputClassesString;
  string inputClassWeightString, inputClassifierWeightString;
  int iOptLoss(1), iMonLoss(1);
  int iMode = 1;
  string indicatorFile;
  string stringVarsDoNotFeed;
  bool split = false;
  double splitFactor = 0;
  bool splitRandomize = false;
  string transformerFile;
  unsigned nCross = 0;
  bool addNRemoveR = false;
  unsigned N(0), R(0);
  bool excludeZeroClasses = false;
  bool applyWeightsToValData = true;
  int mcAlgorithm = 1;
  unsigned valPrint = 0;
  bool mapTrainedVars = false;

  // decode command line
  int c;
  extern char* optarg;
  //  extern int optind;
  while( (c = getopt(argc,argv,"ho:a:Ay:w:qC:W:Q:e:i:jc:g:G:L:bv:f:r:K:Dt:d:x:N:R:V:z:Z:M")) != EOF ) {
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
      case 'w' :
	inputClassWeightString = optarg;
	break;
      case 'q' :
	applyWeightsToValData = false;
	break;
      case 'C' :
	mcAlgorithm = (optarg==0 ? 1 : atoi(optarg));
	break;
      case 'W' :
	inputClassifierWeightString = optarg;
	break;
      case 'Q' :
        transformerFile = optarg;
        break;
      case 'e' :
        iMode = (optarg==0 ? 1 : atoi(optarg));
        break;
      case 'i' :
	indicatorFile = optarg;
	break;
      case 'j' :
	excludeZeroClasses = true;
	break;
      case 'c' :
	configFile = optarg;
	break;
      case 'g' :
        iOptLoss = (optarg==0 ? 1 : atoi(optarg));
        break;
      case 'G' :
        iMonLoss = (optarg==0 ? 1 : atoi(optarg));
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
      case 'f' :
	outFile = optarg;
	break;
      case 'r' :
	resumeFile = optarg;
	break;
      case 'K' :
	split = true;
	splitFactor = (optarg==0 ? 0 : atof(optarg));
	break;
      case 'D' :
	splitRandomize = true;
	break;
      case 't' :
	valFile = optarg;
	break;
      case 'd' :
        valPrint = (optarg==0 ? 0 : atoi(optarg));
        break;
      case 'x' :
	nCross = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'N' :
        N = (optarg==0 ? 1 : atoi(optarg));
	addNRemoveR = true;
        break;
      case 'R' :
        R = (optarg==0 ? 0 : atoi(optarg));
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
      case 'M' :
	mapTrainedVars = true;
	break;
      }
  }

  // sanity check
  if( configFile.empty() && resumeFile.empty() ) {
    cerr << "No classifier configuration file specified." << endl;
    return 1;
  }
  if( !configFile.empty() && !resumeFile.empty() ) {
    cerr << "Cannot train and use saved configuration at the same time." 
	 << endl;
    return 1;
  }
  if( nCross>0 && !valFile.empty() ) {
    cerr << "Cannot validate and cross-validate at the same time." << endl;
    return 1;
  }

  // Must have 2 arguments after all options.
  string trFile = argv[argc-1];
  if( trFile.empty() ) {
    cerr << "No training file is specified." << endl;
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

  // decode input classes
  if( inputClassesString.empty() ) {
    cerr << "No input classes specified." << endl;
    return 2;
  }
  vector<vector<int> > inputIntClasses;
  SprStringParser::parseToInts(inputClassesString.c_str(),inputIntClasses);
  if( inputIntClasses.empty() || inputIntClasses[0].size()<2 ) {
    cerr << "Found less than 2 classes in the input class string." << endl;
    return 2;
  }
  vector<SprClass> inputClasses(inputIntClasses[0].size());
  for( int i=0;i<inputIntClasses[0].size();i++ )
    inputClasses[i] = inputIntClasses[0][i];
  vector<int> classes = inputIntClasses[0];

  // filter training data by class
  filter->chooseClasses(inputClasses);
  if( !filter->filter() ) {
    cerr << "Unable to filter training data by class." << endl;
    return 2;
  }
  cout << "Training data filtered by class." << endl;
  for( int i=0;i<inputClasses.size();i++ ) {
    unsigned npts = filter->ptsInClass(inputClasses[i]);
    if( npts == 0 ) {
      cerr << "Error!!! No points in class " << inputClasses[i] << endl;
      return 2;
    }
    cout << "Points in class " << inputClasses[i] << ":   " << npts << endl;
  }

  // read validation data from file
  auto_ptr<SprAbsFilter> valFilter;
  if( split && !valFile.empty() ) {
    cerr << "Unable to split training data and use validation data " 
	 << "from a separate file." << endl;
    return 2;
  }
  if( split ) {
    cout << "Splitting training data with factor " << splitFactor << endl;
    if( splitRandomize )
      cout << "Will use randomized splitting." << endl;
    vector<double> weights;
    SprData* splitted = filter->split(splitFactor,weights,splitRandomize);
    if( splitted == 0 ) {
      cerr << "Unable to split training data." << endl;
      return 2;
    }
    bool ownData = true;
    valFilter.reset(new SprEmptyFilter(splitted,weights,ownData));
    cout << "Training data re-filtered:" << endl;
    for( int i=0;i<inputClasses.size();i++ ) {
      cout << "Points in class " << inputClasses[i] << ":   " 
	   << filter->ptsInClass(inputClasses[i]) << endl;
    }
  }
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
  }

  // filter validation data by class
  if( valFilter.get() != 0 ) {
    valFilter->chooseClasses(inputClasses);
    if( !valFilter->filter() ) {
      cerr << "Unable to filter validation data by class." << endl;
      return 2;
    }
    cout << "Validation data filtered by class." << endl;
    for( int i=0;i<inputClasses.size();i++ ) {
     unsigned npts = valFilter->ptsInClass(inputClasses[i]);
     if( npts == 0 )
       cerr << "Warning!!! No points in class " << inputClasses[i] << endl;
     cout << "Points in class " << inputClasses[i] << ":   " << npts << endl;
    }
  }

  // decode input class weights
  vector<vector<double> > inputClassWeights;
  SprStringParser::parseToDoubles(inputClassWeightString.c_str(),
				  inputClassWeights);
  if( !inputClassWeights.empty() ) {
    const vector<double>& weights = inputClassWeights[0];
    if( weights.size() != inputClasses.size() ) {
      cerr << "Size of vector of weights not equal to "
	   << "size of vector of classes: " 
	   << weights.size() << " " << inputClasses.size() << endl;
      return 2;
    }
    vector<SprClass> tempClass(1);
    for( int i=0;i<inputClasses.size();i++ ) {
      tempClass[0] = inputClasses[i];
      if( !filter->normalizeWeights(tempClass,weights[i]) 
	  || (valFilter.get()!=0 && applyWeightsToValData
	      && !valFilter->normalizeWeights(tempClass,weights[i])) ) {
	cerr << "Unable to normalize weights for class " 
	     << tempClass[0] << endl;
	return 2;
      }
    }
    cout << "Training data reweighted by class." << endl;
    for( int i=0;i<inputClasses.size();i++ ) {
      double w = filter->weightInClass(inputClasses[i]);
      cout << "Weight in class " << inputClasses[i] << ":   " << w << endl;
    }
  }

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

  // init performance measure
  SprClassificationTable lossTable;
  map<int,double> weightInClass;
  SprValueWithError realLoss;

  // get monitoring loss
  auto_ptr<SprAverageLoss> monLoss;
  switch( iMonLoss )
    {
    case 1 :
      monLoss.reset(new SprAverageLoss(&SprLoss::correct_id));
      cout << "Per-event monitoring loss set to "
	   << "Misid rate int(y==f(x)) " << endl;
      break;
    case 2 :
      monLoss.reset(new SprAverageLoss(&SprLoss::quadratic));
      cout << "Per-event monitoring loss set to "
	   << "Quadratic loss (y-f(x))^2 " << endl;
      break;
    default :
      monLoss.reset(new SprAverageLoss(&SprLoss::correct_id));
      cout << "Per-event monitoring loss set to "
	   << "Misid rate int(y==f(x)) " << endl;
    }

  // prepare trained classifier holder
  auto_ptr<SprAbsTrainedMultiClassLearner> trainedMulti;

  // prepare vectors of objects
  vector<SprAbsTwoClassCriterion*> criteria;
  vector<SprAbsClassifier*> destroyC;// classifiers to be deleted
  vector<SprIntegerBootstrap*> bstraps;
  vector<SprCCPair> useC;// classifiers and cuts to be used

  // open file with classifier configs
  if( !configFile.empty() ) {
    ifstream file(configFile.c_str());
    if( !file ) {
      cerr << "Unable to open file " << configFile.c_str() << endl;
      return 3;
    }
    
    // read classifier params
    unsigned nLine = 0;
    bool discreteTree = false;
    bool mixedNodesTree = false;
    bool fastSort = false;
    bool readOneEntry = true;
    if( !SprClassifierReader::readTrainableConfig(file,nLine,filter.get(),
						  discreteTree,mixedNodesTree,
						  fastSort,criteria,
						  bstraps,destroyC,useC,
						  readOneEntry) ) {
      cerr << "Unable to read classifier configurations from file " 
	   << configFile.c_str() << endl;
      prepareExit(criteria,destroyC,bstraps);
      return 3;
    }
    cout << "Finished reading " << useC.size() << " classifiers from file "
	 << configFile.c_str() << endl;
    assert( useC.size() == 1 );
    SprAbsClassifier* trainable = useC[0].first;

    // find the multi class mode
    SprMultiClassLearner::MultiClassMode multiClassMode 
      = SprMultiClassLearner::OneVsAll;
    switch( iMode )
      {
      case 1 :
        multiClassMode = SprMultiClassLearner::OneVsAll;
        cout << "Multi class learning mode set to OneVsAll." << endl;
        break;
      case 2 :
        multiClassMode = SprMultiClassLearner::OneVsOne;
      	cout << "Multi class learning mode set to OneVsOne." << endl;
  	break;
      case 3:
	if( indicatorFile.empty() ) {
	  cerr << "No indicator matrix specified." << endl;
	  prepareExit(criteria,destroyC,bstraps);
	  return 4;
	}
	multiClassMode = SprMultiClassLearner::User;
	cout << "Multi class learning mode set to User." << endl;
	break;
      default :
        cerr << "No multi class learning mode chosen." << endl;
        prepareExit(criteria,destroyC,bstraps);
        return 4;
      }
    if( excludeZeroClasses ) {
      cout << "Zero elements of the indicator matrix will be ignored "
	   << "for loss computation." << endl;
    }

    // get indicator matrix
    SprMatrix indicator;
    if( multiClassMode==SprMultiClassLearner::User 
	&& !indicatorFile.empty() ) {
      if( !SprMultiClassReader::readIndicatorMatrix(indicatorFile.c_str(),
						    indicator) ) {
	cerr << "Unable to read indicator matrix from file " 
	     << indicatorFile.c_str() << endl;
        prepareExit(criteria,destroyC,bstraps);
	return 4;
      }
    }

    // make a multi class learner
    auto_ptr<SprAbsMultiClassLearner> multi;
    switch( mcAlgorithm )
      {
      case 1 :
	multi.reset(new SprMultiClassLearner(filter.get(),
					     trainable,
					     inputIntClasses[0],
					     indicator,
					     multiClassMode));
	cout << "Multiclass algorithm set to Allwein-Schapire-Singer." << endl;
	break;
      case 2 :
	multi.reset(new SprBinaryEncoder(filter.get(),
					 trainable,
					 inputIntClasses[0]));
	cout << "Multiclass algorithm set to BinaryEncoder." << endl;
	break;
      default :
	cerr << "Unknown multiclass algorithm. Exiting." << endl;
	return 4;
      }
					 
    // ASS operations
    if( multi->name() == "MultiClassLearner" ) {
      SprMultiClassLearner* mcl 
	= static_cast<SprMultiClassLearner*>(multi.get());

      // get the matrix
      mcl->indicatorMatrix(indicator);

      // exclude zero indicator-matrix elements from loss computation
      if( excludeZeroClasses )
	mcl->excludeZeroClasses();

      // decode classifier weights
      vector<vector<double> > inputClassifierWeights;
      SprStringParser::parseToDoubles(inputClassifierWeightString.c_str(),
				      inputClassifierWeights);
      if( !inputClassifierWeights.empty() ) {
	if( inputClassifierWeights[0].size() != indicator.num_col() ) {
	  cerr << "The number of classifier weights is not equal "
	       << "to the number of classifiers: " 
	       <<  inputClassifierWeights[0].size() << " " 
	       << indicator.num_col() << endl;
	  prepareExit(criteria,destroyC,bstraps);
	  return 4;
	}
	mcl->setClassifierWeights(inputClassifierWeights[0]);
      }
    }

    // set print-out frequency for validation data
    auto_ptr<SprAbsFilter> convertedValData;
    if( multi->name()=="BinaryEncoder" && valPrint!=0 ) {
      SprBinaryEncoder converter(valFilter.get(),
				 trainable,
				 inputIntClasses[0]);
      convertedValData.reset(converter.convertData());
      if( convertedValData.get() == 0 ) {
	cerr << "Unable to convert validation data for BinaryEncoder." << endl;
	return 4;
      }
      if(      trainable->name() == "AdaBoost" ) {
	if( !static_cast<SprAdaBoost*>(trainable)
	    ->setValidation(convertedValData.get(),valPrint) ) {
	  cerr << "Unable to set validation loss for BinaryEncoder." << endl;
	  return 4;
	}
      }
      else if( trainable->name()=="Bagger" || trainable->name()=="ArcX4" ) {
	if( !static_cast<SprBagger*>(trainable)
	    ->setValidation(convertedValData.get(),valPrint,0,0) ) {
	  cerr << "Unable to set validation loss for BinaryEncoder." << endl;
	  return 4;
	}
      }
      else if( trainable->name() == "StdBackprop" ) {
	if( !static_cast<SprStdBackprop*>(trainable)
	    ->setValidation(convertedValData.get(),valPrint) ) {
	  cerr << "Unable to set validation loss for BinaryEncoder." << endl;
	  return 4;
	}
      }
    }

    // set missing values for each class individually
    if( setLowCutoff ) {
      bool classBlind = false;
      SprReplaceMissing missingFilter(SprUtils::lowerBound(lowCutoff),
				      classBlind);
      if( !missingFilter.train(filter.get(),verbose) ) {
	cerr << "Cannot compute replacement for class-specific missing values."
	     << endl;
	return 4;
      }
      vector<SprReplaceMissing::ClassAndDefaultValues> dataMissing;
      missingFilter.replacement(dataMissing);
      vector<SprAbsMultiClassLearner::MissingType> multiMissing;
      for( int i=0;i<dataMissing.size();i++ ) {
	vector<int> tempClass;
	if( dataMissing[i].first.value(tempClass) || tempClass.size()!=1 ) {
	  cerr << "Unable to set classes for missing value computation for " 
	       << "MultiClass learner." << endl;
	  prepareExit(criteria,destroyC,bstraps);
	  return 4;
	}
	multiMissing.push_back(SprAbsMultiClassLearner::MissingType(
					       tempClass[0],
					       dataMissing[i].second));
      }
      if( !multi->setDefaultMissing(SprUtils::lowerBound(lowCutoff),
				    multiMissing) ) {
	cerr << "Unable to set class-specific missing values for the "
	     << "MultiClass learner." << endl;
	prepareExit(criteria,destroyC,bstraps);
	return 4;
      }
    }

    // cross-validate
    if( nCross>0 && !addNRemoveR ) {
      SprCrossValidator cv(filter.get(),nCross);
      if( !cv.validate(multi.get(),monLoss.get(),
		       realLoss,lossTable,weightInClass,verbose) ) {
	cerr << "Cannot cross-validate." << endl;
        prepareExit(criteria,destroyC,bstraps);
	return 5;
      }
    }

    // build add N remove R models
    if( addNRemoveR ) {
      std::vector<SprClassifierEvaluator::SetVarsAndValue> vars_and_loss;
      bool integrate = true;
      if( !SprClassifierEvaluator::addNremoveR(filter.get(),
					       valFilter.get(),
					       0,multi.get(),
					       N,R,
					       monLoss.get(),nCross,
					       vars_and_loss,
					       integrate,verbose) ) {
	cerr << "Unable to build add N remove R models." << endl;
	prepareExit(criteria,destroyC,bstraps);
	return 5;
      }

      // print out models
      cout << "=========================================================" 
	   << endl;
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
	       j=vars_and_loss[i].first.begin();
	     j!=vars_and_loss[i].first.end();j++ ) {
	  cout << *j << "   ";
	}
	cout << endl;
      }
      cout << "=========================================================" 
	   << endl;
    }// end if addNRemoveR

    // train
    if( resumeFile.empty() && nCross==0 && !addNRemoveR ) {
      if( !multi->train(verbose) ) {
        cerr << "Unable to train Multi class learner." << endl;
        prepareExit(criteria,destroyC,bstraps);
        return 5;
      }
      else {
        trainedMulti.reset(multi->makeTrained());
        cout << "Multi class learner finished successfully." << endl;
      }

      // save trained multi class learner
      if( !outFile.empty() ) {
	if( !multi->store(outFile.c_str()) ) {
	  cerr << "Cannot store multi class learner in file " 
	       << outFile.c_str() << endl;
	  prepareExit(criteria,destroyC,bstraps);
	  return 6;
	}
      }
    }
  }// end train

  // read saved learner from file
  if( !resumeFile.empty() ) {
    trainedMulti.reset(SprMultiClassReader::readTrained(resumeFile.c_str()));
    if( trainedMulti.get() == 0 ) {
      cerr << "Failed to read saved multi class learner from file " 
	   << resumeFile.c_str() << endl;
      prepareExit(criteria,destroyC,bstraps);
      return 7;
    }
    cout << "Read saved multi class learner from file " 
	 << resumeFile.c_str() << endl;
    if( trainedMulti->name() == "MultiClassLearner" ) {
      /// cast
      SprTrainedMultiClassLearner* tmcl 
	= static_cast<SprTrainedMultiClassLearner*>(trainedMulti.get());

      // get rid of neutral classes
      if( excludeZeroClasses ) tmcl->excludeZeroClasses();

      // print out matrix
      tmcl->printIndicatorMatrix(cout);
    }
  }

  // by now the trained learner should be filled
  if( trainedMulti.get()==0 && nCross==0 && !addNRemoveR ) {
    cerr << "Trained multi learner has not been set." << endl;
    prepareExit(criteria,destroyC,bstraps);
    return 8;
  }

  //
  // process the trained learner
  //
  if( trainedMulti.get() != 0 ) {
    // set loss
    if( trainedMulti->name()=="MultiClassLearner" ) {
      SprTrainedMultiClassLearner* tmcl 
	= static_cast<SprTrainedMultiClassLearner*>(trainedMulti.get());
      switch( iOptLoss )
	{
	case 1 :
	  tmcl->setLoss(&SprLoss::quadratic,
			&SprTransformation::zeroOneToMinusPlusOne);
	  cout << "Per-event loss for rows of the indicator matrix set to "
	       << "Quadratic loss (y-f(x))^2 " << endl;
	  break;
	case 2 :
	  tmcl->setLoss(&SprLoss::exponential,
			&SprTransformation::logitInverse);
	  cout << "Per-event loss for rows of the indicator matrix set to "
	       << "Exponential loss exp(-y*f(x)) " << endl;
	  break;
	default :
	  cerr << "No per-event for rows of the indicator matrix "
	       << "loss specified." << endl;
	  prepareExit(criteria,destroyC,bstraps);
	  return 9;
	}
    }

    // analyze validation data
    const SprAbsFilter* valData
      = ( resumeFile.empty() ? valFilter.get() : filter.get() );
    if( valData != 0 ) {

      // check dimensionality
      if( valData->dim()!=trainedMulti->dim() && !mapTrainedVars ) {
	cerr << "Dimensionality of data and multiclass learner do not match." 
	     << " Have you forgot -M option?" << endl;
	prepareExit(criteria,destroyC,bstraps);
	return 10;
      }
      
      // map variables?
      vector<string> valVars, tableVars;
      valData->vars(valVars);
      trainedMulti->vars(tableVars);
      auto_ptr<SprCoordinateMapper> tableMapper;
      if( mapTrainedVars ) {
	tableMapper.reset(SprCoordinateMapper::createMapper(tableVars,
							    valVars));
	if( tableMapper.get() == 0 ) {
	  cerr << "Unable to create variable mapper." << endl;
	  prepareExit(criteria,destroyC,bstraps);
	  return 10;
	}
      }

      // compute response
      vector<SprMultiClassPlotter::Response> responses(valData->size());
      for( int i=0;i<valData->size();i++ ) {
	if( ((i+1)%1000) == 0 )
	  cout << "Computing response for validation point " << i+1 << endl;
	
	// get point, class and weight
	const SprPoint* p = (*valData)[i];
	const SprPoint* pMapped
	  = ( tableMapper.get()==0 ? p : tableMapper->output(p) );
	int cls = p->class_;
	double w = valData->w(i);
	
	// compute loss
	map<int,double> output;
	int resp = trainedMulti->response(pMapped,output);
	responses[i] = SprMultiClassPlotter::Response(cls,w,resp,output);

	// clear mapper
	if( tableMapper.get() != 0 ) tableMapper->clear();
      }    

      // get the loss table
      SprMultiClassPlotter plotter(responses);
      realLoss = SprValueWithError(plotter.multiClassTable(classes,
							   monLoss.get(),
							   lossTable,
							   weightInClass),0);
    }// end valFilter.get() != 0
  }// end trainedMulti.get() != 0

  // print out
  if( !lossTable.empty() ) {
    cout << "=====================================" << endl;
    cout << "Overall validation loss = " << realLoss.first << endl;
    cout << "=====================================" << endl;
    cout << "Classification table: Fractions of total class weight" << endl;
    char s[200];
    sprintf(s,"True Class \\ Classification |");
    string temp = "------------------------------";
    cout << s;
    for( int i=0;i<classes.size();i++ ) {
      sprintf(s," %5i      |",classes[i]);
      cout << s;
      temp += "-------------";
    }
    sprintf(s,"   Total weight in class |");
    temp += "-------------------------";
    cout << s << endl;
    cout << temp.c_str() << endl;
    for( map<int,vector<double> >::const_iterator
	   i=lossTable.begin();i!=lossTable.end();i++ ) {
      sprintf(s,"%5i                       |",i->first);
      cout << s;
      for( int j=0;j<i->second.size();j++ ) {
	sprintf(s," %10.4f |",i->second[j]);
	cout << s;
      }
      sprintf(s,"              %10.4f |",weightInClass[i->first]);
      cout << s << endl;
    }
    cout << temp.c_str() << endl;
  }// end of !lossTable.empty()

  // make histogram if requested
  if( tupleFile.empty() ) {
    prepareExit(criteria,destroyC,bstraps);
    return 0;
  }

  // make a writer
  auto_ptr<SprAbsWriter> tuple(SprRWFactory::makeWriter(writeMode,"training"));
  if( !tuple->init(tupleFile.c_str()) ) {
    cerr << "Unable to open output file " << tupleFile.c_str() << endl;
    prepareExit(criteria,destroyC,bstraps);
    return 10;
  }

  // determine if certain variables are to be excluded from usage,
  // but included in the output storage file (-Z option)
  string printVarsDoNotFeed;
  vector<vector<string> > varsDoNotFeed;
  SprStringParser::parseToStrings(stringVarsDoNotFeed.c_str(),varsDoNotFeed);
  vector<unsigned> mapper;
  for( int d=0;d<vars.size();d++ ) {
    if( varsDoNotFeed.empty() ||
        (find(varsDoNotFeed[0].begin(),varsDoNotFeed[0].end(),vars[d])
	 ==varsDoNotFeed[0].end()) ) {
      mapper.push_back(d);
    }
    else {
      printVarsDoNotFeed += ( printVarsDoNotFeed.empty() ? "" : ", " );
      printVarsDoNotFeed += vars[d];
    }
  }
  if( !printVarsDoNotFeed.empty() ) {
    cout << "The following variables are not used in the algorithm, " 
         << "but will be included in the output file: " 
         << printVarsDoNotFeed.c_str() << endl;
  }

  // map variables?
  vector<string> trainedVars;
  trainedMulti->vars(trainedVars);
  SprCoordinateMapper* trainedMapper = 0;
  if( mapTrainedVars ) {
    trainedMapper = SprCoordinateMapper::createMapper(trainedVars,vars);
    if( trainedMapper == 0 ) {
      cerr << "Unable to create variable mapper." << endl;
      prepareExit(criteria,destroyC,bstraps);
      return 11;
    }
  }

  // feed
  SprDataFeeder feeder(filter.get(),tuple.get(),mapper);
  if( !feeder.addMultiClassLearner(trainedMulti.get(),
				   "multi",trainedMapper) ) {
    cerr << "Unable to add classifier to feeder." << endl;
    prepareExit(criteria,destroyC,bstraps);
    return 12;
  }
  if( !feeder.feed(1000) ) {
    cerr << "Cannot feed data into file " << tupleFile.c_str() << endl;
    prepareExit(criteria,destroyC,bstraps);
    return 12;
  }

  // cleanup
  prepareExit(criteria,destroyC,bstraps);

  // exit
  return 0;
}
