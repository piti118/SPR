//$Id: SprGEPApp.cc,v 1.5 2008-05-09 21:25:26 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprAdaBoost.hh"
#include "StatPatternRecognition/SprTrainedAdaBoost.hh"
#include "StatPatternRecognition/SprGEP.hh"
#include "StatPatternRecognition/SprTrainedGEP.hh"
#include "StatPatternRecognition/SprEmptyFilter.hh"
#include "StatPatternRecognition/SprAbsReader.hh"
#include "StatPatternRecognition/SprAbsWriter.hh"
#include "StatPatternRecognition/SprDataFeeder.hh"
#include "StatPatternRecognition/SprRWFactory.hh"
#include "StatPatternRecognition/SprStringParser.hh"
#include "StatPatternRecognition/SprIntegerBootstrap.hh"
#include "StatPatternRecognition/SprTwoClassSignalSignif.hh"
#include "StatPatternRecognition/SprTwoClassIDFraction.hh"
#include "StatPatternRecognition/SprTwoClassTaggerEff.hh"
#include "StatPatternRecognition/SprTwoClassPurity.hh"
#include "StatPatternRecognition/SprTwoClassGiniIndex.hh"
#include "StatPatternRecognition/SprTwoClassCrossEntropy.hh"
#include "StatPatternRecognition/SprTwoClassUniformPriorUL90.hh"
#include "StatPatternRecognition/SprTwoClassBKDiscovery.hh"
#include "StatPatternRecognition/SprTwoClassPunzi.hh"
#include "StatPatternRecognition/SprTwoClassFitness.hh"
#include "StatPatternRecognition/SprTransformation.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprClassifierReader.hh"
#include "StatPatternRecognition/SprAbsVarTransformer.hh"
#include "StatPatternRecognition/SprVarTransformerReader.hh"
#include "StatPatternRecognition/SprTransformerFilter.hh"

#include <stdlib.h>

//JJB
#ifndef WIN32
#include <unistd.h>
#else
int getopt(int, char**, const char*);
#endif

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <memory>
#include <algorithm>

using namespace std;


void help(const char* prog) 
{
  cout << "Usage:  " << prog 
       << " training_data_file " << endl;
  cout << "\t Options: " << endl;
  cout << "\t-h --- help                                          " << endl;
  cout << "\t-o output Tuple file                                 " << endl;
  cout << "\t-a input ascii file mode (see SprSimpleReader.hh)    " << endl;
  cout << "\t-C Size of population of Chromosomes (default = 20)  " << endl;
  cout << "\t-H Length of Gene Head (default = 10)                " << endl;
  cout << "\t-N Number of Genes per Chromosome (default=2)        " << endl;
  cout << "\t-D Number of Constants per Chromosome (default=0)    " << endl;
  cout << "\t-F List of functions to add to +-/*                  " << endl;
  cout << "\t\t Can be: A=Abs, Q=Sqrt, L=<, G=>, E=Exp, O=Log     " << endl;
  cout << "\t\t E.g -F \"AQ\"                                     " << endl;
  cout << "\t-M Maximum number of Chromosome Epochs (default=2500)" << endl;
  cout << "\t-O Mutation Rate (default = 0.05)                    " << endl;
  cout << "\t-Q RIS Rate (default = 0.1)                          " << endl;
  cout << "\t-R IS Rate (default = 0.1)                           " << endl;
  cout << "\t-S One Point Rate (default = 0.3)                    " << endl;
  cout << "\t-T Two Point Rate (default = 0.3)                    " << endl;
  cout << "\t-U Whole Gene Rate (default = 0.1)                   " << endl;
  cout << "\t-W Constant Mutation Rate (default = 0.1)            " << endl;
  cout << "\t-X Constant Swap Rate (default = 0.1)                " << endl;
  cout << "\t-Y Constant Range (default = 10.0)                   " << endl;
  cout << "\t-E Fraction of Heads used for Functions (default=0.8)" << endl;
  cout << "\t-c criterion for optimization                        " << endl;
  cout << "\t\t 1 = correctly classified fraction                 " << endl;
  cout << "\t\t 2 = signal significance s/sqrt(s+b)               " << endl;
  cout << "\t\t 3 = purity s/(s+b)                                " << endl;
  cout << "\t\t 4 = tagger efficiency Q                           " << endl;
  cout << "\t\t 5 = Gini index                                    " << endl;
  cout << "\t\t 6 = cross-entropy                                 " << endl;
  cout << "\t\t 7 = 90% Bayesian upper limit with uniform prior   " << endl;
  cout << "\t\t 8 = discovery potential 2*(sqrt(s+b)-sqrt(b))     " << endl;
  cout << "\t\t 9 = Punzi's sensitivity s/(0.5*nSigma+sqrt(b))    " << endl;
  cout << "\t\t -P background normalization factor for Punzi FOM  " << endl;
  cout << "\t\t 10= fitness = specificity*sensitivity           " << endl;
  cout << "\t-t read validation/test data from a file             " << endl;
  cout << "\t\t (must be in same format as input data!!!          " << endl;
  cout << "\t-d frequency of print-outs for validation data       " << endl;
  cout << "\t-K keep this fraction in training set and            " << endl;
  cout << "\t\t put the rest into validation set                  " << endl;
  cout << "\t-V include only these input variables                " << endl;
  cout << "\t-z exclude input variables from the list             " << endl;
  cout << "\t-Z exclude input variables from the list, "
       << "but put them in the output file " << endl;
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

  // GEP Parameters
  unsigned gene_head_length = 10;
  unsigned genes_per_chromosome = 2;
  unsigned population_size = 20;
  unsigned constants_per_chromosome = 0;
  int max_epochs = 2500;
  
  string functions;
  
  double mutation_rate = 0.05;
  double RIS_rate = 0.1;
  double IS_rate = 0.1;
  double One_Point_rate = 0.3;
  double Two_Point_rate = 0.3;
  double Whole_Gene_rate = 0.1;
  double Constant_mutation = 0.1;
  double Constant_swap = 0.1;
  double Constant_range = 10.0;
  double fractionOfHeadsForFuncts = 0.8;

  // init
  string tupleFile;
  int readMode = 0;
  SprRWFactory::DataType writeMode = SprRWFactory::Root;
  int verbose = 0;
  string outFile;
  string valFile;
  unsigned valPrint = 0;
  int iCrit = 0;
  double bW = 1.;
  bool scaleWeights = false;
  double sW = 1.;
  string includeList, excludeList;
  string inputClassesString;
  string stringVarsDoNotFeed;
  bool split = false;
  double splitFactor = 0;
  bool splitRandomize = false;
  string transformerFile;

  // decode command line
  int c;
  extern char* optarg;
  //  extern int optind;
  while((c = getopt(argc,argv,"C:H:N:D:F:M:T:d:c:P:O:Q:R:S:T:U:W:X:Y:E:ho:a:AE:n:l:L:I:i:jy:Q:g:seu:v:f:r:R:S:K:Dt:w:V:z:Z:")) != EOF ) {
    switch( c )
      {
      case 'h' :
	help(argv[0]);
	return 1;
      case 'C':
	population_size = atoi(optarg);
	assert(population_size > 1);
	break;
      case 'D':
	constants_per_chromosome = atoi(optarg);
	break;
      case 'H':
	gene_head_length = atoi(optarg);
	assert(gene_head_length > 1);
	break;
      case 'N':
	genes_per_chromosome = atoi(optarg);
	assert(genes_per_chromosome > 0);
	break;
      case 'F':
	functions = optarg;
	break;
      case 'M':
	max_epochs = atoi(optarg);
	break;
      case 't' :
	valFile = optarg;
	break;
      case 'd' :
	valPrint = (optarg==0 ? 0 : atoi(optarg));
	break;
      case 'c' :
        iCrit = (optarg==0 ? 1 : atoi(optarg));
        break;
      case 'P' :
	bW = (optarg==0 ? 1 : atof(optarg));
	break;
      case 'K' :
	split = true;
	splitFactor = (optarg==0 ? 0 : atof(optarg));
	break;
      case 'O' :
	mutation_rate = (optarg==0 ? 0.05 : atof(optarg));
	break;
      case 'Q' :
	RIS_rate = (optarg==0 ? 0.1 : atof(optarg));
	break;
      case 'R' :
	IS_rate = (optarg==0 ? 0.1 : atof(optarg));
	break;
      case 'S' :
	One_Point_rate = (optarg==0 ? 0.3 : atof(optarg));
	break;
      case 'T' :
	Two_Point_rate = (optarg==0 ? 0.3 : atof(optarg));
	break;
      case 'U' :
	Whole_Gene_rate = (optarg==0 ? 0.1 : atof(optarg));
	break;
      case 'W' :
	Constant_mutation = (optarg==0 ? 0.1 : atof(optarg));
	break;
      case 'X' :
	Constant_swap = (optarg==0 ? 0.1 : atof(optarg));
	break;
      case 'Y' :
	Constant_range = (optarg==0 ? 10. : atof(optarg));
	break;
      case 'E' :
	fractionOfHeadsForFuncts = (optarg==0 ? 0.8 : atof(optarg));
	break;
      case 'o' :
	tupleFile = optarg;
	break;
      case 'a' :
	readMode = (optarg==0 ? 0 : atoi(optarg));
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

  // Get training file.
  string trFile = argv[argc-1];

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

  // scale weights
  if( scaleWeights ) {
    cout << "Signal weights are multiplied by " << sW << endl;
    filter->scaleWeights(inputClasses[1],sW);
  }

  // read validation data from file
  auto_ptr<SprAbsFilter> valFilter;
  if( split && !valFile.empty() ) {
    cerr << "Unable to split training data and use validation data " 
	 << "from a separate file." << endl;
    return 2;
  }
  if( split && valPrint!=0 ) {
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
  if( !valFile.empty() && valPrint!=0 ) {
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

  // make optimization criterion
  auto_ptr<SprAbsTwoClassCriterion> crit;
  switch( iCrit )
    {
    case 1 :
      crit.reset(new SprTwoClassIDFraction);
      cout << "Optimization criterion set to "
           << "Fraction of correctly classified events " << endl;
      break;
    case 2 :
      crit.reset(new SprTwoClassSignalSignif);
      cout << "Optimization criterion set to "
           << "Signal significance S/sqrt(S+B) " << endl;
      break;
    case 3 :
      crit.reset(new SprTwoClassPurity);
      cout << "Optimization criterion set to "
           << "Purity S/(S+B) " << endl;
      break;
    case 4 :
      crit.reset(new SprTwoClassTaggerEff);
      cout << "Optimization criterion set to "
           << "Tagging efficiency Q = e*(1-2w)^2 " << endl;
      break;
    case 5 :
      crit.reset(new SprTwoClassGiniIndex);
      cout << "Optimization criterion set to "
	   << "Gini index  -1+p^2+q^2 " << endl;
      break;
    case 6 :
      crit.reset(new SprTwoClassCrossEntropy);
      cout << "Optimization criterion set to "
	   << "Cross-entropy p*log(p)+q*log(q) " << endl;
      break;
    case 7 :
      crit.reset(new SprTwoClassUniformPriorUL90);
      cout << "Optimization criterion set to "
           << "Inverse of 90% Bayesian upper limit with uniform prior" << endl;
      break;
    case 8 :
      crit.reset(new SprTwoClassBKDiscovery);
      cout << "Optimization criterion set to "
	   << "Discovery potential 2*(sqrt(S+B)-sqrt(B))" << endl;
      break;
    case 9 :
      crit.reset(new SprTwoClassPunzi(bW));
      cout << "Optimization criterion set to "
	   << "Punzi's sensitivity S/(0.5*nSigma+sqrt(B))" << endl;
      break;
    case 10:
      crit.reset(new SprTwoClassFitness);
      cout << "Optimization criterion set to "
           << "Fitness = TP/(TP+FN) * TN/(TN+FP)" << endl;
      break;
    default :
      cout << "Will use exponential loss by default." << endl;
      break;
    }

  // make a single GEP
  auto_ptr<SprGEP> gep(new SprGEP(filter.get(),crit.get(),
				  gene_head_length,
				  genes_per_chromosome,
				  constants_per_chromosome,
				  population_size,
				  max_epochs,
				  functions));
  gep->setRates(mutation_rate,
		RIS_rate,
		IS_rate,
		One_Point_rate,
		Two_Point_rate,
		Whole_Gene_rate,
		Constant_mutation,
		Constant_swap,
		Constant_range,
		fractionOfHeadsForFuncts);

  // set validation
  if( valFilter.get()!=0 && !valFilter->empty() && valPrint>0 )
    gep->setValidation(valFilter.get(),valPrint);

  if( !gep->train(verbose) ) {
    cerr << "Unable to train GEP." << endl;
    return 4;
  }

  auto_ptr<SprTrainedGEP> trainedGEP(gep->makeTrained());
  if( trainedGEP.get() == 0 ) {
    cerr << "Unable to make trained GEP." << endl;
    return 5;
  }

  // exit
  return 0;
}
