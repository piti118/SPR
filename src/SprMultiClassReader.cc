//$Id: SprMultiClassReader.cc,v 1.5 2008-05-08 19:57:43 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprMultiClassReader.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprAbsMultiClassLearner.hh"
#include "StatPatternRecognition/SprBinaryEncoder.hh"
#include "StatPatternRecognition/SprMultiClassLearner.hh"
#include "StatPatternRecognition/SprTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprClassifierReader.hh"
#include "StatPatternRecognition/SprStringParser.hh"

#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;



SprTrainedMultiClassLearner* SprMultiClassReader::readMultiClassLearner(
					   std::istream& input,
					   unsigned& nLine)
{
  // read indicator matrix
  string line;
  nLine++;
  if( !getline(input,line) ) {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }
  nLine++;
  if( !getline(input,line) ) {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }
  if( line.find(':') != string::npos )
    line.erase(0,line.find_first_of(':')+1);
  else {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }
  istringstream ist(line);
  unsigned nClasses(0), nClassifiers(0);
  ist >> nClasses >> nClassifiers;
  if( nClasses == 0 ) {
    cerr << "No classes found on line " << nLine << endl;
    return 0;
  }
  if( nClassifiers == 0 ) {
    cerr << "No classifiers found on line " << nLine << endl;
    return 0;
  }
  nLine++;
  if( !getline(input,line) ) {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }
  vector<int> mapper(nClasses);
  SprMatrix indicator(nClasses,nClassifiers,0);
  for( int i=0;i<nClasses;i++ ) {
    nLine++;
    if( !getline(input,line) ) {
      cerr << "Cannot read from line " << nLine << endl;
      return 0;
    }
    string sclass, srow;
    if( line.find(':') != string::npos ) {
      sclass = line.substr(0,line.find_first_of(':'));
      srow = line.substr(line.find_first_of(':')+1);
    }
    else {
      cerr << "Cannot read from line " << nLine << endl;
      return 0;
    }
    if( sclass.empty() ) {
      cerr << "Cannot read class on line " << nLine << endl;
      return 0;
    }
    if( srow.empty() ) {
      cerr << "Cannot read matrix row on line " << nLine << endl;
      return 0;
    }
    istringstream istclass(sclass), istrow(srow);
    istclass >> mapper[i];
    for( int j=0;j<nClassifiers;j++ )
      istrow >> indicator[i][j];
  }
  nLine++;
  if( !getline(input,line) ) {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }

  // read weights
  nLine++;
  if( !getline(input,line) ) {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }
  if( line.find(':') != string::npos )
    line.erase(0,line.find_first_of(':')+1);
  else {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }
  vector<double> weights(nClassifiers,1.);
  istringstream istw(line);
  for( int j=0;j<nClassifiers;j++ )
    istw >> weights[j];
  
  // read trained classifiers
  vector<pair<const SprAbsTrainedClassifier*,bool> > classifiers(nClassifiers);
  for( int n=0;n<nClassifiers;n++ ) {
    // read index of the current classifier
    nLine++;
    if( !getline(input,line) ) {
      cerr << "Cannot read from line " << nLine << endl;
      return 0;
    }
    if( line.find(':') != string::npos )
      line.erase(0,line.find_first_of(':')+1);
    else {
      cerr << "Cannot read from line " << nLine << endl;
      return 0;
    }
    istringstream istc(line);
    unsigned iClassifiers = 0;
    istc >> iClassifiers;
    if( iClassifiers != n ) {
      cerr << "Wrong classifier index on line " << nLine << endl;
      return 0;
    }

    // read each classifier
    string requested;
    SprAbsTrainedClassifier* trained =
      SprClassifierReader::readTrainedFromStream(input,requested,nLine);
    if( trained == 0 ) {
      cerr << "Unable to read trained classifier " << n << endl;
      return 0;
    }

    // add classifier to the list
    classifiers[n] = pair<const SprAbsTrainedClassifier*,bool>(trained,true);
  }// end of loop over classifiers

  // make classifier
  SprTrainedMultiClassLearner* t 
    = new SprTrainedMultiClassLearner(indicator,mapper,classifiers);

  // exit
  return t;
}


SprTrainedBinaryEncoder* SprMultiClassReader::readBinaryEncoder(
					              std::istream& input,
						      unsigned& nLine)
{
  // read classes
  string line;
  nLine++;
  if( !getline(input,line) ) {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }
  if( line.find(':') != string::npos )
    line.erase(0,line.find_first_of(':')+1);
  else {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }
  istringstream istclass(line);
  unsigned nClasses = 0;
  istclass >> nClasses;
  if( nClasses == 0 ) {
    cerr << "No classes found for BinaryEncoder on line " << nLine << endl;
    return 0;
  }
  nLine++;
  if( !getline(input,line) ) {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }
  vector<int> mapper(nClasses);
  istringstream ist(line);
  for( int i=0;i<nClasses;i++ )
    ist >> mapper[i];

  // read binary classifier
  string requested;
  SprAbsTrainedClassifier* binary 
    = SprClassifierReader::readTrainedFromStream(input,requested,nLine);
  if( binary == 0 ) {
    cerr << "Unable to read binary classifier for BinaryEncoder." << endl;
    return 0;
  }

  // make trained classifier
  SprTrainedBinaryEncoder* t = 
    new SprTrainedBinaryEncoder(mapper,binary,true);

  // exit
  return t;
}


bool SprMultiClassReader::readIndicatorMatrix(const char* filename, 
					      SprMatrix& indicator)
{
  // open file
  string fname = filename;
  ifstream input(fname.c_str());
  if( !input ) {
    cerr << "Unable to open file " << fname.c_str() << endl;
    return false;
  }
  cout << "Reading indicator matrix from file " << fname.c_str() << endl;

  // read indicator matrix dimensionality
  unsigned N(0), M(0);
  string line;
  unsigned nLine = 0;
  while( getline(input,line) ) {
    // update line counter
    nLine++;

    // remove comments
    if( line.find('#') != string::npos )
      line.erase( line.find_first_of('#') );

    // skip empty line
    if( line.find_first_not_of(' ') == string::npos ) continue;

    // make stream
    istringstream ist(line);

    // read matrix dimensions
    ist >> N >> M;
    break;
  }
  if( N==0 || M==0 ) {
    cerr << "Unable to read indicator matrix dimensionality: " 
	 << N << " " << M << "    on line " << nLine << endl;
    return false;
  }

  // read the matrix itself
  SprMatrix temp(N,M,0);
  int n = 0;
  while( n < N ) {
    nLine++;
    if( !getline(input,line) ) {
      cerr << "Unable to read line " << nLine << endl;
      return false;
    }
    if( line.find('#') != string::npos )
      line.erase( line.find_first_of('#') );
    if( line.find_first_not_of(' ') == string::npos ) continue;
    istringstream ist(line);
    for( int m=0;m<M;m++ ) ist >> temp[n][m];
    n++;
  }

  // check columns of indicator matrix
  for( int m=0;m<M;m++ ) {
    unsigned countPlus(0), countMinus(0);
    for( int n=0;n<N;n++ ) {
      int elem = int(temp[n][m]);
      if(      elem == -1 )
	countMinus++;
      else if( elem == +1 )
	countPlus++;
      else if( elem != 0 ) {
	cerr << "Invalid indicator matrix element [" << n+1 << "]" 
	     << "[" << m+1 << "]=" << elem << endl;
	return false;
      }
    }
    if( countPlus==0 || countMinus==0 ) {
      cerr << "Column " << m+1 << " of the indicator matrix does not " 
	   << "have background and signal labels present." << endl;
      return false;
    }
  }

  // check rows
  for( int n=0;n<N;n++ ) {
    unsigned sum = 0;
    for( int m=0;m<M;m++ )
      sum += abs(int(temp[n][m]));
    if( sum == 0 ) {
      cerr << "Row " << n+1 << " of the indicator matrix has nothing "
	   << "but zeros." << endl;
      return false;
    }
  }

  // exit
  indicator = temp;
  return true;
}


SprAbsTrainedMultiClassLearner* SprMultiClassReader::readTrainedFromStream(
                                               std::istream& input,
                                               const std::string& requested,
                                               unsigned& nLine)
{
  // read clasifier name
  nLine++;
  string found = SprClassifierReader::readClassifierName(input);
  if( found.empty() ) {
    cerr << "Unable to read classifier name on line " << nLine << endl;
    return 0;
  }

  // if requested classifier is supplied, make sure it matches
  if( !requested.empty() && (requested!=found) ) {
    cerr << "Requested classifier " << requested.c_str() 
         << " does not match to the actual stored classifier " 
         << found.c_str() << " on line " << nLine << endl;
    return 0;
  }

  // read specific classifier
  return SprMultiClassReader::readSelectedTrained(input,found,nLine);
}


SprAbsTrainedMultiClassLearner* SprMultiClassReader::readSelectedTrained(
                                               std::istream& input,
                                               const std::string& requested,
                                               unsigned& nLine)
{
  // switch between classifier types
  if(      requested == "MultiClassLearner" ) {
    return SprMultiClassReader::readMultiClassLearner(input,nLine);
  }
  else if( requested == "BinaryEncoder" ) {
    return SprMultiClassReader::readBinaryEncoder(input,nLine);
  }
  else {
    cerr << "Unknown multiclass learner requested." << endl;
    return 0;
  }

  // exit
  return 0;
}


SprAbsTrainedMultiClassLearner* SprMultiClassReader::readTrained(
                                                         const char* filename,
							 int verbose)
{
  // open file
  string fname = filename;
  ifstream file(fname.c_str());
  if( !file ) {
    cerr << "Unable to open file " << fname.c_str() << endl;
    return 0;
  }
  if( verbose > 0 ) {
    cout << "Reading classifier configuration from file " 
         << fname.c_str() << endl;
  }

  // read
  return SprMultiClassReader::readTrained(file,verbose);
}


SprAbsTrainedMultiClassLearner* SprMultiClassReader::readTrained(
							 std::istream& input,
							 int verbose)
{
  // make empty classifier name
  string requested;

  // start line counter
  unsigned nLine = 0;

  // read
  SprAbsTrainedMultiClassLearner* t 
    = SprMultiClassReader::readTrainedFromStream(input,requested,nLine);
  if( t == 0 ) return 0;

  // read missing values
  SprCut validRange;
  vector<SprAbsTrainedMultiClassLearner::MissingType> defaultMissing;
  if( !SprMultiClassReader::readMissing(input,validRange,
					defaultMissing,nLine) ) {
    cerr << "Unable to read missing values for MultiClassLearner." << endl;
    delete t;
    return 0;
  }

  // set missing values
  if( !t->setDefaultMissing(validRange,defaultMissing) ) {
    cerr << "Unable to set missing values for MultiClassLearner." << endl;
    delete t;
    return 0;
  }

  // read variables
  vector<string> vars;
  if( !SprClassifierReader::readVars(input,vars,nLine) ) {
    cerr << "Unable to read variables for MultiClassLearner." << endl;
    delete t;
    return 0;
  }

  // set variables
  t->setVars(vars);

  // exit
  return t;
}


bool SprMultiClassReader::readMissing(
     std::istream& input,
     SprCut& validRange,
     std::vector<SprAbsTrainedMultiClassLearner::MissingType>& defaultMissing, 
     unsigned& nLine)
{
  // read '=='
  string line;
  nLine++;
  if( !getline(input,line) ) {
    cerr << "Cannot read from line " << nLine << endl;
    return false;
  }

  // read valid range
  nLine++;
  if( !getline(input,line) ) {
    cerr << "Cannot read from line " << nLine << endl;
    return 0;
  }
  if( line.find(':') != string::npos )
    line.erase(0,line.find_first_of(':')+1);
  else {
    cerr << "Cannot read from line " << nLine << endl;
    return false;
  }
  istringstream istrange(line);
  unsigned nIntervals = 0;
  istrange >> nIntervals;
  validRange.clear();
  validRange.resize(nIntervals);
  for( int i=0;i<nIntervals;i++ )
    istrange >> validRange[i].first >> validRange[i].second;

  // read default values
  nLine++;
  if( !getline(input,line) ) {
    cerr << "Cannot read from line " << nLine << endl;
    return false;
  }
  if( line.find(':') != string::npos )
    line.erase(0,line.find_first_of(':')+1);
  else {
    cerr << "Cannot read from line " << nLine << endl;
    return false;
  }
  istringstream istdefval(line);
  unsigned nDefval = 0;
  istdefval >> nDefval;
  defaultMissing.clear();
  defaultMissing.resize(nDefval);
  for( int i=0;i<nDefval;i++ ) {
    nLine++;
    if( !getline(input,line) ) {
      cerr << "Cannot read from line " << nLine << endl;
      return false;
    }
    istringstream istclass(line);
    string dummy;
    unsigned nValForClass = 0;
    istclass >> dummy >> defaultMissing[i].first >> dummy >> nValForClass;
    if( nValForClass == 0 ) {
      cerr << "No default values found in readMissing for class " 
	   << defaultMissing[i].first << " on line " << nLine << endl;
      return false;
    }
    istclass >> dummy;
    defaultMissing[i].second.resize(nValForClass);
    for( int j=0;j<nValForClass;j++ )
      istclass >> defaultMissing[i].second[j];
  }

  // exit
  return true;
}


bool SprMultiClassReader::makeIndicatorMatrixByClass(
				const std::vector<int>& mapper,
				const std::vector<std::string>& classGroup,
				SprMatrix& indicator)
{
  // get the number of rows and columns
  int nrow = mapper.size();
  int ncol = classGroup.size();

  // sanity check
  if( nrow < 2 ) {
    cerr << "Less than 2 classes are specified." << endl;
    return false;
  }

  // make sure all classes are distinct
  for( int i=0;i<nrow;i++ ) {
    for( int j=i+1;j<nrow;j++ ) {
      if( mapper[i] == mapper[j] ) {
        cerr << "Elements " << i << " and " << j
             << " of the input vector of classes are equal." << endl;
        return false;
      }
    }
  }

  // make a new matrix
  SprMatrix mat(nrow,ncol,0);  

  // loop thru class groups
  for( int j=0;j<ncol;j++ ) {
    // get the groups
    vector<vector<int> > twoGroups;
    SprStringParser::parseToInts(classGroup[j].c_str(),twoGroups);
    if( twoGroups.size()<2 || twoGroups[0].empty() || twoGroups[1].empty() ) {
      cerr << "Unable to decode class groups in makeIndicatorMatrixByClass." 
	   << endl;
      return false;
    }

    // analyze them
    for( int k=0;k<2;k++ ) {
      for( int l=0;l<twoGroups[k].size();l++ ) {
	vector<int>::const_iterator found 
	  = find(mapper.begin(),mapper.end(),twoGroups[k][l]);
	if( found == mapper.end() ) {
	  cerr << "Class " << twoGroups[k][l] 
	       << " not found in the list of classes for "
	       << "makeIndicatorMatrixByClass" << endl;
	  return false;
	}
	int i = found - mapper.begin();
	mat[i][j] = ( k==0 ? -1 : 1 );
      }
    }
  }

  // assign matrix
  indicator = mat;

  // exit
  return true;
}


SprTrainedMultiClassLearner* SprMultiClassReader::readBinaryList(
					     const char* filename,
					     int verbose)
{
  // open file
  string fname = filename;
  ifstream file(fname.c_str());
  if( !file ) {
    cerr << "Unable to open file " << fname.c_str() << endl;
    return 0;
  }
  if( verbose > 0 ) {
    cout << "Reading classifier configuration from file " 
         << fname.c_str() << endl;
  }

  // read
  return SprMultiClassReader::readBinaryList(file,verbose);
}


SprTrainedMultiClassLearner* SprMultiClassReader::readBinaryList(
					      std::istream& input,
					      int verbose)
{
  // init vars
  SprMatrix indicator;
  vector<string> classGroups;
  vector<pair<const SprAbsTrainedClassifier*,bool> > classifiers;
  vector<double> weights;
  vector<int> classes;
  vector<string> vars;

  // read classes from first uncommented line
  string line;
  unsigned nLine = 0;
  while( getline(input,line) ) {
    // update line counter
    nLine++;

    // remove comments
    if( line.find('#') != string::npos )
      line.erase( line.find_first_of('#') );

    // skip empty line
    if( line.find_first_not_of(' ') == string::npos ) continue;

    // make stream
    istringstream ist(line);

    // read classes
    int cls = 0;
    while( ist >> cls ) classes.push_back(cls);
    break;
  }
      
  // check classes
  if( classes.empty() ) {
    cerr << "No classes found on line " << nLine 
	 << " for readBinaryList." << endl;
    return 0;
  }

  // Read row by row.
  // First value on each row is the class groups.
  // Second value is the full path to the binary classifier file. 
  while( getline(input,line) ) {
    // update line counter
    nLine++;

    // remove comments
    if( line.find('#') != string::npos )
      line.erase( line.find_first_of('#') );

    // skip empty line
    if( line.find_first_not_of(' ') == string::npos ) continue;

    // make stream
    istringstream ist(line);

    // read class group, classifier weight and file name
    string group, filename;
    double w = 0;
    ist >> group >> w >> filename;
    if( group.empty() || filename.empty() ) {
      cerr << "Unable to read class group and filename on line " << nLine 
	   << endl;
      for( int i=0;i<classifiers.size();i++ ) delete classifiers[i].first;
      return 0;
    }

    // insert class group and weight
    classGroups.push_back(group);
    weights.push_back(w);

    // read binary classifier from given name
    const SprAbsTrainedClassifier* t 
      = SprClassifierReader::readTrained(filename.c_str(),verbose);
    if( t == 0 ) {
      cerr << "Unable to read binary classifier from file " 
	   << filename.c_str() << endl;
      for( int i=0;i<classifiers.size();i++ ) delete classifiers[i].first;
      return 0;
    }

    // record variables if first classifier and check if not
    if( classifiers.empty() ) {
      t->vars(vars);
      assert( !vars.empty() );
    }
    else {
      vector<string> newvars;
      t->vars(newvars);
      if( vars.size() != newvars.size() ) {
	cerr << "Variable list sizes do not match for readBinaryList." << endl;
	for( int i=0;i<classifiers.size();i++ ) delete classifiers[i].first;
	return 0;
      }
      for( int j=0;j<vars.size();j++ ) {
	if( vars[j] != newvars[j] ) {
	  cerr << "Variable " << j << " does not match for classifier " 
	       << classifiers.size() << endl;
	  for( int i=0;i<classifiers.size();i++ ) delete classifiers[i].first;
	  return 0;
	}
      }
    }

    // insert trained classifier
    classifiers.push_back(pair<const SprAbsTrainedClassifier*,bool>(t,true));
  }

  // make indicator matrix
  assert( classGroups.size() == classifiers.size() );
  if( !SprMultiClassReader::makeIndicatorMatrixByClass(classes,
						       classGroups,
						       indicator) ) {
    cerr << "Unable to make indicator matrix in readBinaryList." << endl;
    for( int i=0;i<classifiers.size();i++ ) delete classifiers[i].first;
    return 0;
  }

  // make trained multiclass learner
  SprTrainedMultiClassLearner* trained 
    = new SprTrainedMultiClassLearner(indicator,classes,classifiers);

  // set classifier weights
  assert( weights.size() == classifiers.size() );
  trained->setClassifierWeights(weights);

  // set vars
  trained->setVars(vars);

  // exit
  return trained;
}
