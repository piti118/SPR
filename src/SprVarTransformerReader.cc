// $Id: SprVarTransformerReader.cc,v 1.3 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprVarTransformerReader.hh"
#include "StatPatternRecognition/SprAbsVarTransformer.hh"
#include "StatPatternRecognition/SprPCATransformer.hh"
#include "StatPatternRecognition/SprInputNormalizer.hh"
#include "StatPatternRecognition/SprReplaceMissing.hh"
#include "StatPatternRecognition/SprVarTransformerSequence.hh"
#include "StatPatternRecognition/SprDefs.hh"
#include "StatPatternRecognition/SprMatrix.hh"

#include <fstream>
#include <sstream>
#include <utility>
#include <cassert>

using namespace std;


SprAbsVarTransformer* SprVarTransformerReader::read(const char* filename)
{
  // open file
  string fname = filename;
  ifstream is(fname.c_str());
  if( !is ) {
    cerr << "Unable to open file " << fname.c_str() << endl;
    return 0;
  }

  // exit
  return SprVarTransformerReader::read(is);
}


SprAbsVarTransformer* SprVarTransformerReader::read(std::istream& is)
{
  unsigned nLine = 0;
  return SprVarTransformerReader::read(is,nLine);
}


SprAbsVarTransformer* SprVarTransformerReader::read(std::istream& is, 
						    unsigned& nLine)
{
  // init
  string line;
  
  // read transformer name
  nLine++;
  if( !getline(is,line) ) {
    cerr << "Unable to read VarTransformer from line " << nLine << endl;
    return 0;
  }
  istringstream ist(line);
  string dummy, transformerName, version;
  ist >> dummy >> transformerName >> version;
  
  // decode name
  if( transformerName.empty() ) {
    cerr << "Unable to read VarTransformer name on line " << nLine << endl;
    return false;
  }
  SprAbsVarTransformer* t = 0;
  if(      transformerName == "PCA" )
    t = SprVarTransformerReader::readPCATransformer(is,nLine);
  else if( transformerName == "InputNormalizer" )
    t = SprVarTransformerReader::readInputNormalizer(is,nLine);
  else if( transformerName == "ReplaceMissing" )
    t = SprVarTransformerReader::readReplaceMissing(is,nLine);
  else if( transformerName == "TransformerSequence" )
    t = SprVarTransformerReader::readTransformerSequence(is,nLine);
  else {
    cerr << "Unknown VarTransformer name specified on line " << nLine << endl;
    return 0;
  }
  if( t == 0 ) return 0;
  
  // read vars
  vector<string> oldVars, newVars;
  if( !SprVarTransformerReader::readVars(is,nLine,oldVars,newVars) ||
      oldVars.empty() || newVars.empty() ) {
    cerr << "Unable to read VarTransformer variables." << endl;
    return 0;
  }
  t->setOldVars(oldVars);
  t->setNewVars(newVars);
  
  // exit
  return t;
}


bool SprVarTransformerReader::readVars(std::istream& is, unsigned& nLine,
				       std::vector<std::string>& oldVars,
				       std::vector<std::string>& newVars)
{
  // read old variables
  oldVars.clear();

  // skip 2 lines
  string line;
  for( int i=0;i<2;i++ ) {
    nLine++;
    if( !getline(is,line) ) {
      cerr << "Unable to read VarTransformer from line " << nLine << endl;
      return false;
    }
  }

  // read all lines skipping those that have nothing but =
  while( getline(is,line) ) {
    nLine++;

    // get rid of spaces
    line.erase( 0, line.find_first_not_of(' ') );
    line.erase( line.find_last_not_of(' ')+1 );

    // get rid of '='
    line.erase( 0, line.find_first_not_of('=') );
    line.erase( line.find_last_not_of('=')+1 );

    // if empty, stop reading
    if( line.empty() ) break;

    // add var
    istringstream ist(line);
    int index = -1;
    string var;
    ist >> index >> var;
    if( index != oldVars.size() ) {
      cerr << "Incorrect VarTransformer variable index on line " 
	   << nLine << endl;
      return false;
    }
    oldVars.push_back(var);
  }

  // read old variables
  newVars.clear();

  // skip 2 lines
  for( int i=0;i<2;i++ ) {
    nLine++;
    if( !getline(is,line) ) {
      cerr << "Unable to read VarTransformer from line " << nLine << endl;
      return false;
    }
  }

  // read all lines skipping those that have nothing but =
  while( getline(is,line) ) {
    nLine++;

    // get rid of spaces
    line.erase( 0, line.find_first_not_of(' ') );
    line.erase( line.find_last_not_of(' ')+1 );

    // get rid of '='
    line.erase( 0, line.find_first_not_of('=') );
    line.erase( line.find_last_not_of('=')+1 );

    // if empty, do stop reading
    if( line.empty() ) break;

    // add var
    istringstream ist(line);
    int index = -1;
    string var;
    ist >> index >> var;
    if( index != newVars.size() ) {
      cerr << "Incorrect VarTransformer variable index on line " 
	   << nLine << endl;
      return false;
    }
    newVars.push_back(var);
  }

  // exit
  return true;
}


SprPCATransformer* SprVarTransformerReader::readPCATransformer(
					    std::istream& is, unsigned& nLine)
{
  // read dimensionality
  string line;
  nLine++;
  if( !getline(is,line) ) {
    cerr << "Unable to read VarTransformer from line " << nLine << endl;
    return 0;
  }
  istringstream ist(line);
  string dummy;
  int dim = -1;
  ist >> dummy >> dim;
  if( dim <= 0 ) {
    cerr << "Unable to read dimensionality from VarTransformer line " 
	 << nLine << endl;
    return 0;
  }
  
  // read eigenvalues
  vector<pair<double,int> > eigenValues(dim);
  nLine++;
  if( !getline(is,line) ) {
    cerr << "Unable to read VarTransformer from line " << nLine << endl;
    return 0;
  }
  istringstream ist_eigen(line);
  ist_eigen >> dummy;
  for( int d=0;d<dim;d++ )
    ist_eigen >> eigenValues[d].first;
  nLine++;
  if( !getline(is,line) ) {
    cerr << "Unable to read VarTransformer from line " << nLine << endl;
    return 0;
  }
  istringstream ist_index(line);
  ist_index >> dummy;
  for( int d=0;d<dim;d++ ) {
    ist_index >> eigenValues[d].second;
    assert( eigenValues[d].second >= 0 );
  }

  // read transformation matrix
  SprMatrix U(dim,dim);
  for( int i=0;i<dim;i++ ) {
    nLine++;
    if( !getline(is,line) ) {
      cerr << "Unable to read VarTransformer from line " << nLine << endl;
      return 0;
    }
    istringstream istU(line);
    int d = -1;
    istU >> d;
    if( d != i ) {
      cerr << "Dimension of VarTransformer does not macth on line " 
	   << nLine << endl;
      return 0;
    }
    for( int j=0;j<dim;j++ )
      istU >> dummy >> dummy >> U[i][j];
  }

  // make PCA transformer
  SprPCATransformer* t = new SprPCATransformer(U,eigenValues);

  // exit
  return t;
}


SprInputNormalizer* SprVarTransformerReader::readInputNormalizer(
					    std::istream& is, unsigned& nLine)
{
  // read dimensionality
  string line;
  nLine++;
  if( !getline(is,line) ) {
    cerr << "Unable to read VarTransformer from line " << nLine << endl;
    return 0;
  }
  istringstream ist(line);
  string dummy;
  int dim = -1;
  ist >> dummy >> dim;
  if( dim <= 0 ) {
    cerr << "Unable to read dimensionality from VarTransformer line " 
	 << nLine << endl;
    return 0;
  }
  
  // read means and sigmas
  vector<double> mean(dim), sigma(dim);
  for( int i=0;i<dim;i++ ) {
    nLine++;
    if( !getline(is,line) ) {
      cerr << "Unable to read VarTransformer from line " << nLine << endl;
      return 0;
    }
    istringstream istms(line);
    int d = -1;
    istms >> d;
    if( d != i ) {
      cerr << "Dimension of VarTransformer does not macth on line " 
	   << nLine << endl;
      return 0;
    }
    istms >> dummy >> mean[i] >> dummy >> sigma[i];
  }

  // make PCA transformer
  SprInputNormalizer* t = new SprInputNormalizer(mean,sigma);

  // exit
  return t;
}


SprReplaceMissing* SprVarTransformerReader::readReplaceMissing(
					  std::istream& is, unsigned& nLine)
{
  // read blindness
  string line;
  nLine++;
  if( !getline(is,line) ) {
    cerr << "Unable to read VarTransformer from line " << nLine << endl;
    return 0;
  }
  string dummy;
  int classBlind = -1;
  istringstream istblind(line);
  istblind >> dummy >> classBlind;
  if( classBlind!=0 && classBlind!=1 ) {
    cerr << "Unable to read classBlind from VarTransformer line " 
	 << nLine << endl;
    return 0;
  }

  // read valid range
  nLine++;
  if( !getline(is,line) ) {
    cerr << "Unable to read VarTransformer from line " << nLine << endl;
    return 0;
  }
  int nValid = 0;
  istringstream istvalid(line);
  istvalid >> dummy >> nValid;
  if( nValid <= 0 ) {
    cerr << "Unable to read size of valid range from VarTransformer line " 
	 << nLine << endl;
    return 0;
  }
  vector<SprInterval> validRange(nValid);
  for( int i=0;i<nValid;i++ ) {
    nLine++;
    if( !getline(is,line) ) {
      cerr << "Unable to read VarTransformer from line " << nLine << endl;
      return 0;
    }
    int index = -1;
    istringstream istpair(line);
    istpair >> index >> validRange[i].first >> validRange[i].second;
    if( index!=i || validRange[i].first>validRange[i].second ) {
      cerr << "Unable to read valid range interval from VarTransformer line " 
	   << nLine << endl;
      return 0;
    }
  }

  // read replacement values
  nLine++;
  if( !getline(is,line) ) {
    cerr << "Unable to read VarTransformer from line " << nLine << endl;
    return 0;
  }
  int nClasses = 0;
  istringstream istclasses(line);
  istclasses >> dummy >> nClasses;
  if( nClasses <= 0 ) {
    cerr << "Unable to read number of classes from VarTransformer line " 
	 << nLine << endl;
    return 0;
  }
  vector<SprReplaceMissing::ClassAndDefaultValues> replacement(nClasses);
  for( int ic=0;ic<nClasses;ic++ ) {
    nLine++;
    if( !getline(is,line) ) {
      cerr << "Unable to read VarTransformer from line " << nLine << endl;
      return 0;
    }
    istringstream istclass(line);
    int dim = 0;
    istclass >> dummy >> replacement[ic].first >> dummy >> dim;
    if( dim <= 0 ) {
      cerr << "Unable to read class dimensionality from VarTransformer line " 
	   << nLine << endl;
      return 0;
    }
    replacement[ic].second.resize(dim);
    for( int j=0;j<dim;j++ ) {
      nLine++;
      if( !getline(is,line) ) {
	cerr << "Unable to read VarTransformer from line " << nLine << endl;
	return 0;
      }
      istringstream istrep(line);
      int index = -1;
      istrep >> index >> replacement[ic].second[j];
      if( index != j ) {     
	cerr << "Unable to read replacement value from VarTransformer line " 
	     << nLine << endl;
	return 0;
      }
    }// end loop over dimensions
  }// end loop over classes

  // make transformer
  SprReplaceMissing* t 
    = new SprReplaceMissing(validRange,(classBlind==1));
  t->setReplacement(replacement);

  // exit
  return t;
}


SprVarTransformerSequence* SprVarTransformerReader::readTransformerSequence(
					    std::istream& is, unsigned& nLine)
{
  // read transformer names
  string line;
  nLine++;
  if( !getline(is,line) ) {
    cerr << "Unable to read VarTransformer from line " << nLine << endl;
    return 0;
  }
  istringstream istn(line);
  unsigned nTransformers = 0;
  istn >> nTransformers;
  if( nTransformers == 0 ) {
    cerr << "Unable to read the number of transformers on line " 
	 << nLine << endl;
    return 0;
  }
  vector<string> transformerNames(nTransformers);
  for( int n=0;n<nTransformers;n++ )
    istn >> transformerNames[n];
  
  // read transformers (all owned by the sequencer)
  vector<pair<SprAbsVarTransformer*,bool> > 
    transformers(nTransformers,pair<SprAbsVarTransformer*,bool>(0,true));
  for( int n=0;n<nTransformers;n++ ) {
    SprAbsVarTransformer* t = SprVarTransformerReader::read(is,nLine);
    if( t == 0 ) return 0;
    transformers[n].first = t;
  }
    
  // make the sequencer
  SprVarTransformerSequence* seq 
    = new SprVarTransformerSequence(transformers);

  // exit
  return seq;
}
