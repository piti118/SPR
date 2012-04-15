// $Id: SprReplaceMissing.cc,v 1.1 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh" 
#include "StatPatternRecognition/SprReplaceMissing.hh" 
#include "StatPatternRecognition/SprAbsFilter.hh" 
#include "StatPatternRecognition/SprPoint.hh" 
#include "StatPatternRecognition/SprUtils.hh" 

#include <cassert>
#include <algorithm>
#include <functional>

using namespace std;


struct SRMCmpPairFirst
  : public binary_function<pair<double,double>,pair<double,double>,bool> {
  bool operator()(const pair<double,double>& l, const pair<double,double>& r)
    const {
    return (l.first < r.first);
  }
};


SprReplaceMissing::SprReplaceMissing()
  :
  SprAbsVarTransformer(),
  mode_(Median),
  validRange_(),
  classBlind_(true),
  replacement_(),
  chosenClass_(0)
{}


SprReplaceMissing::SprReplaceMissing(
		    Mode mode,
		    const SprCut& validRange,
		    bool classBlind)
  :
  SprAbsVarTransformer(),
  mode_(mode),
  validRange_(validRange),
  classBlind_(classBlind),
  replacement_(),
  chosenClass_(0)
{}


SprReplaceMissing::SprReplaceMissing(
		    const SprCut& validRange,
		    bool classBlind)
  :
  SprAbsVarTransformer(),
  mode_(Median),
  validRange_(validRange),
  classBlind_(classBlind),
  replacement_(),
  chosenClass_(0)
{}


SprReplaceMissing::SprReplaceMissing(const SprReplaceMissing& other)
  :
  SprAbsVarTransformer(other),
  mode_(other.mode_),
  validRange_(other.validRange_),
  classBlind_(other.classBlind_),
  replacement_(other.replacement_),
  chosenClass_(other.chosenClass_)
{}


bool SprReplaceMissing::allVarsIndependent() const
{
  if( oldVars_.size() != newVars_.size() ) {
    cerr << "Variable lists sizes do not match in "
	 << "SprReplaceMissing::allVarsIndependent()" << endl;
    return false;
  }
  for( int i=0;i<oldVars_.size();i++ ) {
    if( oldVars_[i] != newVars_[i] ) {
      cerr << "Variable mismatch in SprReplaceMissing::allVarsIndependent(): "
	   << oldVars_[i] << " " << newVars_[i] << endl;
      return false;
    }
  }
  return true;
}


bool SprReplaceMissing::reduceVars(const std::vector<std::string>& vars)
{
  // sanity check
  if( !this->allVarsIndependent() ) {
    cerr << "SprReplaceMissing cannot reduce variable list. "
	 << "Independence check fails." << endl;
    return false;
  }

  // make sure dimensionality matches everywhere
  int nClass = replacement_.size();
  for( int ic=0;ic<nClass;ic++ )
    assert( replacement_[ic].second.size() == oldVars_.size() );

  // prepare new vectors
  vector<string> newVars;
  vector<ClassAndDefaultValues> replacement(nClass);
  for( int ic=0;ic<nClass;ic++ )
    replacement[ic].first = replacement_[ic].first;

  // keep only supplied vars
  for( int d=0;d<oldVars_.size();d++ ) {
    if( find(vars.begin(),vars.end(),oldVars_[d]) != vars.end() ) {
      newVars.push_back(oldVars_[d]);
      for( int ic=0;ic<nClass;ic++ )
	replacement[ic].second.push_back(replacement_[ic].second[d]);
    }
  }

  // reassign
  oldVars_ = newVars;
  newVars_ = newVars;
  replacement_ = replacement;

  // exit
  return true;
}

bool SprReplaceMissing::train(const SprAbsFilter* data, int verbose)
{
  // get vars
  data->vars(oldVars_);
  assert( !oldVars_.empty() );
  newVars_ = oldVars_;
  int dim = oldVars_.size();

  // init
  replacement_.clear();

  // sanity check
  if( validRange_.empty() ) return true;

  // get classes
  vector<SprClass> classes;
  data->allClasses(classes);
  int nClasses = classes.size();
  assert( nClasses > 0 );

  // prepare vector
  replacement_.resize(nClasses);
  for( int ic=0;ic<nClasses;ic++ ) {
    replacement_[ic] = ClassAndDefaultValues(classes[ic],
					     vector<double>(dim,0));
  }

  // loop thru dimensions
  for( int d=0;d<dim;d++ ) {

    // vector of valid values used for replacement computation
    // first number in a pair is coordinate, 2nd number is weight
    vector<ValuesWithWeights> in_range(nClasses);

    // points outside the valid range
    vector<SprPoint*> out_range;

    // split data into points with and without missing values
    for( int i=0;i<data->size();i++ ) {
      SprPoint* p = (*data)[i];
      bool valid = false;
      double r = (p->x_)[d];
      for( int j=0;j<validRange_.size();j++ ) {
	if( r>validRange_[j].first && r<validRange_[j].second ) {
	  valid = true;
	  break;
	}
      }
      if( valid ) {
	if( classBlind_ ) {
	  for( int ic=0;ic<nClasses;ic++ )
	    in_range[ic].push_back(pair<double,double>(r,data->w(i)));
	}
	else {
	  vector<SprClass>::const_iterator found 
	    = ::find(classes.begin(),classes.end(),p->class_);
	  assert( found != classes.end() );
	  int ic = found - classes.begin();
	  in_range[ic].push_back(pair<double,double>(r,data->w(i)));
	}
      }
      else {
	out_range.push_back(p);
      }
    }

    // loop thru classes to compute replacement
    for( int ic=0;ic<nClasses;ic++ ) {
      // print out
      if( verbose > 1 ) {
	cout << "Found " << in_range[ic].size() << " valid values for class "
	     << classes[ic] << endl;
      }

      // get vector and class
      const SprClass& cls = classes[ic];
      ValuesWithWeights& in = in_range[ic];

      // sort valid values
      int validSize = in.size();
      if( validSize < 1 ) {
	cerr << "No points with valid values found for variable " 
	     << oldVars_[d].c_str() << " for class " << cls << endl;
	return false;
      }
      stable_sort(in.begin(),in.end(),SRMCmpPairFirst());

      // compute cumulative weights
      double wtot = 0;
      vector<double> w(validSize);
      for( int i=0;i<validSize;i++ ) {
	w[i] = in[i].second;
	wtot += w[i];
      }
      if( wtot < SprUtils::eps() ) {
	cerr << "Total weight of valid values found for variable " 
	     << oldVars_[d].c_str() << " for class " << cls << endl;
	return false;
      }

      // compute replacements
      double forReplacement = 0;
      if(      mode_ == Median ) {
	// get index
	int medIndex = -1;
	w[0] /= wtot;
	if( w[0] > 0.5 ) {
	  medIndex = 0;
	}
	else {
	  for( int i=1;i<validSize;i++ ) {
	    w[i] = w[i-1] + w[i]/wtot;
	    if( w[i] > 0.5 ) {
	      medIndex = i;
	      break;
	    }
	  }
	}
	if( medIndex<0 || medIndex>(validSize-1)) {
	  cerr << "Unable to find median index in dimension " 
	       << d << " for class " << cls << endl;
	  return false;
	}

	// compute median
	double med = 0;
	if( medIndex == 0 ) 
	  med = in[medIndex].first;
	else
	  med = 0.5*(in[medIndex].first+in[medIndex-1].first);

	// store the median
	forReplacement = med;
      }
      else if( mode_ == Average ) {
	double ave = 0;
	for( int i=0;i<validSize;i++ )
	  ave += w[i]*in[i].first;
	ave /= wtot;

	// store the average
	forReplacement = ave;
      }

      // store the replacement value
      replacement_[ic].second[d] = forReplacement;

      // print out
      if( verbose > 0 ) {
	cout << "Computed replacement for class " << cls
	     << " in variable " << oldVars_[d].c_str() 
	     << "      with value " << forReplacement << endl;
      }
    }// end loop thru classes to compute replacements
  }// end loop thru dimensions

  // exit
  return true;
}


void SprReplaceMissing::transform(const std::vector<double>& in, 
				  std::vector<double>& out) const
{
  out = in;
  int ic = ( classBlind_ ? 0 : chosenClass_ );
  if( !replacement_.empty() ) {
    for( int d=0;d<in.size();d++ ) {
      bool valid = false;
      double r = in[d];
      for( int j=0;j<validRange_.size();j++ ) {
	if( r>validRange_[j].first && r<validRange_[j].second ) {
	  valid = true;
	  break;
	}
      }
      if( !valid )
	out[d] = replacement_[ic].second[d];
    }// end loop over d
  }// end if !replacement_.empty()
}


void SprReplaceMissing::inverse(const std::vector<double>& in, 
				std::vector<double>& out) const
{
  out = in;
}


void SprReplaceMissing::print(std::ostream& os) const
{
  // save name
  os << "VarTransformer: " << this->name().c_str() 
     << " " << SprVersion.c_str() << endl;

  // init
  int dim = oldVars_.size();
  vector<string> vars(oldVars_);

  // protect againt spaces in var names
  for( int i=0;i<dim;i++ ) {
    if( vars[i].find(' ') != string::npos )
	vars[i].erase(vars[i].find_first_of(' '));
  }

  // store blindness
  os << "ClassBlind: " << int(classBlind_) << endl;

  // store valid range
  os << "ValidRange: " << validRange_.size() << endl;
  for( int i=0;i<validRange_.size();i++ ) {
    os << i << "     " << validRange_[i].first 
       << " " << validRange_[i].second << endl;
  }

  // store classes, mode and replacement values
  os << "Classes: " << replacement_.size() 
     << " Mode: " << int(mode_) << endl;
  for( int ic=0;ic<replacement_.size();ic++ ) {
    os << "Class: " << replacement_[ic].first
       << " Size: " << replacement_[ic].second.size() << endl;
    for( int j=0;j<replacement_[ic].second.size();j++ )
      os << j << " " << replacement_[ic].second[j] << endl;
  }
}


bool SprReplaceMissing::chooseClass(const SprClass& cls)
{
  // find class
  for( int ic=0;ic<replacement_.size();ic++ ) {
    if( cls == replacement_[ic].first ) {
      chosenClass_ = ic;
      return true;
    }
  }

  // class not found
  cerr << "SprReplaceMissing cannot find requested class." << endl;
  return false;
}
