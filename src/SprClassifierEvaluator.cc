// $Id: SprClassifierEvaluator.cc,v 1.7 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprClassifierEvaluator.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprAbsClassifier.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprAbsMultiClassLearner.hh"
#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprCoordinateMapper.hh"
#include "StatPatternRecognition/SprClass.hh"
#include "StatPatternRecognition/SprAverageLoss.hh"
#include "StatPatternRecognition/SprIntegerPermutator.hh"
#include "StatPatternRecognition/SprPoint.hh"
#include "StatPatternRecognition/SprStringParser.hh"
#include "StatPatternRecognition/SprUtils.hh"
#include "StatPatternRecognition/SprVarSelectorFilter.hh"
#include "StatPatternRecognition/SprFomCalculator.hh"
#include "StatPatternRecognition/SprCrossValidator.hh"

#include <stdio.h>
#include <iostream>
#include <memory>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iomanip>
#include <map>

using namespace std;


struct SCECmpPairSDDSecondFirst
  : public binary_function<SprClassifierEvaluator::NameAndValue,
			   SprClassifierEvaluator::NameAndValue,bool> {
  bool operator()(const SprClassifierEvaluator::NameAndValue& l, 
		  const SprClassifierEvaluator::NameAndValue& r)
    const {
    return (l.second.first < r.second.first);
  }
};


bool SprClassifierEvaluator::variableImportance(
			       const SprAbsFilter* data,
			       SprAbsTrainedClassifier* trained,
			       SprAbsTrainedMultiClassLearner* mcTrained,
			       SprAverageLoss* loss,
			       SprCoordinateMapper* mapper,
			       unsigned nPerm,
			       std::vector<NameAndValue>& lossIncrease)
{
  // sanity check
  if( data == 0 ) {
    cerr << "No data supplied for variableImportance." << endl;
    return false;
  }
  if( trained==0 && mcTrained==0 ) {
    cerr << "No classifiers provided for variableImportance." << endl;
    return false;
  }
  if( trained!=0 && mcTrained!=0 ) {
    cerr << "variableImportance cannot process both two-class " 
	 << "and multi-class learners." << endl;
    return false;
  }
  if( nPerm == 0 ) {
    cout << "No permutations requested. Will use one by default." << endl;
    nPerm = 1;
  }

  // check classes
  vector<SprClass> classes; 
  data->classes(classes); 
  if( classes.size() < 2 ) {
    cerr << "Classes have not been set." << endl;
    return false; 
  }
  vector<int> mcClasses;
  if( mcTrained != 0 ) {
    mcTrained->classes(mcClasses);
    assert( !mcClasses.empty() );
  }

  // enforce normalization
  if( trained != 0 )
    trained->useNormalized();

  //
  // pass through all variables
  //
  vector<string> testVars;
  if( mcTrained != 0 )
    mcTrained->vars(testVars);
  else
    trained->vars(testVars);
  int N = data->size();
  SprIntegerPermutator permu(N);

  // make first pass without permutations
  for( int n=0;n<N;n++ ) {
    const SprPoint* p = (*data)[n];
    const SprPoint* mappedP = p;
    int icls = p->class_;
    if( mcTrained != 0 ) {
      if( find(mcClasses.begin(),mcClasses.end(),icls) == mcClasses.end() )
	continue;
    }
    else {
      if(      icls == classes[0] )
	icls = 0;
      else if( icls == classes[1] )
	icls = 1;
      else
	continue;
    }
    if( mapper != 0 ) mappedP = mapper->output(p);
    if( mcTrained != 0 )
      loss->update(icls,mcTrained->response(mappedP),data->w(n));
    else
      loss->update(icls,trained->response(mappedP),data->w(n));
    if(  mapper != 0 ) mapper->clear();
  }
  double nominalLoss = loss->value();

  //
  // loop over permutations
  //
  cout << "Will perform " << nPerm << " permutations per variable." << endl;
  int nVars = testVars.size();
  lossIncrease.clear();
  lossIncrease.resize(nVars);
  for( int d=0;d<nVars;d++ ) {
    cout << "Permuting variable " << testVars[d].c_str() << endl;

    // map this var
    int mappedD = d;
    if( mapper != 0 )
      mappedD = mapper->mappedIndex(d);
    assert( mappedD>=0 && mappedD<data->dim() );

    // pass through all points permuting them
    vector<double> losses(nPerm);
    double aveLoss = 0;
    for( int i=0;i<nPerm;i++ ) {

      // permute this variable
      vector<unsigned> seq;
      if( !permu.sequence(seq) ) {
        cerr << "variableImportance is unable to permute points." << endl;
        return false;
      }

      // pass through points
      loss->reset();
      for( int n=0;n<N;n++ ) {
        SprPoint p(*(*data)[n]);
        p.x_[mappedD] = (*data)[seq[n]]->x_[mappedD];
        const SprPoint* mappedP = &p;
        int icls = p.class_;
	if( mcTrained != 0 ) {
	  if( find(mcClasses.begin(),mcClasses.end(),icls) == mcClasses.end() )
	    continue;
	}
	else {
	  if(      icls == classes[0] )
	    icls = 0;
	  else if( icls == classes[1] )
	    icls = 1;
	  else
	    continue;
	}
        if( mapper != 0 ) mappedP = mapper->output(&p);
	if( mcTrained != 0 )
	  loss->update(icls,mcTrained->response(mappedP),data->w(n));
	else
	  loss->update(icls,trained->response(mappedP),data->w(n));
        if( mapper != 0 ) mapper->clear();
      }

      // store loss
      losses[i] = loss->value();
      aveLoss += losses[i];
    }// end loop over permutations

    // compute error
    aveLoss /= nPerm;
    double err = 0;
    for( int i=0;i<nPerm;i++ )
      err += (losses[i]-aveLoss)*(losses[i]-aveLoss);
    if( nPerm > 1 )
      err /= (nPerm-1);
    err = sqrt(err);

    // store values
    lossIncrease[d] = NameAndValue(testVars[d],
				   SprValueWithError(aveLoss-nominalLoss,err));
  }// end loop over variables

  // exit
  return true;
}


bool SprClassifierEvaluator::twoSubsetInteraction(
				  const SprAbsFilter* data,
				  SprAbsTrainedClassifier* trained,
				  SprAbsTrainedMultiClassLearner* mcTrained,
				  SprCoordinateMapper* mapper,
				  const std::vector<std::string>& subset1,
				  const std::vector<std::string>& subset2,
				  unsigned nPoints,
				  SprValueWithError& interaction,
				  int verbose)
{
  // sanity check
  if( data == 0 ) {
    cerr << "No data supplied for twoSubsetInteraction." << endl;
    return false;
  }
  if( trained==0 && mcTrained==0 ) {
    cerr << "No classifiers provided for twoSubsetInteraction." << endl;
    return false;
  }
  if( trained!=0 && mcTrained!=0 ) {
    cerr << "variableInteraction cannot process both two-class " 
	 << "and multi-class learners." << endl;
    return false;
  }
  if( nPoints > data->size() ) {
    cerr << "Number of points for integration " << nPoints 
	 << " cannot exceed the data size " << data->size() << endl;
    return false;
  }

  // set number of points and passes
  int nPass = 2;
  if( nPoints == 0 ) {
    nPoints = data->size();
    nPass = 1;
  }

  // get trained vars
  vector<string> testVars;
  if(      trained != 0 )
    trained->vars(testVars);
  else if( mcTrained != 0 )
    mcTrained->vars(testVars);
  assert( !testVars.empty() );
  unsigned dim  = testVars.size();

  // get dimensionality of subsets
  unsigned dim1 = subset1.size();
  unsigned dim2 = subset2.size();
  if( dim<=dim1 || dim<=dim2 ) {
    cerr << "Too many variables requested in one of subsets: " 
	 << dim << " " << dim1 << " " << dim2 << endl;
    return false;
  }

  // check that the subsets do not overlap
  for( int d=0;d<subset1.size();d++ ) {
    if( find(subset2.begin(),subset2.end(),subset1[d]) != subset2.end() ) {
      cerr << "Cannot process the two subsets because they have "
	   << "common variable " << subset1[d].c_str() << endl;
      return false;
    }
  }

  // map subset vars onto classifier vars
  vector<int> subsetIndex1(dim1,-1), subsetIndex2(dim2,-1);
  for( int d=0;d<dim1;d++ ) {
    vector<string>::const_iterator found 
      = find(testVars.begin(),testVars.end(),subset1[d]);
    if( found == testVars.end() ) {
      cerr << "Variable " << subset1[d].c_str() 
	   << " not found among trained variables in " 
	   << "twoSubsetInteraction." << endl;
      return false;
    }
    subsetIndex1[d] = found - testVars.begin();
  }
  for( int d=0;d<dim2;d++ ) {
    vector<string>::const_iterator found 
      = find(testVars.begin(),testVars.end(),subset2[d]);
    if( found == testVars.end() ) {
      cerr << "Variable " << subset2[d].c_str() 
	   << " not found among trained variables in " 
	   << "twoSubsetInteraction." << endl;
      return false;
    }
    subsetIndex2[d] = found - testVars.begin();
  }

  // use normalized output
  if( trained != 0 )
    trained->useNormalized();

  // init interaction
  interaction = SprValueWithError(0,0);

  // reduce data size for integration by randomly choosing
  //   the number of points defined by the user
  unsigned N = data->size();
  SprIntegerPermutator permu(N);

  //
  // two passes to get a rough error estimate
  //
  vector<double> pass(nPass);
  for( int ipass=0;ipass<nPass;ipass++ ) {
    if( verbose > 1 )
      cout << "Pass " << ipass << endl;

    // permute indices
    vector<unsigned> indices;
    if( !permu.sequence(indices) ) {
      cerr << "Unable to permute input indices." << endl;
      return false;
    }
    double wtot = 0;
    for( int i=0;i<nPoints;i++ )
      wtot += data->w(indices[i]);

    //
    // compute correlation between F(S1) and F(S2)
    //
      
    // compute mean F1 and F2 at each point
    vector<double> F1(nPoints,0), F2(nPoints,0);
    for( int i=0;i<nPoints;i++ ) {
      if( verbose>1 && (i+1)%1000==0 )
	cout << "Processing point " << i+1 << endl;
      int ii = indices[i];
      double wi = data->w(ii);
      const SprPoint* pi = (*data)[ii];
      const SprPoint* mappedPi = pi;
      if( mapper != 0 ) mappedPi = mapper->output(pi);
      const vector<double>& xi = mappedPi->x_;
      vector<double> xi_subset_1(dim1), xi_subset_2(dim2);
      for( int k=0;k<dim1;k++ )
	xi_subset_1[k] = xi[subsetIndex1[k]];
      for( int k=0;k<dim2;k++ )
	xi_subset_2[k] = xi[subsetIndex2[k]];
	
      for( int j=0;j<nPoints;j++ ) {
	int jj = indices[j];
	const SprPoint* pj = (*data)[jj];
	vector<double> x_S1(pj->x_), x_S2(pj->x_);
	double wj = data->w(jj);
	if( mapper != 0 ) {
	  mapper->map(pj->x_,x_S1);
	  mapper->map(pj->x_,x_S2);
	}
	for( int k=0;k<dim1;k++ )
	  x_S1[subsetIndex1[k]] = xi_subset_1[k];
	for( int k=0;k<dim2;k++ )
	  x_S2[subsetIndex2[k]] = xi_subset_2[k];
	if(      trained != 0 ) {
	  F1[i] += wj * trained->response(x_S1);
	  F2[i] += wj * trained->response(x_S2);
	}
	else if( mcTrained != 0 ) {
	  F1[i] += wj * mcTrained->response(x_S1);
	  F2[i] += wj * mcTrained->response(x_S2);
	}
      }// end loop over j
      F1[i] /= wtot;
      F2[i] /= wtot;
	
      // cleanup
      if( mapper != 0 ) mapper->clear();
    }// end loop over i
      
    // compute correlation between F1 and F2
    double F1_mean(0), F2_mean(0);
    for( int i=0;i<nPoints;i++ ) {
      int ii = indices[i];
      double wi = data->w(ii);
      F1_mean += wi*F1[i];
      F2_mean += wi*F2[i];
    }
    F1_mean /= wtot;
    F2_mean /= wtot;
    double var1(0), var2(0), cov(0);
    for( int i=0;i<nPoints;i++ ) {
      int ii = indices[i];
      double wi = data->w(ii);
      var1 += wi*pow(F1[i]-F1_mean,2);
      var2 += wi*pow(F2[i]-F2_mean,2);
      cov  += wi*(F1[i]-F1_mean)*(F2[i]-F2_mean);
    }
    var1 /= wtot;
    var2 /= wtot;
    cov  /= wtot;
      
    // set correlation
    if( var1<SprUtils::eps() || var2<SprUtils::eps() ) {
      cerr << "Variance too small: " << var1 << " " << var2 
	   << ". Unable to compute variable interaction." << endl;
      cout << "Subset 1:   ";
      for( int m=0;m<subset1.size();m++ ) 
	cout << " " << subset1[m].c_str();
      cout << endl;
      cout << "Subset 2:   ";
      for( int m=0;m<subset2.size();m++ ) 
	cout << " " << subset2[m].c_str();
      cout << endl;
      pass[ipass] = 0;
    }
    else      
      pass[ipass] = cov/(sqrt(var1)*sqrt(var2));
  }// end loop over ipass

  //
  // compute average over passes
  //
  double mean = 0;
  for( int ipass=0;ipass<nPass;ipass++ ) mean += pass[ipass];
  mean /= nPass;
  double var = 0;
  for( int ipass=0;ipass<nPass;ipass++ )
    var += pow(pass[ipass]-mean,2);
  if( nPass > 1 )
    var /= (nPass-1);
  double sigma = ( var>0 ? sqrt(var) : 0 );
  interaction = SprValueWithError(mean,sigma);

  // exit
  return true;
}


bool SprClassifierEvaluator::variableInteraction(
				  const SprAbsFilter* data,
				  SprAbsTrainedClassifier* trained,
				  SprAbsTrainedMultiClassLearner* mcTrained,
				  SprCoordinateMapper* mapper,
				  const char* vars,
				  unsigned nPoints,
				  std::vector<NameAndValue>& interaction,
				  int verbose)
{
  // sanity check
  if( trained==0 && mcTrained==0 ) {
    cerr << "No classifiers specified for variableInteraction." << endl;
    return false;
  }

  // init
  interaction.clear();

  // decode input string
  vector<vector<string> > v_of_vars;
  SprStringParser::parseToStrings(vars,v_of_vars);

  // determine subsets
  const vector<string> *subset_1(0), *subset_2(0);
  if( !v_of_vars.empty() ) {
    if( !v_of_vars[0].empty() )
      subset_1 = &v_of_vars[0];
    if( !v_of_vars[1].empty() )
      subset_2 = &v_of_vars[1];
  }

  // make sure subset_1 gets filled first
  if( subset_1==0 && subset_2!=0 ) {
    subset_1 = subset_2;
    subset_2 = 0;
  }

  // get trained vars
  vector<string> testVars;
  if(      trained != 0 )
    trained->vars(testVars);
  else if( mcTrained != 0 )
    mcTrained->vars(testVars);
  assert( !testVars.empty() );

  // if one of subsets is set to ".", replace with all variables
  // from testVars not in the other subset
  vector<string> complementary;
  if( subset_1!=0 && (*subset_1)[0]=="." ) {
    if(      subset_2 == 0 ) {
      cerr << "Unable to find a complementary set to an empty set." << endl;
      return false;
    }
    else if( (*subset_2)[0] == "." ) {
      cerr << "Unable to estimate interaction between . and ." << endl;
      return false;
    }
    else
      swap(subset_1,subset_2);
  }
  if( subset_1 != 0 ) {
    if( subset_2!=0 && (*subset_2)[0]=="." ) {
      for( int d=0;d<testVars.size();d++ ) {
	if( find(subset_1->begin(),subset_1->end(),testVars[d])
	    == subset_1->end() ) {
	  complementary.push_back(testVars[d]);
	}
      }
    }
  }
  if( !complementary.empty() )
    subset_2 = &complementary;

  // make sure all requested vars are present
  if( subset_1 != 0 ) {
    for( int d=0;d<subset_1->size();d++ ) {
      if( find(testVars.begin(),testVars.end(),(*subset_1)[d])
	  == testVars.end() ) {
	cerr << "Variable " << (*subset_1)[d].c_str() 
	     << " not found among classifier variables." << endl;
	return false;
      }
    }
  }
  if( subset_2!=0 && complementary.empty() ) {
    for( int d=0;d<subset_2->size();d++ ) {
      if( find(testVars.begin(),testVars.end(),(*subset_2)[d])
	  == testVars.end() ) {
	cerr << "Variable " << (*subset_2)[d].c_str() 
	     << " not found among classifier variables." << endl;
	return false;
      }
    }
  }

  // if both subsets are filled, estimate specified interaction
  if( subset_2 != 0 ) {
    string interactionName;
    for( int d=0;d<subset_1->size()-1;d++ ) {
      interactionName += (*subset_1)[d];
      interactionName += ",";
    }
    interactionName += (*subset_1)[subset_1->size()-1];
    interactionName += ":";
    for( int d=0;d<subset_2->size()-1;d++ ) {
      interactionName += (*subset_2)[d];
      interactionName += ",";
    }
    interactionName += (*subset_2)[subset_2->size()-1];
    if( verbose > 0 )
      cout << "Estimating interaction     " << interactionName.c_str() << endl;
    interaction.resize(1);
    SprValueWithError value;
    if( !SprClassifierEvaluator::twoSubsetInteraction(data,
						      trained,
						      mcTrained,
						      mapper,
						      *subset_1,
						      *subset_2,
						      nPoints,
						      value,
						      verbose) ) {
      cerr << "Unable to compute interaction." << endl;
      return false;
    }
    interaction[0] = NameAndValue(interactionName,value);
  }

  // if neither subset is filled, estimate interaction between each variable
  // and the rest
  if( subset_1 == 0 ) {
    int dim = testVars.size();
    interaction.resize(dim);
    SprValueWithError value;
    for( int d=0;d<dim;d++ ) {
      string interactionName = testVars[d];
      interactionName += ":.";
      if( verbose > 0 ) {
	cout << "Estimating interaction     " 
	     << interactionName.c_str() << endl;
      }
      vector<string> one_var(1,testVars[d]);
      vector<string> other_vars;
      for( int k=0;k<dim;k++ ) {
	if( k == d ) continue;
	other_vars.push_back(testVars[k]);
      }
      if( !SprClassifierEvaluator::twoSubsetInteraction(data,
							trained,
							mcTrained,
							mapper,
							one_var,
							other_vars,
							nPoints,
							value,
							verbose) ) {
	cerr << "Unable to compute interaction." << endl;
	return false;
      }
      interaction[d] = NameAndValue(interactionName,value);
    }
  }

  // if one subset is filled and another one is not, 
  // compute interaction between first subset and each variable
  // not included in this subset
  if( subset_2==0 && subset_1!=0 ) {
    int dim = testVars.size();
    int nCompute = dim - subset_1->size();
    assert( nCompute >= 0 );
    interaction.resize(nCompute);
    string baseInteractionName;
    for( int d=0;d<subset_1->size()-1;d++ ) {
      baseInteractionName += (*subset_1)[d];
      baseInteractionName += ",";
    }
    baseInteractionName += (*subset_1)[subset_1->size()-1];
    baseInteractionName += ":";
    SprValueWithError value;
    int index = 0;
    for( int d=0;d<dim;d++ ) {
      if( find(subset_1->begin(),subset_1->end(),testVars[d])
	  != subset_1->end() ) continue;
      string interactionName = baseInteractionName;
      interactionName += testVars[d];
      if( verbose > 0 ) {
	cout << "Estimating interaction     " 
	     << interactionName.c_str() << endl;
      }
      vector<string> one_var(1,testVars[d]);
      if( !SprClassifierEvaluator::twoSubsetInteraction(data,
							trained,
							mcTrained,
							mapper,
							*subset_1,
							one_var,
							nPoints,
							value,
							verbose) ) {
	cerr << "Unable to compute interaction." << endl;
	return false;
      }
      interaction[index++] = NameAndValue(interactionName,value);	
    }
  }

  // exit
  return true;
}


bool SprClassifierEvaluator::addNremoveR(
		   const SprAbsFilter* trainData,
		   const SprAbsFilter* testData,
		   SprAbsClassifier* trainable,
		   SprAbsMultiClassLearner* mcTrainable,
		   unsigned N, unsigned R,
		   SprAverageLoss* aveLoss,
		   unsigned nCross,
		   std::vector<SetVarsAndValue>& vars_and_loss,
		   bool integrate,
		   int verbose)
{
  // sanity check
  if( trainData == 0 ) {
    cerr << "Training data not supplied for addNremoveR." << endl;
    return false;
  }
  if( testData==0 && nCross==0 ) {
    cerr << "Must supply either test data or the number of cross-validation "
	 << "pieces for addNremoveR." << endl;
    return false;
  }
  if( N == 0 ) {
    cerr << "Must specify a positive number of variables for addition." 
	 << endl;
    return false;
  }
  if( trainable!=0 && mcTrainable!=0 ) {
    cerr << "addNremoveR cannot process both binary classifier and " 
	 << "multiclass learner."<< endl;
    return false;
  }
  if( trainable==0 && mcTrainable==0 ) {
    cerr << "No classifier supplied for addNremoveR." << endl;
    return false;
  }
  if( aveLoss == 0 ) {
    cerr << "No loss specified for addNremoveR." << endl;
    return false;
  }

  // check classes
  vector<SprClass> classes; 
  trainData->classes(classes); 
  if( classes.size() < 2 ) {
    cerr << "Classes have not been set." << endl;
    return false; 
  }

  // get a list of full vars
  vector<string> allVars;
  trainData->vars(allVars);
  assert( !allVars.empty() );
  int dim = allVars.size();

  // init optimal loss
  vars_and_loss.clear();
  vars_and_loss.resize(dim,SetVarsAndValue(set<string>(),
					   SprValueWithError(1,0)));

  // init containers and counters
  int d = 0;// current dimension
  set<string> best_model;// current optimal set of variables

  // print out
  if( verbose > 0 ) {
    cout << setw(40) << "Variable" << " Add(+)/Rem(-) " 
	 << setw(10) << "Loss" << endl;
  }

  // loop through dimensions
  while( d < dim ) {

    //
    // add
    //
    for( int n=0;n<N;n++ ) {
      // check dimension index
      if( d >= dim ) break;

      // make a vector of losses
      vector<NameAndValue> loss_per_var;
      
      // loop over vars
      for( int v=0;v<allVars.size();v++ ) {
	
	// if the var is not included yet, try it
	if( best_model.find(allVars[v]) == best_model.end() ) {

	  // add to the attempted model
	  set<string> attempt_model(best_model);
	  attempt_model.insert(allVars[v]);

	  // get the best variable and associated loss
	  pair<SprValueWithError,bool> lossAndStatus 
	    = SprClassifierEvaluator::lossPerVar(attempt_model,
						 trainData,testData,nCross,
						 trainable,mcTrainable,
						 aveLoss,classes,
						 integrate,verbose-1);
	  if( !lossAndStatus.second ) {
	    cerr << "Unable to compute loss at add step for variable " 
		 << allVars[v] << endl;
	    return false;
	  }

	  // add loss to the list
	  loss_per_var.push_back(NameAndValue(allVars[v],lossAndStatus.first));
	}
      }// end of loop over vars

      // get the best loss
      NameAndValue bestLossPair 
	= *min_element(loss_per_var.begin(),loss_per_var.end(),
		       SCECmpPairSDDSecondFirst());

      // show losses
      if( verbose > 1 ) {
	cout << "Losses:   ";
	for( int i=0;i<loss_per_var.size();i++ ) {
	  cout << loss_per_var[i].first << ": " 
	       << loss_per_var[i].second.first << " ";
	}
	cout << endl;
      }

      // update the best model
      if( verbose > 0 ) {
	cout << setw(40) << bestLossPair.first << "       +       "
	     << setw(10) << bestLossPair.second.first << endl;
      }
      best_model.insert(bestLossPair.first);      
      vars_and_loss[d++] = SetVarsAndValue(best_model,bestLossPair.second);
    }// end of add

    // check dimensionality
    if( d >= dim ) break;

    //
    // remove
    //
    for( int r=0;r<R;r++ ) {
      // check dimension
      if( d < 1 ) break;

      // make a vector of losses
      vector<NameAndValue> loss_per_var;
      
      // loop over vars
      for( set<string>::const_iterator 
	     v=best_model.begin();v!=best_model.end();v++ ) {

	// remove from attempted model
	set<string> attempt_model(best_model);
	attempt_model.erase(attempt_model.find(*v));
	
	// get the best variable and associated loss
	pair<SprValueWithError,bool> lossAndStatus 
	  = SprClassifierEvaluator::lossPerVar(attempt_model,
					       trainData,testData,nCross,
					       trainable,mcTrainable,
					       aveLoss,classes,
					       integrate,verbose-1);
	if( !lossAndStatus.second ) {
	  cerr << "Unable to compute loss at remove step for variable " 
	       << *v << endl;
	  return false;
	}

	// add loss to the list
	loss_per_var.push_back(NameAndValue(*v,lossAndStatus.first));
      }// end loop over vars

      // get the best loss
      NameAndValue bestLossPair 
	= *min_element(loss_per_var.begin(),loss_per_var.end(),
		       SCECmpPairSDDSecondFirst());

      // show losses
      if( verbose > 1 ) {
	cout << "Losses:   ";
	for( int i=0;i<loss_per_var.size();i++ ) {
	  cout << loss_per_var[i].first << ": " 
	       << loss_per_var[i].second.first << " ";
	}
	cout << endl;
      }

      // update the best model
      if( verbose > 0 ) {
	cout << setw(40) << bestLossPair.first << "       -       " 
	     << setw(10) << bestLossPair.second.first << endl;
      }
      best_model.erase(best_model.find(bestLossPair.first));
      vars_and_loss[--d] 
	= SetVarsAndValue(set<string>(),SprValueWithError(1,0));

      // check dimension
      if( d < 1 ) break;
      vars_and_loss[d-1] = SetVarsAndValue(best_model,bestLossPair.second); 
    }// end of remove loop
  }// end of while d < dim

  // exit
  return true;
}
      

std::pair<SprValueWithError,bool> SprClassifierEvaluator::lossPerVar(
			     const std::set<std::string>& attempt_model,
			     const SprAbsFilter* trainData,
			     const SprAbsFilter* testData,
			     unsigned nCross,
			     SprAbsClassifier* trainable,
			     SprAbsMultiClassLearner* mcTrainable,
			     SprAverageLoss* aveLoss,
			     const std::vector<SprClass>& classes,
			     bool integrate,
			     int verbose)
{
  // sanity check
  assert( trainable!=0 || mcTrainable!=0 );
  assert( trainable==0 || mcTrainable==0 );

  // init
  pair<SprValueWithError,bool> 
    lossAndStatus(make_pair(SprUtils::max(),0.),false);

  // make a selector
  SprVarSelectorFilter train_selector(trainData);
      
  // make test data
  auto_ptr<SprVarSelectorFilter> test_selector; 
  if( testData != 0 )
    test_selector.reset(new SprVarSelectorFilter(testData));
      
  // select vars
  if( !train_selector.chooseVars(attempt_model) 
      || (test_selector.get()!=0 &&
	  !test_selector->chooseVars(attempt_model)) ) {
    cerr << "Unable to choose variables." << endl;
    return lossAndStatus;
  }
  if( !train_selector.filter() ) {
    cerr << "Unable to filter training variable selector." << endl;
    return lossAndStatus;
  }
  if( test_selector.get()!=0 && !test_selector->filter() ) {
    cerr << "Unable to filter test variable selector." << endl;
    return lossAndStatus;
  }

  // 
  // test data supplied by the user
  //
  if( test_selector.get() != 0 ) {
    // set data
    if( (trainable!=0 && !trainable->setData(&train_selector)) 
	|| (mcTrainable!=0 && !mcTrainable->setData(&train_selector)) ) {
      cerr << "Unable to reset data for classifier in addNremoveR." << endl;
      return lossAndStatus;
    }
    
    // train
    if( (trainable!=0 && !trainable->train(verbose))
	|| (mcTrainable!=0 && !mcTrainable->train(verbose)) ) {
      cerr << "Unable to train classifier in addNremoveR." << endl;
      return lossAndStatus;
    }
    auto_ptr<SprAbsTrainedClassifier> trained;
    auto_ptr<SprAbsTrainedMultiClassLearner> mcTrained;
    if( trainable != 0 )
      trained.reset(trainable->makeTrained());
    if( mcTrainable != 0 )
      mcTrained.reset(mcTrainable->makeTrained());
    if( trained.get()==0 && mcTrained.get()==0 ) {
      cerr << "Unable to make trained classifier in addNremoveR." << endl;
      return lossAndStatus;
    }
    
    // switch output range
    if( trained.get() != 0 ) trained->useNormalized();
    
    // evaluate classifier performance
    vector<const SprAbsFilter*> fom_data(1,test_selector.get());
    if(      trained.get() != 0 ) {
      vector<const SprAbsTrainedClassifier*> fom_trained(1,trained.get());
      lossAndStatus.first = SprFomCalculator::fom(fom_data,
						  fom_trained,
						  0,aveLoss,
						  classes[0],classes[1],
						  false,verbose-1);
    }
    else if( mcTrained.get() != 0 ) {
      vector<const SprAbsTrainedMultiClassLearner*> 
	fom_trained(1,mcTrained.get());
      SprClassificationTable classificationTable;
      map<int,double> weightInClass;
      lossAndStatus.first = SprFomCalculator::loss(fom_data,
						   fom_trained,
						   aveLoss,
						   classificationTable,
						   weightInClass,
						   verbose-1);
    }
    lossAndStatus.second = true;
  }
  //
  // test data not supplied by the user
  //
  else {
    // cross validate
    SprCrossValidator cv(&train_selector,nCross);
    if(      trainable != 0 ) {
      vector<SprAbsClassifier*> classifiers(1,trainable);
      vector<SprValueWithError> cvFOM;
      if( !cv.validate(0,aveLoss,classifiers,classes[0],classes[1],
		       cvFOM,integrate,verbose-1) ) {
	cerr << "Unable to cross-validate for addNremoveR." << endl;
	return lossAndStatus;
      }
      assert( cvFOM.size() == 1 );

      // assign loss
      lossAndStatus.first = cvFOM[0];
      lossAndStatus.second = true;
    }// end binary classifier
    else if( mcTrainable != 0 ) {
      SprValueWithError realLoss;
      SprClassificationTable classificationTable;
      map<int,double> weightInClass;
      if( !cv.validate(mcTrainable,aveLoss,
		       realLoss,classificationTable,
		       weightInClass,verbose-1) ) {
	cerr << "Unable to cross-validate for addNremoveR." << endl;
	return lossAndStatus;
      }
      
      // assign loss
      lossAndStatus.first = realLoss;
      lossAndStatus.second = true;
    }// end multiclass learner
  }

  // exit
  return lossAndStatus;
}


bool SprClassifierEvaluator::sortByInteraction(
				const SprAbsFilter* data,
				SprAbsTrainedClassifier* trained,
				SprAbsTrainedMultiClassLearner* mcTrained,
				SprCoordinateMapper* mapper,
				unsigned nPoints,
				std::vector<ListVarsAndValue>& interaction,
				int verbose)
{
  // sanity check
  if( trained==0 && mcTrained==0 ) {
    cerr << "No classifiers specified for sortByInteraction." << endl;
    return false;
  }

  // get trained vars
  vector<string> testVars;
  if(      trained != 0 )
    trained->vars(testVars);
  else if( mcTrained != 0 )
    mcTrained->vars(testVars);
  assert( !testVars.empty() );
  int dim = testVars.size();

  // init
  interaction.clear();

  // init vector of best vars
  vector<string> bestVars, notBestVars(testVars);

  // loop through variables
  while( bestVars.size() < (dim-1) ) {
    // make a vector of vars to attempt at this step
    vector<NameAndValue> attempted;

    // add vars one by one and estimate interaction
    for( int d=0;d<dim;d++ ) {
      if( find(bestVars.begin(),bestVars.end(),testVars[d])
	  != bestVars.end() ) continue;
      if( verbose > 1 )
	cout << "Attempting variable   " << testVars[d].c_str() << endl;
      //      vector<string> attemptVars = bestVars;
      vector<string> attemptVars(1,testVars[d]);
      vector<string> notAttemptVars = notBestVars;
      //      attemptVars.push_back(testVars[d]);
      notAttemptVars.erase(remove(notAttemptVars.begin(),
				  notAttemptVars.end(),testVars[d]),
			   notAttemptVars.end());
      SprValueWithError value;
      if( !SprClassifierEvaluator::twoSubsetInteraction(data,
							trained,
							mcTrained,
							mapper,
							attemptVars,
							notAttemptVars,
							nPoints,
							value,
							verbose) ) {
	cerr << "Unable to compute interaction." << endl;
	return false;
      }
      attempted.push_back(NameAndValue(testVars[d],value));
    }// end loop over dimensions

    // find max interaction
    int indMax = -1;
    double valMax = SprUtils::min();
    for( int i=0;i<attempted.size();i++ ) {
      if( attempted[i].second.first > valMax ) {
	valMax = attempted[i].second.first;
	indMax = i;
      }
    }
    assert( indMax >= 0 );

    // record max interaction
    bestVars.push_back(attempted[indMax].first);
    notBestVars.erase(remove(notBestVars.begin(),
			     notBestVars.end(),attempted[indMax].first),
		      notBestVars.end());
    interaction.push_back(ListVarsAndValue(bestVars,attempted[indMax].second));
    if( verbose > 0 ) {
      char s [200];
      sprintf(s,"Adding variable  %35s  with interaction   %15.10f +- %15.10f",
	      attempted[indMax].first.c_str(),
	      attempted[indMax].second.first,attempted[indMax].second.second);
      cout << s << endl;
    }
  }

  // exit
  return true;
}
