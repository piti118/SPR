//$Id: SprFomCalculator.cc,v 1.8 2008-04-02 23:36:45 narsky Exp $

#include "StatPatternRecognition/SprExperiment.hh"
#include "StatPatternRecognition/SprFomCalculator.hh"
#include "StatPatternRecognition/SprAbsFilter.hh"
#include "StatPatternRecognition/SprAbsTwoClassCriterion.hh"
#include "StatPatternRecognition/SprAverageLoss.hh"
#include "StatPatternRecognition/SprLoss.hh"
#include "StatPatternRecognition/SprAbsTrainedClassifier.hh"
#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprUtils.hh"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <algorithm>
#include <sstream>

using namespace std;


SprValueWithError SprFomCalculator::fom(
          const std::vector<std::vector<SprPlotter::Response> >& responses,
	  const SprAbsTwoClassCriterion* crit, 
	  SprAverageLoss* loss, 
	  bool integrate, int verbose)
{
  // sanity check
  assert( !responses.empty() );
  assert( crit!=0 || loss!=0 );

  // check if there are filled out responses
  int N = responses.size();
  bool filledResponse = false;
  for( int n=0;n<N;n++ ) {
    if( !responses[n].empty() ) {
      filledResponse = true;
      break;
    }
  }
  if( !filledResponse ) {
    cerr << "No filled responses found in SprFomCalculator::fom." << endl;
    return SprValueWithError(SprUtils::min(),0);
  }

  // apply classifiers to samples
  double wcor0(0), wcor1(0), wmis0(0), wmis1(0);
  vector<double> v_fom(N);
  for( int n=0;n<N;n++ ) {
    // if no responses for classifier n, skip
    if( responses[n].empty() ) continue;
 
    // reset if not integrating
    if( !integrate ) {
      if( loss != 0 ) loss->reset();
      wcor0 = 0;
      wcor1 = 0;
      wmis0 = 0;
      wmis1 = 0;
    }

    // loop thru data
    const vector<SprPlotter::Response>& response = responses[n];
    for( int i=0;i<response.size();i++ ) {
      int cls = response[i].cls;
      double w = response[i].weight;
      assert( !response[i].accepted.empty() );
      int accepted = response[i].accepted.begin()->second;
      if( loss == 0 ) {
	if( accepted ) {
	  if(      cls == 0 )
	    wmis0 += w;
	  else if( cls == 1 )
	    wcor1 += w;
	}
	else {
	  if(      cls == 0 )
	    wcor0 += w;
	  else if( cls == 1 )
	    wmis1 += w;
	}
      }
      else {
	assert( !response[i].response.empty() );
	double r = response[i].response.begin()->second;
	loss->update(cls,r,w);
      }
    }// end loop through data

    // compute fom
    if( !integrate ) {
      if( loss == 0 ) 
	v_fom[n] = crit->fom(wcor0,wmis0,wcor1,wmis1);
      else
	v_fom[n] = loss->value();
    }
  }// end loop over subsamples

  // if integrating over subsamples, compute overall fom
  if( integrate ) {
    if( verbose > 0 ) {
      cout << "Computed integrated weights:   "
	   << "    wcor0=" << wcor0
	   << "    wmis0=" << wmis0
	   << "    wcor1=" << wcor1
	   << "    wmis1=" << wmis1
	   << endl;
    }
    double integrated_fom = 0;
    if( loss == 0 ) 
      integrated_fom = crit->fom(wcor0,wmis0,wcor1,wmis1);
    else
      integrated_fom = loss->value();
    return SprValueWithError(integrated_fom,0);
  }

  // if not integrating, estimate error
  if( verbose > 0 ) {
    cout << "Computed FOMs for " << v_fom.size() << " subsamples:    ";
    for( int i=0;i<v_fom.size();i++ )
      cout << v_fom[i] << " ";
    cout << endl;
  }
  double ave = 0;
  double err = 0;
  for( int i=0;i<v_fom.size();i++ )
    ave += v_fom[i];
  ave /= v_fom.size();
  for( int i=0;i<v_fom.size();i++ )
    err += (v_fom[i]-ave)*(v_fom[i]-ave);
  if( v_fom.size() > 1 ) err /= (v_fom.size()-1);
  err = ( err>0 ? sqrt(err) : 0 );
  return SprValueWithError(ave,err);
}


SprValueWithError SprFomCalculator::fom(
		 const std::vector<const SprAbsFilter*>& data, 
		 const std::vector<const SprAbsTrainedClassifier*>& trained,
		 const SprAbsTwoClassCriterion* crit, 
		 SprAverageLoss* loss,
		 const SprClass& cls0, const SprClass& cls1,
		 bool integrate, int verbose)
{
  // sanity check
  assert( !data.empty() );
  assert( !trained.empty() );
  assert( data.size() == trained.size() );
  assert( cls0 != cls1 );
  assert( crit!=0 || loss!=0 );

  // check if there are non-null classifiers
  int N = data.size();
  bool nonNullTrained = false;
  for( int n=0;n<N;n++ ) {
    if( trained[n] != 0 ) {
      nonNullTrained = true;
      break;
    }
  }
  if( !nonNullTrained ) {
    cerr << "No trained classifiers found in SprFomCalculator::fom." << endl;
    return SprValueWithError(SprUtils::min(),0);
  }

  // fill out responses
  vector<vector<SprPlotter::Response> > responses(N);
  for( int n=0;n<N;n++ ) {
    // check classifier and data
    if( trained[n]==0 || data[n]==0 || data[n]->empty() ) continue;

    // loop thru subsample and fill out trained responses
    const SprAbsFilter* sub = data[n];
    for( int i=0;i<sub->size();i++ ) {
      const SprPoint* p = (*sub)[i];
      int cls = -1;
      if(      p->class_ == cls0 )
	cls = 0;
      else if( p->class_ == cls1 )
	cls = 1;
      else
	continue;
      double w = sub->w(i);
      SprPlotter::Response resp(cls,w);
      double r = 0;
      int accept = ( trained[n]->accept(p,r) ? 1 : 0 );
      ostringstream ost;
      ost << i;
      resp.set(ost.str().c_str(),r,accept);
      responses[n].push_back(resp);
    }// end loop over points in this subsample
  }// end loop over subsets

  // exit
  return SprFomCalculator::fom(responses,crit,loss,integrate,verbose-1);
}


SprValueWithError SprFomCalculator::loss(
		 const std::vector<int>& classes,
                 const std::vector<SprMultiClassPlotter::Response>& responses,
		 SprAverageLoss* loss,
		 SprClassificationTable& classificationTable, 
		 std::map<int,double>& weightInClass,
		 int verbose)
{
  // sanity check
  assert( !responses.empty() );

  // get the loss table
  SprMultiClassPlotter plotter(responses);
  double totalLoss = plotter.multiClassTable(classes,
					     loss,
					     classificationTable,
					     weightInClass);

  // exit
  return SprValueWithError(totalLoss,0);
}


SprValueWithError SprFomCalculator::loss(
	   const std::vector<const SprAbsFilter*>& data, 
	   const std::vector<const SprAbsTrainedMultiClassLearner*>& trained,
	   SprAverageLoss* loss,
	   SprClassificationTable& classificationTable, 
	   std::map<int,double>& weightInClass,
	   int verbose)
{
  // sanity check
  assert( !data.empty() );
  assert( !trained.empty() );
  assert( data.size() == trained.size() );

  // check if there are non-null classifiers and get classes
  vector<int> classes;
  int N = data.size();
  bool nonNullTrained = false;
  for( int n=0;n<N;n++ ) {
    if( trained[n] != 0 ) {
      trained[n]->classes(classes);
      nonNullTrained = true;
      break;
    }
  }
  if( !nonNullTrained ) {
    cerr << "No trained classifiers found in SprFomCalculator::loss." << endl;
    return SprValueWithError(SprUtils::min(),0);
  }
  assert( !classes.empty() );

  // apply classifiers to samples
  vector<SprMultiClassPlotter::Response> responses;
  for( int n=0;n<N;n++ ) {
    // if trained classifier is 0, skip
    if( trained[n] == 0 ) continue;

    // check number of classes
    // a serious atempt to make sure that classes match for all classifiers
    // is not made to save cpu
    vector<int> local_classes;
    trained[n]->classes(local_classes);
    assert( local_classes.size() == classes.size() );
 
    // loop thru data
    for( int i=0;i<data[n]->size();i++ ) {
      const SprPoint* p = (*data[n])[i];
      int cls = p->class_;
      if( find(classes.begin(),classes.end(),cls) == classes.end() )
	continue;
      double w = data[n]->w(i);
      map<int,double> output;
      int r = trained[n]->response(p,output);
      responses.push_back(SprMultiClassPlotter::Response(cls,w,r,output));
    }// end loop through data
  }// end loop over subsamples
  assert( !responses.empty() );

  // get the loss table
  return SprFomCalculator::loss(classes,
				responses,
				loss,
				classificationTable,
				weightInClass,
				verbose);
}
