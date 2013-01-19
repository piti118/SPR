lib_LTLIBRARIES = libSPR.la

AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)
AM_LDFLAGS  = -L/usr/lib -L$(top_srcdir)/src
AM_LIBS     = -lSPR 
APPLDADD    = -lSPR
if WITH_ROOT
   AM_CXXFLAGS += $(ROOTCXXFLAGS) -O2
   AM_CPPFLAGS += $(ROOTCPPFLAGS) -DSPRROOTTUPLE
   AM_LDFLAGS  += $(ROOTLDFLAGS)
   AM_LIBS     += $(ROOTLIBS)
   APPLDADD    += $(ROOTLIBS)
else
   AM_CXXFLAGS += -O4
endif

libSPR_la_SOURCES  = \
  SprDataMoments.cc                 SprDecisionTree.cc \
  SprFisher.cc                      SprFomCalculator.cc \
  SprAbsClassifier.cc               SprCombiner.cc \
  SprIntegerBootstrap.cc            SprAbsFilter.cc \
  SprIntegerPermutator.cc           SprAbsTrainedClassifier.cc \
  SprTrainedCombiner.cc             SprMultiClassReader.cc \
  SprMultiClassLearner.cc           SprRandomNumber.cc \
  SprAdaBoost.cc                    SprTrainedAdaBoost.cc \
  SprSimpleReader.cc                SprArcE4.cc \
  SprAsciiWriter.cc                 SprStringParser.cc \
  SprBagger.cc                      SprTopdownTree.cc \
  SprTrainedBagger.cc               SprClassifierReader.cc \
  SprBinarySplit.cc                 SprTrainedBinarySplit.cc \
  SprBootstrap.cc                   SprTrainedDecisionTree.cc \
  SprBoxFilter.cc                   SprTrainedFisher.cc \
  SprTrainedMultiClassLearner.cc    SprBumpHunter.cc \
  SprTrainedRBF.cc                  SprClass.cc \
  SprTrainedStdBackprop.cc          SprCrossValidator.cc \
  SprTrainedTopdownTree.cc          SprData.cc \
  SprTreeNode.cc                    SprDataFeeder.cc \
  SprLogitR.cc                      SprTrainedLogitR.cc \
  SprStdBackprop.cc                 SprTransformation.cc \
  SprAbsTrainedMultiClassLearner.cc SprAbsMultiClassLearner.cc \
  SprPreFilter.cc                   SprCoordinateMapper.cc \
  SprPlotter.cc                     SprMultiClassPlotter.cc \
  SprClassifierEvaluator.cc         SprAbsVarTransformer.cc \
  SprPCATransformer.cc              SprTransformerFilter.cc \
  SprVarTransformerReader.cc        SprVarSelectorFilter.cc \
  SprVarTransformerSequence.cc      SprInputNormalizer.cc \
  SprRangeBooster.cc                SprTrainedRangeBooster.cc \
  SprGEP.cc                         SprTrainedGEP.cc \
  SprChromosome.cc                  SprTrainedBinaryEncoder.cc \
  SprBinaryEncoder.cc               SprReplaceMissing.cc \
  SprIndicatorMatrix.cc             SprConfig.cc \
  SprChiCdf.cc                      SprGenMatrix.cc \
  SprMatrix.cc                      SprRanluxEngine.cc \
  SprSymMatrix.cc                   SprSymMatrixInvert.cc \
  SprVector.cc

if WITH_ROOT
   libSPR_la_SOURCES  += SprRootReader.cc SprRootWriter.cc SprAdapterDict.C SprRootAdapter.cc
endif



#####################################################################



noinst_PROGRAMS = exampleBootstrapVarianceEst  exampleCorrTest \
                  exampleUserCuts  exampleMatrixTest

exampleBootstrapVarianceEst_SOURCES = exampleBootstrapVarianceEst.cc
exampleBootstrapVarianceEst_LDADD   = $(APPLDADD)
exampleCorrTest_SOURCES             = exampleCorrTest.cc
exampleCorrTest_LDADD               = $(APPLDADD)
exampleMatrixTest_SOURCES           = exampleMatrixTest.cc
exampleMatrixTest_LDADD             = $(APPLDADD)

bin_PROGRAMS = SprFisherLogitApp             SprDecisionTreeApp         SprExploratoryAnalysisApp \
               SprMultiClassApp              SprGoFDecisionTreeApp      SprAdaBoostBinarySplitApp \
               SprRBFNetApp                  SprAdaBoostDecisionTreeApp SprStdBackpropApp \
               SprCombinerApp                SprBaggerDecisionTreeApp   SprBumpHunterApp  \
               SprInteractiveAnalysisApp     SprOutputAnalyzerApp       SprOutputWriterApp \
               SprBoosterApp                 SprBaggerApp               exampleUserCuts \
               SprVariableImportanceApp      SprIOTestApp               SprAddBaggersApp \
               SprTransformationApp          SprSplitterApp             SprAddNRemoveRApp \
               SprGEPApp                     SprCrossValidatorApp       SprAddColumnsForMCLApp \
               SprIndicatorMatrixApp

SprFisherLogitApp_SOURCES = SprFisherLogitApp.cc
SprFisherLogitApp_LDADD = $(APPLDADD)
SprDecisionTreeApp_SOURCES = SprDecisionTreeApp.cc
SprDecisionTreeApp_LDADD = $(APPLDADD)
SprExploratoryAnalysisApp_SOURCES = SprExploratoryAnalysisApp.cc
SprExploratoryAnalysisApp_LDADD = $(APPLDADD)
SprMultiClassApp_SOURCES = SprMultiClassApp.cc
SprMultiClassApp_LDADD = $(APPLDADD)
SprGoFDecisionTreeApp_SOURCES = SprGoFDecisionTreeApp.cc
SprGoFDecisionTreeApp_LDADD = $(APPLDADD)
SprAdaBoostBinarySplitApp_SOURCES = SprAdaBoostBinarySplitApp.cc
SprAdaBoostBinarySplitApp_LDADD = $(APPLDADD)
SprRBFNetApp_SOURCES = SprRBFNetApp.cc
SprRBFNetApp_LDADD = $(APPLDADD)
SprAdaBoostDecisionTreeApp_SOURCES = SprAdaBoostDecisionTreeApp.cc
SprAdaBoostDecisionTreeApp_LDADD = $(APPLDADD)
SprStdBackpropApp_SOURCES = SprStdBackpropApp.cc
SprStdBackpropApp_LDADD = $(APPLDADD)
SprCombinerApp_SOURCES = SprCombinerApp.cc
SprCombinerApp_LDADD = $(APPLDADD)
SprBaggerDecisionTreeApp_SOURCES = SprBaggerDecisionTreeApp.cc
SprBaggerDecisionTreeApp_LDADD = $(APPLDADD)
SprBumpHunterApp_SOURCES = SprBumpHunterApp.cc
SprBumpHunterApp_LDADD = $(APPLDADD)
SprInteractiveAnalysisApp_SOURCES = SprInteractiveAnalysisApp.cc
SprInteractiveAnalysisApp_LDADD = $(APPLDADD)
SprOutputAnalyzerApp_SOURCES = SprOutputAnalyzerApp.cc
SprOutputAnalyzerApp_LDADD = $(APPLDADD)
SprOutputWriterApp_SOURCES = SprOutputWriterApp.cc
SprOutputWriterApp_LDADD = $(APPLDADD)
SprBoosterApp_SOURCES = SprBoosterApp.cc
SprBoosterApp_LDADD = $(APPLDADD)
SprBaggerApp_SOURCES = SprBaggerApp.cc
SprBaggerApp_LDADD = $(APPLDADD)
exampleUserCuts_SOURCES = exampleUserCuts.cc
exampleUserCuts_LDADD = $(APPLDADD)
SprVariableImportanceApp_SOURCES = SprVariableImportanceApp.cc
SprVariableImportanceApp_LDADD = $(APPLDADD)
SprIOTestApp_SOURCES = SprIOTestApp.cc
SprIOTestApp_LDADD = $(APPLDADD)
SprAddBaggersApp_SOURCES = SprAddBaggersApp.cc
SprAddBaggersApp_LDADD = $(APPLDADD)
SprTransformationApp_SOURCES = SprTransformationApp.cc
SprTransformationApp_LDADD = $(APPLDADD)
SprSplitterApp_SOURCES = SprSplitterApp.cc
SprSplitterApp_LDADD = $(APPLDADD)
SprAddNRemoveRApp_SOURCES = SprAddNRemoveRApp.cc
SprAddNRemoveRApp_LDADD = $(APPLDADD)
SprGEPApp_SOURCES = SprGEPApp.cc
SprGEPApp_LDADD = $(APPLDADD)
SprCrossValidatorApp_SOURCES = SprCrossValidatorApp.cc
SprCrossValidatorApp_LDADD = $(APPLDADD)
SprAddColumnsForMCLApp_SOURCES = SprAddColumnsForMCLApp.cc
SprAddColumnsForMCLApp_LDADD = $(APPLDADD)
SprIndicatorMatrixApp_SOURCES = SprIndicatorMatrixApp.cc
SprIndicatorMatrixApp_LDADD = $(APPLDADD)
