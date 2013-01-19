// File and Version Information:
//      $Id: SprMultiClassReader.hh,v 1.5 2008-05-08 19:57:43 narsky Exp $
//
// Description:
//      Class SprMultiClassReader :
//          Reads trained multi class learner from a file.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Ilya Narsky                     Original author
//
// Copyright Information:
//      Copyright (C) 2005, 2006        California Institute of Technology
//
//------------------------------------------------------------------------
 
#ifndef _SprMultiClassReader_HH
#define _SprMultiClassReader_HH

#include "StatPatternRecognition/SprAbsTrainedMultiClassLearner.hh"
#include "StatPatternRecognition/SprMatrix.hh"

#include <vector>
#include <string>
#include <istream>

class SprAbsMultiClassLearner;
class SprTrainedMultiClassLearner;
class SprTrainedBinaryEncoder;


class SprMultiClassReader
{
public:
  virtual ~SprMultiClassReader();

  SprMultiClassReader() {}

  // Read trained classifier from a file or from an input stream.
  static SprAbsTrainedMultiClassLearner* readTrained(const char* filename,
						     int verbose=0);
  static SprAbsTrainedMultiClassLearner* readTrained(std::istream& input,
						     int verbose=0);

  // Read indicator matrix.
  static bool readIndicatorMatrix(const char* filename, SprMatrix& indicator);

  /*
    Prepare indicator matrix based on supplied class info
    for individual classifiers.
  */
  static bool makeIndicatorMatrixByClass(
				const std::vector<int>& mapper,
				const std::vector<std::string>& classGroup,
				SprMatrix& indicator);

  // Construct a trained Allwein-Schapire-Singer learner from a file
  // with class groups and binary classifier file names.
  static SprTrainedMultiClassLearner* readBinaryList(
					 const char* filename,
					 int verbose=0);
  static SprTrainedMultiClassLearner* readBinaryList(
					 std::istream& input,
					 int verbose=0);

private:
  /*
    Reads trained classifier from random position in a stream.
  */
  static SprAbsTrainedMultiClassLearner* readTrainedFromStream(
					       std::istream& input,
                                               const std::string& classifier,
					       unsigned& lineCounter);

  /*
    Reads trained classifier of known type from random position in a stream.
  */
  static SprAbsTrainedMultiClassLearner* readSelectedTrained(
					       std::istream& input,
                                               const std::string& classifier,
					       unsigned& lineCounter);

  // Read missing values.
  static bool readMissing(
    std::istream& input,
    SprCut& validRange,
    std::vector<SprAbsTrainedMultiClassLearner::MissingType>& defaultMissing, 
    unsigned& nLine);

  // specific reader implementations
  static SprTrainedBinaryEncoder* readBinaryEncoder(
					    std::istream& input,
					    unsigned& lineCounter);
  static SprTrainedMultiClassLearner* readMultiClassLearner(
						std::istream& input,
                                                unsigned& lineCounter);
};

#endif
