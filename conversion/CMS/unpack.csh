#!/bin/tcsh -fvx

gunzip *.gz
tar -xvf *tar
set CTSPR=CommonTools/StatPatternRecognition
mkdir CommonTools
mkdir ${CTSPR}
mkdir ${CTSPR}/interface/
mkdir ${CTSPR}/src
mkdir ${CTSPR}/test
mkdir ${CTSPR}/bin

# move everything from include/StatPatternRecognition to interface
foreach i (StatPatternRecognition/src/*cc)
  ./jigger.csh $i >! ${CTSPR}/src/$i:t
end

foreach i (${CTSPR}/src/example*cc)
  mv $i ${CTSPR}/test/
end

foreach i (${CTSPR}/src/*App.cc)
  mv $i ${CTSPR}/bin/
end

foreach i (StatPatternRecognition/include/StatPatternRecognition/*hh)
  ./jigger.csh $i >! ${CTSPR}/interface/$i:t
end

foreach i (StatPatternRecognition/math/*hh StatPatternRecognition/math/*cc)
  ./jigger.csh $i >! ${CTSPR}/src/$i:t
end
