#!/bin/bash
#
# find rootcint
#
which rootcint
if [ $? != 0 ]; then
  echo "rootcint not found. Exiting."
  exit 1
fi
#
# run rootcint
#
rootcint SprAdapterDict.C -c -p -v -I./include include/StatPatternRecognition/SprRootAdapter.hh
#
# fix include paths
#
awk '{ gsub("include/StatPatternRecognition/SprRootAdapter.hh","StatPatternRecognition/SprRootAdapter.hh"); print $0 }' SprAdapterDict.C > SprAdapterDict.C.new
mv -f SprAdapterDict.C.new SprAdapterDict.C
awk '{ gsub("SprAdapterDict.h","StatPatternRecognition/SprAdapterDict.h"); print $0 }' SprAdapterDict.C > SprAdapterDict.C.new
mv -f SprAdapterDict.C.new SprAdapterDict.C
#
awk '{ gsub("include/StatPatternRecognition/SprRootAdapter.hh","StatPatternRecognition/SprRootAdapter.hh"); print $0 }' SprAdapterDict.h > SprAdapterDict.h.new
mv -f SprAdapterDict.h.new SprAdapterDict.h
awk '{ gsub("SprAdapterDict.h","StatPatternRecognition/SprAdapterDict.h"); print $0 }' SprAdapterDict.h > SprAdapterDict.h.new
mv -f SprAdapterDict.h.new SprAdapterDict.h
#
# move files to proper dirs
#
mv -f SprAdapterDict.h include/StatPatternRecognition/
mv -f SprAdapterDict.C src/
#
# exit
#
exit 0