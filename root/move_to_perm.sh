#!/bin/bash
cp -f /raid1/narsky/SPRroot/share/lambda-{train,test}.pat /home/narsky/StatPatternRecognition
cp -f /raid1/narsky/SPRroot/lib/libSPR.so /home/narsky/StatPatternRecognition
awk '{ gsub("/raid1/narsky/SPRroot/lib/",""); gsub("/raid1/narsky/SPRroot/share/",""); print $0 }' root/spr_tutorial.C > spr_tutorial.C.new
mv -f spr_tutorial.C.new /home/narsky/StatPatternRecognition/spr_tutorial.C
awk '{ gsub("/raid1/narsky/SPRroot/lib/",""); gsub("/raid1/narsky/SPRroot/share/",""); print $0 }' root/spr_multiclass.C > spr_multiclass.C.new
mv -f spr_multiclass.C.new /home/narsky/StatPatternRecognition/spr_multiclass.C
cp -f root/spr_plot.C /home/narsky/StatPatternRecognition/spr_plot.C
cp -f include/StatPatternRecognition/SprRootAdapter.hh /home/narsky/StatPatternRecognition
exit 0
