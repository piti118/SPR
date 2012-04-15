#!/bin/bash
for i in *.cc; do awk '{ gsub("StatPatternRecognition/Spr","PhysicsTools/StatPatternRecognition/interface/Spr"); print $0 }' $i > $i.new; mv -f $i.new $i; done
for i in *.cc; do awk '{ gsub("math/Spr","PhysicsTools/StatPatternRecognition/src/Spr"); print $0 }' $i > $i.new; mv -f $i.new $i; done
