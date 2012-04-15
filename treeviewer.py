#!/usr/bin/env python
import os
import os.path
import math
import sys
import getopt
hasPyLab = False
try: #import matplotlib if exists for drawing piechart
    import matplotlib
    matplotlib.use('Agg')
    import pylab
    hasPyLab = True
except ImportError:
    hasPyLab = False

class Setting:
    #contains all the setting the initial values are default
    #can be overwritten via command line or edit this file
    transparentBackground = False
    verbrose = 0
    numTree = 0 #limit on number of tree
    
    #list of node shape
    #see full list here: http://www.graphviz.org/doc/info/shapes.html
    shapes = ['polygon','triangle','diamond','trapezium','house',
        'pentagon','hexagon','septagon','octagon',
        'invtriangle','invtrapezium','invhouse'
        ]
    signalEdgeColor = '#006837' #edge color for signal
    backgroundEdgeColor ='#a50026' #edge color for background

    #margin for graph sometime graph doesn't show all the label 
    #increase this number to fix it  
    margin=1
    
    #color list for red-->green node
    rgcolor = ["#a50026","#d73027","#f46d43","#fdae61",
        "#fee08b","#ffffbf","#d9ef8b","#a6d96a",
        "#66bd63","#1a9850","#006837"]
        
    #color list for blue yellow node (great for red green color blind)
    bycolor =   ["#543005","#8c510a","#bf812d","#dfc27d",
        "#f6e8c3","#f5f5f5","#c7eae5","#80cdc1",
        "#35978f","#01665e","#003c30"]
    
    #color list for black and white
    bwcolor =   ["#ffffff","#f0f0f0","#d9d9d9",
        "#bdbdbd","#969696","#737373","#525252"]
    
    #current color list for node
    colors =    ["#a50026","#d73027","#f46d43","#fdae61",
        "#fee08b","#ffffbf","#d9ef8b","#a6d96a",
        "#66bd63","#1a9850","#006837"]

    piecolor = [ #rotating color used in drawing piechart add more if you wish
    '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462',
    '#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f',
    '#6699FF','#7F66FF','#CC66FF','#FF66E6','#66E6FF','#2970FF',
    '#004EEB','#FF6699','#66FFCC','#EB9C00','#FFB829','#FF7F66',
    '#66FF7F','#99FF66','#E6FF66','#FFCC66']
    pieOtherColor = '#cccccc'

def help():
    #print usage
    print "Usage %s [-n] sprfile"%sys.argv[0]
    options = {
        '-n':"number: number of tree to parse to dot file (default = 0(all))",
        '-t':':use transparent background instead of white background',
        '-s':':use only diamond for all the variable',
        '-c':':use other color set for tree (default 0: green red 1:blue yellow 2: black and whtie)',
        '-v':'number: verbrosity (1,2,3,4,...)',
        '-h':'print this message',
        '-m':'set the margin of the graph(incase some label are cut off)'
    }
    keys = options.keys()
    keys.sort()
    for option in keys:
        print '\t'+option+' '+options[option]

    
def main():
    opts, args = getopt.getopt(sys.argv[1:],'n:tsv:m:c:bsh:')
    numTree = 0
    global hasPyLab
    success = True
    for opt, arg in opts:
        if opt == '-n':
            Setting.numTree = int(arg)
        if opt == '-t':
            Setting.transparentBackground=True
        if opt == '-s':
            Setting.shapes = ['diamond']
        if opt == '-v':
            Setting.verbrose = int(arg)
        if opt == '-c':
            success &= pickColor(opt)
        if opt == '-m':
            Setting.margin = float(arg)
        if opt == '-h':
            help()
            sys.exit()
    sprfiles = args
    
    if not success or len(sprfiles)==0:
        help()      
        sys.exit(1)

    for sprfile in sprfiles:
        Util.debug('Working On: '+sprfile,-1)
        writeDotFiles(sprfile)
        
def pickColor(opt):
    colorArray = [
        Setting.rgcolor,
        Setting.bycolor,
        Setting.bwcolor
    ]
    errmsg = "Invalid color"
    if not Util.isInteger(opt):
        print errmsg
        return False
    option = int(opt)
    if option >= len(colorArray) or option < 0:
        print errmsg
        return False
    else:
        Setting.colors = colorArray[int(opt)]

def testMultiClass(): #test multiclass
    Setting.verbrose = 10
    sprfile = "multi.spr"
    
    parser = MulticlassParser(sprfile)
    parser.parse()
    basename,ext = os.path.splitext(sprfile)
    dirname=basename+'_tree'
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    
    for multiclassKey,multiclassClassifiers in parser._classifiers.iteritems():
        counter = 0
        for index, classifier in multiclassClassifiers.iteritems():
            outputDotFile = "%s/multi_%d_tree_%d.dot"%(dirname,multiclassKey,index)
            classifier.writeDotFile(outputDotFile)
            counter += 1
            if Setting.numTree != 0 and counter >=Setting.numTree:
                break

def writeDotFiles(sprfile):
    global hasPyLab
    tree = None
    isMulticlass = Util.isMulticlass(sprfile)
    prefix = ''
    if isMulticlass: #check if it is multiclass
        tree = MulticlassParser(sprfile)
    else:
        tree = TreeMaker(sprfile)
    tree.parse()
    #print tree._colname

    basename,ext = os.path.splitext(sprfile)
    dirname=basename+'_tree'
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    classifiersDict = tree._classifiers
    if not isMulticlass:
        #wrap it around another ditionary to mimic the behavior of multiclass
        classifiersDict = { 0: classifiersDict }
    
    for multiIndex, classifiers in classifiersDict.iteritems():
        counter = 0
        if isMulticlass:
            prefix = "multi_%d_"%(multiIndex)
        #for each classifier in multiclass
        #make a dot tree
        for thisid,classifier in classifiers.iteritems():
            outputfile = "%s/%stree_%d.dot"%(dirname,prefix,counter)
            classifier.writeDotFile(outputfile)
            counter+=1
            #check number of tree limit (0 means no limit)
            if Setting.numTree!=0 and counter >= Setting.numTree:
                break
        #now write legend
        if hasPyLab:
            legendWriter = LegendWriter(classifiers)
            weightedDeltaFomFile = "%s/%sweightedDeltaFom.eps"%(dirname,prefix)
            legendWriter.writeWeightedDeltaFomPieChart(weightedDeltaFomFile,0.03)
            splitCountFile = "%s/%ssplitCount.pdf"%(dirname,prefix)
            legendWriter.writeSplitCountPieChart(splitCountFile)
            #legendFile = "%s/%slegend.eps"%(dirname,prefix)
            #legendWriter.writeLegendFile(legendFile)
    del tree
    
class Util:
    
    def normalize(self, value, minSource, maxSource, minTarget, maxTarget):
        """
            linearly transform the range min-maxSource to min-maxTarget
            and compute correspoding value of *value* in the new range
        """
        if maxSource == minSource:
            return (minTarget+maxTarget)/2
        return float(value-minSource)*float(maxTarget-minTarget)/float(maxSource-minSource) + minTarget
    normalize = classmethod(normalize)
    def isInteger(self,value):
        """
            check if the value can be transform into an integer
        """
        try:
            temp = int(value)
            return True
        except ValueError:
            return False
    isInteger = classmethod(isInteger)
    def debug(self, msg, verbrose=0):
        """
            print debug message if setting.verbrose > verbrose
        """
        if Setting.verbrose > verbrose:
            print msg
    debug = classmethod(debug)
    
    def toRGB(self,html):
        """
            transform html color notation to rgb array
        """
        split = (html[1:3], html[3:5], html[5:7])
        return [int(x, 16) for x in split]
    toRGB = classmethod(toRGB)
    
    def colorPSArray(self):
        """
            making index postscript color array(I think)...
        """
        colorTuple = map(Util.toRGB, Setting.colors)
        mc = 255.0
        colorString = map(lambda x: '[ %5.2f %5.2f %5.2f ]'%(x[0]/mc,x[1]/mc,x[2]/mc),colorTuple)
        result = ' '.join(colorString)
        result = '[ '+result+' ]'
        return result
    colorPSArray = classmethod(colorPSArray)
    def isMulticlass(self,sprfile):
        """
            check if this sprfile is a multiclass
        """
        f = open(sprfile)
        s = f.readline()
        a = s.split(' ')
        f.close()
        return len(a) == 3 and a[1]=='MultiClassLearner'
    isMulticlass = classmethod(isMulticlass)

class MulticlassParser:
    #finite state machine for parsing spr file
    _readIndicatorMatrixMode = False
    _currentMulticlassClassifierId = -1
    _currentClassifierId = -1
    _filename  = ''
    _classifierType = ''
    _colreadmode = False
    _classifierReadMode = False
    _classifierCreateMode = False
    
    def __init__(self,filename):
        self._indicatorMatrix = []
        self._classifiers = {}#map of tuple to classifier
        self._filename = filename
        self._colname = {}
    
    def parse(self):
        f = open(self._filename)
        firstline = True
        
        for sline in f:
            sline = sline.strip()
            line = sline.split()
            if firstline:
                _classifierType = line[2]
                firstline = False
            Util.debug(sline,4 )
            if self._isStartReadIndicatorMatrix(line):
                self._startReadIndicatorMatrix(line)
            elif self._isIndicatorMatrix(line):
                self._readIndicatorMatrix(line)
            elif self._isEndReadIndicatorMatrix(line):
                self._endReadIndicatroMatrix(line)
            elif self._isChangeMultiClassiferId(line):
                self._changeMulticlassClassifierId(line)
            elif self._isColNameTrigger(line):
                self._triggerColName(line)
            elif self._isColNameEnd(line):
                self._endColName(line)
            elif self._isClassifierStart(line):
                self._classifierStart(line)
            elif self._isClassifierStop(line):
                self._classifierStop(line)
            elif self._isChangeClassifierId(line):
                self._changeClassifierId(line)
            elif self._isNodeLine(line):
                self._addNode(line)
            else:
                if self._colreadmode:
                    self._addColName(line)
        #set colname
        for multiclasskey,multiclassClassifiers in self._classifiers.iteritems():
            for index , classifier in multiclassClassifiers.iteritems():
                classifier._translation = self._colname
                classifier.prepare()
        f.close()
    
    def _isStartReadIndicatorMatrix(self,line):
        return line[0] == 'Indicator' and line[1]=='matrix:'
    
    def _startReadIndicatorMatrix(self,line):
        Util.debug("start reading indicator matrix",3)
        self._readIndicatorMatrixMode=True
    
    def _isIndicatorMatrix(self,line):
        return self._readIndicatorMatrixMode and Util.isInteger(line[0]) and line[1] == ':'
    
    def _readIndicatorMatrix(self,line):
        Util.debug("read indicator matrix",3)
        self._indicatorMatrix.append(map(lambda x: int(x),line[2:]))
    
    def _isEndReadIndicatorMatrix(self,line):
        return line[0] == 'Weights:'
    
    def _endReadIndicatroMatrix(self,line):
        Util.debug("end indicator matrix",3)
        self._readIndicatorMatrixMode=False
    
    def _isChangeMultiClassiferId(self,line):
        return line[0]=='Multi' and line[1] =='class' and line[2]=='learner' and line[3]=='classifier:'
    
    def _changeMulticlassClassifierId(self,line):
        Util.debug('change classifier id',3)
        self._currentMulticlassClassifierId = int(line[4])
        if not self._currentMulticlassClassifierId in self._classifiers.keys():
            self._classifiers[self._currentMulticlassClassifierId] = {}
    
    def _isColNameTrigger(self,line):
        return line[0]=='Dimensions:'
    
    def _triggerColName(self,line):
        Util.debug("start col name",3)
        self._colreadmode=True
    
    def _addColName(self,line):
        Util.debug("add column name",3)
        self._colname[int(line[0])]=line[1]
    
    def _isColNameEnd(self,line):
        return self._colreadmode and line[0]=='=================================================='
    
    def _endColName(self,line):
        Util.debug("end colname",3)
        self._colreadmode = False
    
    def _isClassifierStart(self,line):
        return line[0]=="Nodes:"
    
    def _classifierStart(self,line):
        Util.debug("start classifier",3)
        self._classifierReadMode=True
    
    def _isClassifierStop(self,line):
        return line[0]=='=================================================='
    
    def _classifierStop(self,line):
        Util.debug("stop classifier",3)
        self._classifierReadMode=False
    
    def _isChangeClassifierId(self,line):
        return line[0]=="Classifier"
    
    def _changeClassifierId(self,line):
        self._createClassifier(line)
        Util.debug("change classifier id",3)
        self._currentClassifierId=int(line[1])
    
    def _isNodeLine(self,line):
        return line[0]=='Id:'
    
    def _addNode(self,line):
        Util.debug("addNode",3)
        self._classifiers[self._currentMulticlassClassifierId][self._currentClassifierId].addNode(line)
    
    def _createClassifier(self,line):
        Util.debug("create",3)
        cid = int(line[1])
        if not cid in self._classifiers[self._currentMulticlassClassifierId].keys():
            self._classifiers[self._currentMulticlassClassifierId][cid] = DecisionTree(line,self._classifierType)
    
class TreeMaker:
    """
        finite state machine to parse a decision tree
    """
    _filename  = ''
    _classifierType = ''
    _colreadmode = False
    _classifierReadMode = False
    _currentClassifierId = -1
    _classifierCreateMode = False
    def __init__(self,filename):
        self._colname = {}
        self._classifiers = {}
        self._filename = filename
    
    #finite state machine implementation
    def parse(self):
        f = open(self._filename)
        firstline = True
        
        for sline in f:
            sline = sline.strip()
            line = sline.split()
            if firstline:
                _classifierType = line[2]
                firstline = False
            Util.debug(sline,4 )
            if self._isColNameTrigger(line):
                self._triggerColName(line)
            elif self._isColNameEnd(line):
                self._endColName(line)
            elif self._isClassifierStart(line):
                self._classifierStart(line)
            elif self._isClassifierStop(line):
                self._classifierStop(line)
            elif self._isChangeClassifierId(line):
                self._changeClassifierId(line)
            elif self._isNodeLine(line):
                self._addNode(line)
            else:
                if self._colreadmode:
                    self._addColName(line)
        #set colname
        for key,classifier in self._classifiers.iteritems():
            classifier._translation = self._colname
            classifier.prepare()
        f.close()
            
    # def getLegendWriter(self):
    #   return [ LegendWriter(self._classifiers) ]
    
    def _isColNameTrigger(self,line):
        return line[0]=='Dimensions:'
    
    def _triggerColName(self,line):
        Util.debug("start col name",3)
        self._colreadmode=True
    
    def _addColName(self,line):
        Util.debug("add column name",3)
        self._colname[int(line[0])]=line[1]
    
    def _isColNameEnd(self,line):
        return self._colreadmode and line[0]=='=================================================='
    
    def _endColName(self,line):
        Util.debug("end colname",3)
        self._colreadmode = False
    
    def _isClassifierStart(self,line):
        return line[0]=="Nodes:"
    
    def _classifierStart(self,line):
        Util.debug("start classifier",3)
        self._classifierReadMode=True
    
    def _isClassifierStop(self,line):
        return line[0]=='=================================================='
    
    def _classifierStop(self,line):
        Util.debug("stop classifier",3)
        self._classifierReadMode=False
    
    def _isChangeClassifierId(self,line):
        return line[0]=="Classifier"
    
    def _changeClassifierId(self,line):
        self._createClassifier(line)
        Util.debug("change classifier id",3)
        self._currentClassifierId=int(line[1])
    
    def _isNodeLine(self,line):
        return line[0]=='Id:'
    
    def _addNode(self,line):
        Util.debug("addNode",3)
        self._classifiers[self._currentClassifierId].addNode(line)
    
    def _createClassifier(self,line):
        Util.debug("create",3)
        cid = int(line[1])
        if not cid in self._classifiers.keys():
            self._classifiers[int(line[1])]=DecisionTree(line,self._classifierType)


class DecisionTreeNode:
    """
        Decision node parsing a line and make the node
    """
    _id=0
    _score=0
    _dim=0
    _cut=0
    _lchild=-1
    _rchild=-1
    _n0=0
    _n1=0
    _w0=0
    _w1=0
    _sizeScale = 1
    _isNewVersionNode=False
    _cutOff = 0.5
    _deltaFom = 0.0
    
    def __init__(self,line):
        #Id: 0 Score: 0 Dim: 1 Cut: 0.894776 Daughters: 1 2
        self._id = int(line[1])
        self._score = float(line[3])
        self._dim = int(line[5])
        self._cut = float(line[7])
        self._lchild = int(line[9])
        self._rchild = int(line[10])
        if self.isNewVersionNode(line):
            self._isNewVersionNode=True
            self._n0 = int(line[12])
            self._n1 = int(line[14])
            self._w0 = float(line[16])
            self._w1 = float(line[18])
            self._deltaFom = float(line[20])
    
    def isNewVersionNode(self,line):
        return len(line)>11
    
    def label(self,translation):
        if self.isTerminal():
            return 'signal' if self._score>0 else 'background'
        else:
            return translation[self._dim]
                
    def shape(self):
        return "ellipse" if self.isTerminal() else "diamond"
    
    def isTerminal(self):
        return self._dim<0
    def isSignal(self):
        return self._score>0.5
    def isNewVersion(self):
        return self._isNewVersionNode
    
    def sumN(self):
        return self._n0 + self._n1
    def sumW(self):
        return self._w0 + self._w1

class DecisionTree: 
    """
        a class for decision tree
        for most of the class's field either look up the meaning from
        SPR file format or graphviz dot reference 
    """
    _beta = 1.0 
    _type = "" #what kind of classifier is this one
    _minN = 0
    _maxN = 0
    _minW = 0
    _maxW = 0
    _minW0 = 0
    _maxW0 = 0
    _minW1 = 0
    _maxW1 = 0
    _minScore = 0
    _maxScore = 0
    _maxHeight = 1.5
    _minHeight = 0.3
    _maxWidth = 1.5
    _minWidth = 0.3
    _minPenWidth = 0.5
    _maxPenWidth = 6.0
    _minArrowSize = 1
    _maxArrowSize = 3
    # _colorscheme = 'rdylgn11'
    # _numColor = 11
    _minNodeLineWidth = 1.0
    _maxNodeLineWidth = 5.0
    # _shapeArray = ['polygon','triangle','diamond','trapezium','house',
    #   'pentagon','hexagon','septagon','octagon',
    #   'invtriangle','invtrapezium','invhouse']
    _minDimUsed = 0
    _maxDimUsed = 0
    _minColUsed = 0
    _maxColUsed = 0
    _minDeltaFom = 0
    _maxDeltaFom = 0
    #_borderColor = ["#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858"]
    #_borderColor = ["#7f3b08","#b35806","#e08214","#fdb863","#fee0b6","#f7f7f7","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b"]
    _borderColor = 'gray16'
    def __init__(self,line=None,ttype=None):
        #Classifier      0 TopdownTree Beta: 1.4094040719
        if not line is None:
            if len(line)>=4:
                self._beta = float(line[4])
        self._nodes = {}
        self._type = ttype
        self._translation = {}
        self._colUsed = {}
        self._deltaFom = None
        
    def addNode(self,line):
        """
            add a node to this tree
        """
        thisId = int(line[1])
        self._nodes[thisId] = DecisionTreeNode(line)
        
    def prepare(self):
        """
            after adding all the node compute all the constant needed
            to draw a tree
        """
        #find normalization factors
        nList = map(lambda x: self._nodes[x].sumN(),self._nodes)
        wList = map(lambda x: self._nodes[x].sumW(),self._nodes)
        deltaFomList = map(lambda x:self._nodes[x]._deltaFom,self._nodes)
        weightedDeltaFomList = map(lambda x,y: x*y, wList, deltaFomList)
        scoreList = map(lambda x: self._nodes[x]._score,self._nodes)
        self._minN = min(nList)
        self._maxN = max(nList)
        self._minW = min(wList)
        self._maxW = max(wList)
        self._minScore = min(scoreList)
        self._maxScore = max(scoreList)
        self._minDeltaFom = min(weightedDeltaFomList)
        self._maxDeltaFom = max(weightedDeltaFomList)
        
        striphead = filter(lambda x: x!=0,self._nodes)#head screwed up the normalization for edges
        w0List = map(lambda x: self._nodes[x]._w0,striphead)
        w1List = map(lambda x: self._nodes[x]._w1,striphead)
        self._minW0 = min(w0List)
        self._maxW0 = max(w0List)
        self._minW1 = min(w1List)
        self._maxW1 = max(w1List)
        
        #count column used
        dimList = map(lambda x: self._nodes[x]._dim,self._nodes)
        dimList = filter(lambda x: x>=0, dimList)
        self._maxDimUsed = max(dimList)
        self._minDimUsed = min(dimList)
        self._colUsed = map(lambda x: dimList.count(x),range(len(self._translation)))
        self._minColUsed = min(self._colUsed)
        self._maxColUsed = max(self._colUsed)
        
    def writeDotFile(self,fname):
        f = open(fname,'w')
        Util.debug('Writing: '+fname,-1)
        f.write(self.toDot())
        f.close()
    
    def toDot(self):
        toReturn = 'digraph title{\n'
        toReturn += 'graph [ rankdir="LR"'
        toReturn += ', bgcolor="%s"'%self.toBackground()
        toReturn += ', pad=%f'%Setting.margin
        toReturn += ' ]\n'
        for thisId, node in self._nodes.iteritems():
            #nodelabel
            toReturn += self.toNodeLine(node)
            
        toReturn += '}\n'
        return toReturn
    
    def toBackground(self):
        return 'transparent' if Setting.transparentBackground else 'white'
    
    def toNodeLine(self,node):
        toReturn = ''
        toReturn += '%d ['%node._id
        toReturn += 'label="%s"'%self.toLabel(node)
        toReturn += ', style="%s"'%self.toNodeStyle(node)
        toReturn += ', fillcolor="%s"'%self.toColor(node)
        toReturn += ', color="%s"'%self.toNodeColor(node)
        toReturn += ', shape=%s'%self.toShape(node)
        if node.isNewVersion():
            toReturn += ', fixedsize=true'
            toReturn += ', height=%f'%self.toHeight(node)
            toReturn += ', width=%f'%self.toWidth(node)
        toReturn += ']\n'
        
        if not node.isTerminal():#edges
            if not node.isNewVersion():
                toReturn += '%d -> %d [label="<%5.2f"]\n'%(node._id,node._lchild,node._cut)
                toReturn += '%d -> %d [label=">%5.2f"]\n'%(node._id,node._rchild,node._cut)
            else:
                toReturn += self.toEdgeLine(node)
        return toReturn
    
    def toEdgeLine(self,node):
        toReturn = '';
        toReturn += self.toEdgeLineHelper(node,True,True);
        toReturn += self.toEdgeLineHelper(node,False,True);
        toReturn += self.toEdgeLineHelper(node,False,False);
        toReturn += self.toEdgeLineHelper(node,True,False);
        return toReturn
        
    def toEdgeLineHelper(self,node,signal,left):
        toReturn = ''
        toReturn += '%d -> %d ['%(node._id,node._lchild if left else node._rchild)
        toReturn += ' color="%s"'%(Setting.signalEdgeColor if signal else Setting.backgroundEdgeColor)
        if signal and left:
            toReturn += ' label="<%5.2f"'%node._cut
        elif signal and not left:
            toReturn += ' label=">%5.2f"'%node._cut
        toReturn += ', penwidth="%f"'%self.toPenWidth(node,signal,left)
        toReturn += ', arrowhead="%s"'%('normal' if signal else 'dot')
        toReturn += ', arrowsize=%f'%self.toArrowSize(node,signal,left)
        toReturn += ' ]\n'
        return toReturn
    
    def toNodeColor(self,node):
        return '#292929'
        # if node.isTerminal(): return 'black'
        #       toReturn = int(round(Util.normalize(node._deltaFom,self._minDeltaFom,self._maxDeltaFom,0,len(self._borderColor)-1)))
        #       return self._borderColor[toReturn]
        #       
    def toShape(self,node):
        if node.isTerminal(): 
            return 'circle'
        else:
            numShape = len(Setting.shapes)
            return Setting.shapes[node._dim%numShape]
    
    def toNodeStyle(self,node):
        toReturn = 'filled'
        toReturn +=',setlinewidth(%4.2f)'%self.toNodeLineWidth(node)
        toReturn +=''
        return toReturn
    
    def toNodeLineWidth(self,node):
        if node.isTerminal(): return 1
        weightedDeltaFom = node._deltaFom * node.sumW()
        toReturn = Util.normalize(weightedDeltaFom,self._minDeltaFom,self._maxDeltaFom,self._minNodeLineWidth,self._maxNodeLineWidth)
        return toReturn
        # return 1
    def toArrowSize(self,node,signal,left):
        interestedNode = self._nodes[node._lchild] if left else self._nodes[node._rchild]
        value = interestedNode._w1 if signal else interestedNode._w0
        minw = self._minW1 if signal else self._minW0
        maxw = self._maxW1 if signal else self._maxW0
        toReturn = Util.normalize(value,minw,maxw,self._minArrowSize,self._maxArrowSize)
        return toReturn
    
    def toPenWidth(self,node,signal,left):
        interestedNode = self._nodes[node._lchild] if left else self._nodes[node._rchild]
        value = interestedNode._w1 if signal else interestedNode._w0
        minw = self._minW1 if signal else self._minW0
        maxw = self._maxW1 if signal else self._maxW0
        toReturn = Util.normalize(value,minw,maxw,self._minPenWidth,self._maxPenWidth)
        return toReturn 
    
    def toColor(self,node):
        normalizedScore = Util.normalize(node._score,self._minScore,self._maxScore,0,len(Setting.colors)-1)
        index = int(round(normalizedScore))
        return Setting.colors[index]
    
    def toHeight(self,node):
        w = node.sumW()
        toReturn = Util.normalize(w, self._minW, self._maxW, self._minHeight, self._maxHeight)
        return toReturn;
    
    def toWidth(self,node):
        w = node.sumW()
        toReturn = Util.normalize(w, self._minW, self._maxW, self._minWidth, self._maxWidth)
        return toReturn
    
    def toLabel(self,node):
        if not node.isTerminal():
            return self._translation[node._dim]
        else:
            return 's' if self.isSignal(node) else 'b'
    
    def isSignal(self,node):
        return node._score > (self._maxScore - self._minScore)/2

class LegendWriter:
    #this one doesn't look nice but I'll just leave it as is
    minHeight = 0.3
    maxHeight = 1.5
    minWidth = 0.3
    maxWidth = 1.5
    minCount = 0
    maxCount = 0
    nodePerCol = 5
    numColor = 8
    #colorScheme = 'bugn9'

    def __init__(self, trees):
        self.translation = trees[0]._translation
        self.count=[0]*len(self.translation)
        self.weightedDeltaFom=[0]*len(self.translation)
        #count number of how many time each variable was used
        
        for key,tree in trees.iteritems():
            colUsed = tree._colUsed
            self.count = map(lambda x,y: x+y, self.count,colUsed )
        self.minCount = min(self.count)
        self.maxCount = max(self.count)
        self.initWeightedDeltaFom(trees)
    
    # def writeDotFile(self,fname):
    #   f = open(fname,'w')
    #   f.write(self.toDot())
    #   f.close()
    
    def _writePieChart(self,fname,valueList,minOther = 0.05):
        class tuplehelper:
                value=0
                name =''
                def __init__(self,value,name):
                    self.value = value
                    self.name = name
        class triplehelper:
                value=0
                name =''
                color=''
                def __init__(self,value,name,color):
                    self.value = value
                    self.name = name
                    self.color = color
        pylab.clf()
        countObj = map(lambda index,value: triplehelper(value,self.translation[index],Setting.piecolor[index%len(Setting.piecolor)]), \
            xrange(len(valueList)), valueList)
        #now try to collect all the small one into "others"
        vList = map(lambda x:x.value,countObj)
        sumValue = sum(vList)
        cutoff = minOther*sumValue
        failCutoff = filter(lambda x: x.value < cutoff, countObj) #this go into others
        passCutoff = filter(lambda x: x.value >= cutoff, countObj) #this will be shown individually
        passCutoff.sort(lambda a,b: 1 if b.value > a.value else -1) #now sort
        if len(failCutoff) != 0: #there are some small stuff to collect
            sumFail = reduce(lambda x,y:x+y.value,failCutoff,0.0)
            other = triplehelper(sumFail,'Other',Setting.pieOtherColor)
            passCutoff.append(other)
        
        myLabels = map(lambda x: x.name, passCutoff)
        fracs = map(lambda x: float(x.value), passCutoff)
        colors = map(lambda x: x.color, passCutoff)
        #print labels
        #print fracs
        pylab.pie(fracs,labels=myLabels, colors=colors, autopct='%1.1f%%', shadow=True)
        pylab.savefig(fname)
    
    def initWeightedDeltaFom(self,trees):
        self.weightedDeltaFom
        for index , tree in trees.iteritems():
            for nodeIndex, node in tree._nodes.iteritems():
                if node.isTerminal():
                    continue
                self.weightedDeltaFom[node._dim] += node.sumW()*node._deltaFom
    
    def writeWeightedDeltaFomPieChart(self,fname,minOther=0.05):
        self._writePieChart(fname,self.weightedDeltaFom,minOther)
    
    def writeLegendFile(self,fname):
        f = open(fname,'w')
        f.write(self.toLegendFile())
        f.close()
    
    def writeSplitCountPieChart(self,fname,minOther=0.05):
        self._writePieChart(fname,self.count,minOther)
    
    # def toDot(self):
    #   numCol = len(self.count)
    #   opened = False
    #   first = False
    #   toReturn = "graph legend{\n"
    #   toReturn += 'graph [rankdir="LR"]'
    #   nodePerCol = self.nodePerCol
    #   for index in xrange(numCol):
    #       # if index%nodePerCol == 0:
    #       #   opened = True
    #       #   first = True
    #       #   toReturn += "subgraph cluster%d{"%index
    #       toReturn += self.toLine(index)
    #       # if index%nodePerCol == nodePerCol-1:
    #       #   opened = False
    #       #   toReturn += "}"
    #   # if opened:
    #   #   opened = False
    #   #   toReturn += "}"         
    #   toReturn +="}\n"
    #   return toReturn
    
    # def toLine(self,index):
    #   toReturn = '"%s" ['%(self.translation[index])
    #   toReturn += "fixedsize=true"
    #   toReturn += ", height=%f"%self.toHeight(index)
    #   toReturn += ", width=%f"%self.toWidth(index)
    #   toReturn += ', shape=%s'%self.toShape(index)
    #   toReturn += ', style="filled"'
    #   toReturn += ', fillcolor=%d'%self.toColor(index)
    # 
    #   toReturn += "]\n"
    #   return toReturn
    
    # def toColor(self,index):
    #   index = round(Util.normalize(self.count[index],self.minCount,self.maxCount,1,self.numColor))
    #   return toReturn
    
    def toShape(self,index):
        return Setting.shapes[index%len(Setting.shapes)]
    
    def toWidth(self,index):
        return Util.normalize(self.count[index],self.minCount,self.maxCount,self.minWidth,self.maxWidth)
    
    def toHeight(self,index):
        return Util.normalize(self.count[index],self.minCount,self.maxCount,self.minHeight,self.maxHeight)
    
    def toLegendFile(self): #primitive legend in postscript format
        s = '%!PS-Adobe-2.0 EPSF-2.1\n\
%%BoundingBox: 0 0 550 200\n\
%%EndComments\n'
        s += '/colors ' + Util.colorPSArray() + ' def\n'
        s += '/strokes ' + '[ 1.0 2.0 4.0 6.0 ]' + ' def\n'
        s += '/linelength 70 def\n\
/cwidth 550 def\n\
/size 20 def\n\
/fontsize 15 def\n\
/offset 20 def\n\
/debugint{\n\
    50 string cvs 20 20 moveto show\n\
} def\n\
/box {% [r g b]\n\
    gsave\n\
    aload\n\
    pop\n\
    setrgbcolor\n\
    currentpoint\n\
    /orgy exch def\n\
    /orgx exch def\n\
    newpath\n\
    orgx orgy moveto\n\
    orgx size add orgy lineto\n\
    orgx size add orgy size add lineto\n\
    orgx orgy size add lineto\n\
    orgx orgy lineto\n\
    fill\n\
    closepath\n\
    grestore\n\
    orgx size add orgy moveto\n\
} def\n\
/movetocolorcenter{\n\
    cwidth\n\
    size\n\
    colors length\n\
    mul\n\
    sub\n\
    2\n\
    div\n\
    50\n\
    moveto\n\
} def\n\
/centershow { %show center text\n\
    dup stringwidth pop 2 div neg 0 rmoveto show\n\
} def\n\
/centershowreturn {\n\
    currentpoint\n\
    /orgy exch def\n\
    /orgx exch def\n\
    centershow\n\
    orgx orgy moveto\n\
} def\n\
/showoffsetlabeldown {\n\
    0 offset neg rmoveto\n\
    centershowreturn\n\
    0 offset rmoveto\n\
} def\n\
/showoffsetlabelup {\n\
    dup \n\
    stringwidth /offsetcorrection exch def % now offsetcorrection\n\
    pop %cleanup\n\
    0 offset offsetcorrection sub rmoveto\n\
    centershowreturn\n\
    0 offset offsetcorrection sub neg rmoveto\n\
} def\n\
/drawflows{ %[r g b]\n\
    aload\n\
    pop\n\
    gsave\n\
    setrgbcolor\n\
    strokes\n\
    {\n\
        setlinewidth\n\
        currentpoint\n\
        /orgy exch def\n\
        /orgx exch def\n\
        newpath\n\
            orgx orgy moveto\n\
            linelength 0 rlineto\n\
            stroke\n\
        closepath\n\
        linelength orgx add linelength 3 div add orgy moveto\n\
    }\n\
    forall\n\
    currentpoint\n\
    /tempy exch def\n\
    /tempx exch def\n\
    grestore\n\
    tempx tempy moveto\n\
} def\n\
/drawsignalflow {\n\
    (less signal flow)showoffsetlabeldown\n\
    colors colors length 1 sub get drawflows\n\
    (more signal flow)showoffsetlabeldown\n\
} def\n\
/drawbackgroundflow {\n\
    (less background flow)showoffsetlabelup\n\
    colors 0 get drawflows\n\
    (more background flow)showoffsetlabelup \n\
} def\n\
/Times-Roman findfont fontsize scalefont setfont\n\
movetocolorcenter\n\
(background) showoffsetlabeldown\n\
colors {box} forall\n\
(signal) showoffsetlabeldown\n\
90 150 moveto\n\
drawbackgroundflow\n\
90 140 moveto\n\
drawsignalflow\n\
showpage\n'
        return s

main()
# testMultiClass()