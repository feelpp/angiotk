#!/usr/bin/env python                                                                       

import argparse
import ConfigParser
import datetime
import errno
import glob
import os
import shutil
import stat
import subprocess
import sys
import threading
import time
# to generate results database (json file)                                                  
import json
from collections import OrderedDict
import string

import testAngioTKPipeline

def main():
    
    testAngioTK = os.path.expandvars('$INSTALL/bin/testAngioTKPipeline.py')
    asciidocpy = os.path.expandvars('$INSTALL/bin/createAsciidocDB.py')

    initialTime = time.time()
    # arguments management
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputpath', default=os.path.expandvars("/data/vivabrain"), help='path to all input files in .mha format')
    parser.add_argument('--outputpath', default=os.path.expandvars("$RESULTS"), help='path to write results to')
    parser.add_argument('-n', default=1, type=int, help='limit actual processing to n mha files')
    args = parser.parse_args()
    
    inputFile = args.inputpath
    resultsPath = args.outputpath
    #fileID = string.lstrip(args.inputfile, '/data/vivabrain/IRM/')
    mhaFiles = []
    osWalkInitTime = time.time()
    for path, dirs, files in os.walk(inputFile):
        for f in files:
            if f.endswith("art.mha"):
                mhaFiles.append(os.path.join(path,f))
    osWalkEndTime = time.time()
    print "The following .mha files were found"
    for f in mhaFiles:
        print f
        
    print "os.walk total time: " + str(round(osWalkEndTime-osWalkInitTime,1))
    
    i = 0
    limit = min(args.n, len(mhaFiles))
    for mhaFile in mhaFiles: # ex: /data/vivabrain/IRM/Base2/TOF/natives/MHA/GEO26/TOF_art.mha
        # pathNoExt, ext = os.path.splitext(mhaFile) # ex: separate /data/vivabrain/IRM/Base2/TOF/natives/MHA/GEO26/TOF_art and .mha 
        pathToMha = string.rstrip(string.lstrip(mhaFile, '/data/vivabrain/IRM/')) # only keep Base2/TOF/natives/MHA/GEO26/TOF_art.mha
        fileIDNoSlash = pathToMha.replace('/','_').replace('.','_') # turn it into Base2_TOF_natives_MHA_GEO26_TOF_art_mha
        inputPath = os.path.join(os.path.expandvars('$INSTALL/share/AngioTK/Examples/MeshFromMRI'), fileIDNoSlash) #
        testAngioTKPipeline.makedir(inputPath)
        cmd = "python " + testAngioTK + ' --inputfile=' + mhaFile + ' --inputpath='+ inputPath  + ' --outputpath=' + os.path.join(resultsPath, fileIDNoSlash) 
        print 'Processing file ' + str(i+1) + '/' + str(limit) + ': ' + mhaFile + '...'
        print 'cmd: ' + cmd
        print '----------------------------------------' + testAngioTK + '--------------------------------------------'
        testAngioTKPipeline.executeCommand("", cmd, "")
        i=i+1
        if(i==args.n):
            break

    resultsDBPath = os.path.join(resultsPath, 'resultsDataBase')
    cmd = 'python ' + asciidocpy + ' -i ' + resultsDBPath + '/resultsDB.json'
    print '----------------------------------------' + asciidocpy + '--------------------------------------------'
    testAngioTKPipeline.executeCommand("", cmd, "")
    cmd = 'asciidoctor ' + resultsDBPath + '/resultsDB.adoc'
    print '----------------------------------------' + 'asciidoctor' + '--------------------------------------------'
    testAngioTKPipeline.executeCommand("", cmd, "")
        
# Only do this, if we do not import as a module
if __name__ == "__main__":
    main()
