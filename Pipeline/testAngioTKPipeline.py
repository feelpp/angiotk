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

debug = 1



# execute a command line
def executeCommand(text, cmd, logPrefix=""):
    out = ""
    returnCode = 0

    if(debug > 0):
        print cmd

    # if we specify an absolute command path
    # we check that the file exists
    if(os.path.isabs(cmd[0])):
        if(not os.path.exists(cmd[0])):
            return 1, ""

    # Ensure that we are not overwriting previous files by appending an index
    logIndex = 0
    while(os.path.exists(logPrefix + "." + str(logIndex) + ".log")):
        logIndex = logIndex + 1

    # Log execution time
    tstart = time.time()
    try:
        # Launch the current command and also redirect its output to log files
        out = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        tee = subprocess.Popen(['tee', logPrefix + "." + str(logIndex) + ".log"], stdin=out.stdout)
        out.stdout.close()
        tee.communicate()
    except subprocess.CalledProcessError, e:
        out = e.output
        ret = e.returncode
    # Log execution time
    tend = time.time()
    tFinal = (tend - tstart)

    # Log command parameters if needed
    if(logPrefix != 0):
        execlog = open(logPrefix + "." + str(logIndex) + ".exec.log", "w")
        execlog.write("cmd=\"" + cmd + "\"\n")
        execlog.write("elapsedTime=" + str(tFinal) + "\n")
        execlog.close()

    return returnCode, tFinal

# make a directory (don't fail if it exists)
def makedir(path):
    print path
    print " "
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def sanityCheckFile(filepath, msg=""):
    if( not os.path.exists(filepath) ):
        if(msg != ""):
            print msg + "\nAborting ..."
        else:
            print "The file " + filepath + " was not found\nAborting ..."
        exit(1)
# add a new entry to the results database (an ordered dictionary)
def addStepEntry(resDB, currentStepName, success, time, screenshotfile, comments):
    manualCommentHint = "Users can edit and add comments here."
    if(not resDB.has_key(currentStepName)):
        resDB[currentStepName] = OrderedDict()
    resDB[currentStepName]['Success'] =  success
    resDB[currentStepName]['Time'] = round(time,3)
    resDB[currentStepName]['Screenshot'] = screenshotfile
    resDB[currentStepName]['Comments']= comments
    if(not resDB[currentStepName].has_key('User comments')):
        resDB[currentStepName]['Users comments'] = manualCommentHint

# make a screenshot of model and save it to screenshotPath
def saveScreenshot(modelPath, screenshotPath):
    # path to screenshot script:
    screenshotScriptPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "screenshot.py")
    screenshotCmd = "pvbatch --use-offscreen-rendering " + screenshotScriptPath + " " + modelPath + " " + screenshotPath
    print '--------------------------------------------------------------------------------'
    executeCommand("", screenshotCmd, "")
    
    
def main():
    # Check that the meshing and filtering modules are enabled
    hasMeshingModule = "ON"
    hasFilteringModule = "ON"

    if((hasMeshingModule == "ON" or hasMeshingModule == "TRUE") \
        and hasFilteringModule == "ON" or hasFilteringModule == "TRUE"):
        print "Modules OK"
    else:
        print "The Filtering and/or Meshing modules were not enabled in AngioTK (-DBUILD_MODULE_*)"
        print "Aborting ..."
        exit(1)

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputfile', required=True, help='Path to the original input file in .mha format')
    parser.add_argument('--inputpath', required=True, help='Path where the configuration files for the current pipeline are stored')
    parser.add_argument('--outputpath', default=os.path.expandvars("/AngioTK/ProcessedData"), help='Path where intermediary data is written')
    args = parser.parse_args()

    #if len(sys.argv) != 3:
        #print "usage: " + os.path.basename(sys.argv[0]) + " <path_to_config_files> <output_path;default=" + outpath + ">"
        #exit(1)

    # Stores useful paths
    # cfgpath = "/ssd/derhovsepian/angiotk/install/share/AngioTK/"
    # outpath = os.path.expandvars("/AngioTK/")

    #cfgpath = os.path.abspath(sys.argv[1])
    cfgpath = args.inputpath
    sanityCheckFile(cfgpath)

    #outpath = os.path.abspath(sys.argv[2])
    baseOutputPath = args.outputpath
    makedir(baseOutputPath)

    script = open(baseOutputPath + "/relaunch-script.sh", "w")
    script.write("#!/bin/bash\n\n")
    script.write(' '.join(sys.argv) + "\n\n")
    script.close()

    # make script executable
    st = os.stat(baseOutputPath + "/relaunch-script.sh")
    os.chmod(baseOutputPath + "/relaunch-script.sh", st.st_mode | stat.S_IEXEC)

    AngiotTKOutputPath = baseOutputPath + "/angiotk"
    makedir(AngiotTKOutputPath + "/data")

    print "Using configuration:"
    print "- Location of configuration files: " + cfgpath
    print "- Output path: " + AngiotTKOutputPath

    # Setup environment to find the correct paths and set output paths
    my_env = os.environ
    my_env["PATH"] = "/ssd/derhovsepian/angiotk/install/Modules/Meshing/bin:" + my_env["PATH"]
    my_env["PATH"] = "/ssd/derhovsepian/angiotk/install/Modules/Filtering/RORPO/bin:" + my_env["PATH"]
    my_env["PATH"] = "/ssd/derhovsepian/angiotk/install/Modules/Filtering/ITK/bin:" + my_env["PATH"]
    my_env["PATH"] = "/ssd/derhovsepian/angiotk/install/bin/AngioTK/Data:" + my_env["PATH"]
    my_env["LD_LIBRARY_PATH"] = "/ssd/derhovsepian/feelpp/install/lib:" + my_env["LD_LIBRARY_PATH"]
    my_env["LD_LIBRARY_PATH"] = "/ssd/derhovsepian/angiotk/install/Modules/Filtering/RORPO/lib:" + my_env["LD_LIBRARY_PATH"]
    my_env["LD_LIBRARY_PATH"] = "/ssd/derhovsepian/angiotk/install/Modules/Meshing/lib:" + my_env["LD_LIBRARY_PATH"]
    my_env["OUTPUT_PATH"] = baseOutputPath
    my_env["FEELPP_WORKDIR"] = my_env["OUTPUT_PATH"] + "/feel"

    print "- PATH=" + my_env["PATH"]
    print ""

    # Create a log directory for the current execution
    dt = datetime.datetime.now()
    logPath = baseOutputPath + "/log/" + dt.strftime("%Y%m%d-%H%M%S")
    makedir(logPath)

    # Start timing
    tinit = time.time()

    inputFileNoExt, inputFileExt = os.path.splitext(os.path.basename(args.inputfile))

    # Check that input file exists
    sanityCheckFile(args.inputfile)
    
    # ------------ Results Database initialization -------------
    # create results paths
    resultsRoot, datasetDir = os.path.split(baseOutputPath) # ex: $RESULTS/Phantom -> $RESULTS and Phantom
    if baseOutputPath.endswith('/'): # in case baseOutpuPath ends with a slash '/', split again
        resultsRoot, datasetDir = os.path.split(resultsRoot)
    print 'resultsRoot: ' + resultsRoot # ex: $RESULTS
    print 'datasetDir: ' + datasetDir # ex: Phantom
    resultsPath = os.path.join(resultsRoot, "resultsDataBase") # ex: $RESULTS/resultsDataBase
    resultsDBName = "resultsDB.json"
    resultsDBPath = os.path.join(resultsPath, resultsDBName) # ex: $RESULTS/resultsDB.json
    # initialize database by reading existing file or creating one if needed.
    resultsDB = OrderedDict() # Data structure: ordered dictionary
    try: # if .json already exists, read it 
        dbFile = open(resultsDBPath,'r')
        resultsDB = json.load(dbFile, object_pairs_hook=OrderedDict)
        print json.dump(resultsDB, indent=4)
        dbFile.close()
    except:
        print "****RESULTS DATABASE:"
        print '\t > The file ' + resultsDBPath + ' was not found, it will be created'
    # adds an entry for the inputfile in the database if needed.
    fileID = string.lstrip(args.inputfile, '/data/vivabrain/IRM/')
    if(not resultsDB.has_key(fileID)):
        print "\t > Adding " + args.inputfile + " entry..."
        resultsDB[fileID] = OrderedDict()
    resDBfileEntry = resultsDB[fileID]
    
    stepsNames = ['1 - RORPO processing', '2 - Mesh extraction from image', '3 - Centerline computing', '4 - Image generation from centerlines', '5 - Second mesh extraction from image', '6 - Volume mesh processing']
    for stepName in stepsNames:
        if(not resDBfileEntry.has_key(stepName)):
           addStepEntry(resDBfileEntry, stepName, False, 0.0, 'none yet', 'Step not reached.')
    # write changes to the database file and closes it.                                                                                                                                                                                      
    dbFile = open(resultsDBPath, 'w')
    json.dump(resultsDB, dbFile, indent=4)
    dbFile.close()

    # ------------ Screenshots settings -----------
    # file extension to use
    screenshotExtension = ".jpg"
    # path to store important step images
    fileIDNoSlash = fileID.replace('/','_') # ex: CEMRACS_2015/Phantom -> CEMRACS_2015_Phantom to keep dataset id and avoid deep directory trees
    screenshotPath = os.path.join("screenshots/", fileIDNoSlash) # relative link, for json/asciidoc/html generation
    resultsScreenshotPath = os.path.join(resultsPath, screenshotPath) # relative link, to write screenshot files 
    makedir(resultsScreenshotPath) # ex: create $RESULTS/resultsDataBase/screenshots if needed
    print 'screenshotPath: ' + screenshotPath
    print 'resultsScreenshotPath: ' + resultsScreenshotPath
    
    # Ensure that the input data file has the correct format and copy it to the data directory
    # otherwise convert it
    print "Checking input file " + args.inputfile + " ... "
    if(inputFileExt != ".mha"):
        cmd = "vmtkimagereader -ifile " + args.inputfile + " " \
            "--pipe vmtkimagewriter -ofile " + AngiotTKOutputPath + "/data/" + inputFileNoExt + ".mha"
        retcode, _ = executeCommand("", cmd, logPath + "/Convert_format_mha")
    else:
        shutil.copy(args.inputfile, AngiotTKOutputPath + "/data/")
    sanityCheckFile(AngiotTKOutputPath + "/data/" + inputFileNoExt + ".mha")
    print ""

    # Get element spacing from input dataset
    # We need isotropic data for RORPO
    print "Checking if input data is isotropic ... "
    ps = subprocess.Popen(('head', '-n', '12', AngiotTKOutputPath + "/data/" + inputFileNoExt + ".mha"), stdout=subprocess.PIPE)
    spacing = subprocess.check_output(('grep', 'ElementSpacing'), stdin=ps.stdout)
    ps.wait()
    
    # clean up spacing line
    spacing = spacing.translate(None, "\n")
    lst = spacing.split()
    # Check that spacing has been correctly read
    if(len(lst) == 5):
        # Ensure spacing for isotropic data
        if(float(lst[2]) == float(lst[3]) and float(lst[3]) == float(lst[4])):
            print "Current spacing OK (" + spacing + ")"
            if(not os.path.exists(AngiotTKOutputPath + "/data/" + inputFileNoExt + "_iso.mha")):
                os.symlink(AngiotTKOutputPath + "/data/" + inputFileNoExt + ".mha", AngiotTKOutputPath + "/data/" + inputFileNoExt + "_iso.mha")
        else:
            print "Incorrect spacing (" + spacing + ")"
            print "Attempting to correct spacing"
            cmd = "ResampleVolumesToBeIsotropic " + AngiotTKOutputPath + "/data/" + inputFileNoExt + ".mha " + \
                  AngiotTKOutputPath + "/data/" + inputFileNoExt + "_iso.mha 0 255"
            if(not os.path.exists(AngiotTKOutputPath + "/data/" + inputFileNoExt + "_iso.mha")):
                retcode, _ = executeCommand("", cmd, logPath + "/Convert_iso")
            sanityCheckFile(AngiotTKOutputPath + "/data/" + inputFileNoExt + "_iso.mha")
            #print "Aborting ..."
            #exit(1)
    else:
        print "Unable to get spacing"
        exit(1)
    print ""
    
    # Convert mha data to nii for RORPO
    # Print "Converting input file to .nii ... "
    cmd = "vmtkimagereader -ifile " + AngiotTKOutputPath + "/data/" + inputFileNoExt + "_iso.mha" + " " \
        "--pipe vmtkimagewriter -ofile " + AngiotTKOutputPath + "/data/" + inputFileNoExt + ".nii"
    if(not os.path.exists(AngiotTKOutputPath + "/data/" + inputFileNoExt + ".nii")):
        retcode, _ = executeCommand("", cmd, logPath + "/Convert_mha_nii")
    sanityCheckFile(AngiotTKOutputPath + "/data/" + inputFileNoExt + ".nii")
    print ""

    # Filter data with RORPO
    scaleMin = 25
    factor = 1.34
    nbScales = 7
    nbCores = 7

    # Parses the rorpo.cfg file if it exists
    tmpTime = time.time()
    print "Process data with RORPO ... "
    if(os.path.exists(cfgpath + "/rorpo.cfg")):
        Config = ConfigParser.ConfigParser()
        Config.read(cfgpath + "/rorpo.cfg")
        try: scaleMin = int(Config.get("rorpo", "scaleMin"))
        except: print("Option scaleMin not found in " + cfgpath + "/rorpo.cfg")
        try: factor = float(Config.get("rorpo", "factor"))
        except: print("Option factor not found in " + cfgpath + "/rorpo.cfg")
        try: nbScales = int(Config.get("rorpo", "nbScales"))
        except: print("Option nbScales not found in " + cfgpath + "/rorpo.cfg")

    makedir(AngiotTKOutputPath + "/RORPO/")
    RORPOOutputPrefix = inputFileNoExt + "_RORPO_" + str(scaleMin) + "_" + str(factor) + "_" + str(nbScales)
    cmd = "RORPO_multiscale_usage " + AngiotTKOutputPath + "/data/" + inputFileNoExt + ".nii " + \
        AngiotTKOutputPath + "/RORPO/" + RORPOOutputPrefix + ".nii " + \
        str(scaleMin) + " " + str(factor) + " " + str(nbScales) + " --core " + str(nbCores) + " --verbose"
    if(not os.path.exists(AngiotTKOutputPath + "/RORPO/" + RORPOOutputPrefix + ".nii")):
        retcode, _ = executeCommand("", cmd, logPath + "/RORPO")
        tmpComment = 'RORPO launched.'
    else:
        tmpComment = "RORPO not launched, target file already exists."
    sanityCheckFile(AngiotTKOutputPath + "/RORPO/" + RORPOOutputPrefix + ".nii")
    print ""
    
    # Convert filtered data back to mha
    cmd = "vmtkimagereader -ifile " + AngiotTKOutputPath + "/RORPO/" + RORPOOutputPrefix + ".nii " + \
        "--pipe vmtkimagewriter -ofile " + AngiotTKOutputPath + "/RORPO/" + RORPOOutputPrefix + ".mha"
    if(not os.path.exists(AngiotTKOutputPath + "/RORPO/" + RORPOOutputPrefix + ".mha")):
        retcode, _ = executeCommand("", cmd, logPath + "/Convert_nii_mha")
    sanityCheckFile(AngiotTKOutputPath + "/RORPO/" + RORPOOutputPrefix + ".mha")

    # ++++++ Save a screenshot of the converted  mha
    model = AngiotTKOutputPath + "/RORPO/" + RORPOOutputPrefix + ".mha"
    screenshotFile = "/1_RORPO" + screenshotExtension
    saveScreenshot(model, resultsScreenshotPath + screenshotFile)

    # ------ Add Step Entry
    addStepEntry(resDBfileEntry, '1 - RORPO processing', True, time.time() - tmpTime, screenshotPath + screenshotFile, tmpComment)
    # write changes to the database file and closes it.                                                                                                                                                                                      
    dbFile = open(resultsDBPath, 'w')
    json.dump(resultsDB, dbFile, indent=4)
    dbFile.close()

    # Extract mesh from images
    tmpTime = time.time()
    print "Extracting surface from image ... "
    sanityCheckFile(cfgpath + "/surfacefromimage.cfg")
    cmd = "meshing_surfacefromimage --config-file " + cfgpath + "/surfacefromimage.cfg " + \
        "--input.image.filename " + AngiotTKOutputPath + "/RORPO/" + RORPOOutputPrefix + ".mha " + \
        "--pre-process.resize-from-reference-image.path " + AngiotTKOutputPath + "/data/" + inputFileNoExt + "_iso.mha" + " " \
        "--output.path " + AngiotTKOutputPath + "/surfacefromimage/model.stl"
    #cmd = "meshing_surfacefromimage --config-file " + cfgpath + "/surfacefromimage.cfg " + \
    #"--input.image.filename=" + baseOutputPath + "/data/TOF_art_FOA42_RORPO_30_1.34_7.mha " + \
    #"--pre-process.resize-from-reference-image.path " + baseOutputPath + "/data/TOF_art.mha"
    retcode, _ = executeCommand("", cmd, logPath + "/meshing_surfacefromimage")
    sanityCheckFile(AngiotTKOutputPath + "/surfacefromimage/model.stl")

    # ++++++ Save a screenshot of the extracted  mesh
    model = AngiotTKOutputPath + "/surfacefromimage/model.stl"
    screenshotFile = "/2_SurfaceFromImage" + screenshotExtension
    saveScreenshot(model, resultsScreenshotPath + screenshotFile)

    # ------ Add Step Entry
    addStepEntry(resDBfileEntry, '2 - Mesh extraction from image', True, time.time() - tmpTime, screenshotPath + screenshotFile, 'Mesh extracted.')
    # write changes to the database file and closes it.                                                                                                                                                                                      
    dbFile = open(resultsDBPath, 'w')
    json.dump(resultsDB, dbFile, indent=4)
    dbFile.close()

    # Compute centerlines
    # Iterate over pointpair data
    tmpTime = time.time()
    print "Computing centerlines ... "

    sanityCheckFile(cfgpath + "/centerlines.cfg")
    clid = 0
    cltpool = []

    #cmd = "meshing_centerlines --config-file " + cfgpath + "/centerlines1.cfg " \
        #"--input.surface.filename=$repository/angiotk/surfacefromimage/model.stl " \
        #"--output.directory=angiotk/centerlines/part" + str(clid)
    #t = threading.Thread(target=executeCommand, args=("", cmd, logPath + "/meshing_centerlines9999"))
    #clid = clid + 1
    #t.start()
    #cltpool.append(t)

    cldata = glob.glob(cfgpath + "/data/pointset*.data")
    for cl in cldata:
        cmd = "meshing_centerlines --config-file " + cfgpath + "/centerlines.cfg " \
            "--input.surface.filename " + AngiotTKOutputPath + "/surfacefromimage/model.stl " \
            "--input.pointset.filename " + cl + " " \
            "--output.directory " + AngiotTKOutputPath + "/centerlines/part" + str(clid) + " " \
            "--delaunay-tessellation.output.directory " + AngiotTKOutputPath + "/centerlines/part" + str(clid) + " " \
            "--source-ids 0 --target-ids 1"

        # Create a new thread
        t = threading.Thread(target=executeCommand, args=("", cmd, logPath + "/meshing_centerlines" + str(clid)))
        t.start()
        cltpool.append(t)

        clid = clid + 1


    cldata = glob.glob(cfgpath + "/data/pointpair*.data")
    for cl in cldata:
        cmd = "meshing_centerlines --config-file " + cfgpath + "/centerlines.cfg " \
            "--input.surface.filename " + AngiotTKOutputPath + "/surfacefromimage/model.stl " \
            "--input.pointpair.filename " + cl + " " \
            "--output.directory " + AngiotTKOutputPath + "/centerlines/part" + str(clid) + " " \
            "--delaunay-tessellation.force-rebuild 0 " \
            "--delaunay-tessellation.output.directory " + AngiotTKOutputPath + "/centerlines/part" + str(clid) + "/delaunay-tessellation" 
        # Create a new thread
        t = threading.Thread(target=executeCommand, args=("", cmd, logPath + "/meshing_centerlines" + str(clid)))
        t.start()
        cltpool.append(t)

        clid = clid + 1

    # Iterate over geocenterline data
    cldata = glob.glob(cfgpath + "/data/geocenterlines*.data")
    for cl in cldata:
        cmd = "meshing_centerlines --config-file " + cfgpath + "/centerlines.cfg " \
            "--input.surface.filename " + AngiotTKOutputPath + "/surfacefromimage/model.stl " \
            "--input.geo-centerlines.filename " + cl + " " \
            "--output.directory " + AngiotTKOutputPath + "/centerlines/part" + str(clid) + " " \
            "--delaunay-tessellation.output.directory " + AngiotTKOutputPath + "/centerlines/part" + str(clid)
        # Create a new thread
        t = threading.Thread(target=executeCommand, args=("", cmd, logPath + "/meshing_centerlines" + str(clid)))
        t.start()
        cltpool.append(t)

        clid = clid + 1

    # wait for threads to finish
    for t in cltpool:
        t.join()    

    # Merge centerlines
    print "Merging centerlines ... "
    sanityCheckFile(cfgpath + "/centerlinesmanager.cfg")
    cmd = "meshing_centerlinesmanager --config-file " + cfgpath + "/centerlinesmanager.cfg "
    cmd = cmd + "--input.surface.filename " + AngiotTKOutputPath + "/surfacefromimage/model.stl "
    cmd = cmd + "--output.directory " + AngiotTKOutputPath + "/centerlinesmanager "
    for i in range(len(cltpool)):
        cmd = cmd + "--input.centerlines.filename " + AngiotTKOutputPath + "/centerlines/part" + str(i) + "/model_centerlines.vtk "
    executeCommand("", cmd, logPath + "/meshing_centerlinesmanager")

    # ++++++ Save a screenshot of the centerlines
    model = AngiotTKOutputPath + "/centerlinesmanager/model_centerlines_up.vtk"
    screenshotFile = "/3_Centerlines" + screenshotExtension
    saveScreenshot(model, resultsScreenshotPath + screenshotFile)


    # ------ Add Step Entry
    addStepEntry(resDBfileEntry, '3 - Centerline computing', True, time.time() - tmpTime, screenshotPath + screenshotFile, 'Centerlines computed.')
    # write changes to the database file and closes it.
    dbFile = open(resultsDBPath, 'w')
    json.dump(resultsDB, dbFile, indent=4)
    dbFile.close()
    
    # Regenerate image from centerlines
    tmpTime = time.time()
    print "Generating image from centerlines ... "
    sanityCheckFile(cfgpath + "/imagefromcenterlines.cfg")
    cmd = "meshing_imagefromcenterlines --config-file " + cfgpath + "/imagefromcenterlines.cfg " \
        "--input.centerlines.filename " + AngiotTKOutputPath + "/centerlinesmanager/model_centerlines_up.vtk " \
        "--output.directory " +  AngiotTKOutputPath + "/imagefromcenterlines"
    executeCommand("", cmd, logPath + "/meshing_imagefromcenterlines")
    
    # Get the output of the imagefromcenterlines algorithm
    imgGlob = glob.glob(AngiotTKOutputPath + "/imagefromcenterlines/*.mha")
    if(len(imgGlob) == 0):
        print "Cannot find output of imagefromcenterlines. Aborting."
        exit(1)
    newImage = max(imgGlob, key=os.path.getctime)

    # ++++++ Save a screenshot of the output
    model = newImage
    screenshotFile = "/4_ImageFromCenterlines" + screenshotExtension
    saveScreenshot(model, resultsScreenshotPath + screenshotFile)

    # ------ Add Step Entry
    addStepEntry(resDBfileEntry, '4 - Image generation from centerlines', True, time.time() - tmpTime, screenshotPath + screenshotFile, 'Image generated.')
    # write changes to the database file and closes it.
    dbFile = open(resultsDBPath, 'w')
    json.dump(resultsDB, dbFile, indent=4)
    dbFile.close()

    # Extract the mesh
    tmpTime = time.time()
    print "Extracting surface from image ... "
    sanityCheckFile(cfgpath + "/surfacefromimage2.cfg")
    cmd = "meshing_surfacefromimage --config-file " + cfgpath + "/surfacefromimage2.cfg " \
        "--input.image.filename " + newImage + " " + \
        "--output.directory " +  AngiotTKOutputPath + "/surfacefromimage2"
    executeCommand("", cmd, logPath + "/meshing_surfacefromimage2")

    newSTL = max(glob.iglob(AngiotTKOutputPath + "/surfacefromimage2/*.stl"), key=os.path.getctime)

    # ++++++ Save a screenshot of the extracted surface
    model = newSTL
    screenshotFile = "/5_SurfaceFromImage2" + screenshotExtension
    saveScreenshot(model, resultsScreenshotPath + screenshotFile)

    # ------ Add Step Entry
    addStepEntry(resDBfileEntry, '5 - Second mesh extraction from image', True, time.time() - tmpTime, screenshotPath + screenshotFile, 'Mesh extracted.')
    # write changes to the database file and closes it.
    dbFile = open(resultsDBPath, 'w')
    json.dump(resultsDB, dbFile, indent=4)
    dbFile.close()

    # Process surface
    tmpTime = time.time()
    print "Processing mesh ... "
    sanityCheckFile(cfgpath + "/remeshstlgmsh.cfg")
    cmd = "meshing_remeshstl --config-file " + cfgpath + "/remeshstlgmsh.cfg " \
        "--input.surface.filename " + newSTL + " "  \
        "--gmsh.centerlines.filename " + AngiotTKOutputPath + "/centerlinesmanager/model_centerlines_up.vtk " \
        "--output.directory " +  AngiotTKOutputPath + "/remeshgmsh"
    executeCommand("", cmd, logPath + "/meshing_remeshstlgmsh")

    newSTL = max(glob.iglob(AngiotTKOutputPath + "/remeshgmsh/*"), key=os.path.getctime)

    # Prepare volume data for simulation
    print "Generating volume mesh ... "
    sanityCheckFile(cfgpath + "/volumefromstlandcenterlines.cfg")
    cmd = "meshing_volumefromstlandcenterlines --config-file " + cfgpath + "/volumefromstlandcenterlines.cfg " \
        "--input.surface.filename " + newSTL + " " \
        "--input.centerlines.filename " + AngiotTKOutputPath + "/centerlinesmanager/model_centerlines_up.vtk " \
        "--output.directory " + AngiotTKOutputPath + "/volumefromstlandcenterlines"
    executeCommand("", cmd, logPath + "/meshing_volumefromstlandcenterlines")

    # ++++++ Save a screenshot of the stl model.
    model = max(glob.iglob(AngiotTKOutputPath + "/remeshgmsh/*.stl"), key=os.path.getctime)
    screenshotFile = "/6_VolumeMesh" + screenshotExtension
    saveScreenshot(model, resultsScreenshotPath + screenshotFile)

    # ------ Add Step Entry
    addStepEntry(resDBfileEntry, '6 - Volume mesh processing', True, time.time() - tmpTime, screenshotPath + screenshotFile, 'Volume mesh generated.')
    # write changes to the database file and closes it.
    dbFile = open(resultsDBPath, 'w')
    json.dump(resultsDB, dbFile, indent=4)
    dbFile.close()

    # Print timing
    print "Total time elapsed in script: " + str(time.time() - tinit)

# Only do this, if we do not import as a module
if __name__ == "__main__":
    main()
