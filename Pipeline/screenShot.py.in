#!/usr/bin/env python

# This scripts takes an offline screenshot using ParaView
# usage: pvpython <script_name> <path_to_data_file>

# To use this, you can either use:
# - The classical ParaView version with X exported through ssh (DISPLAY)
# - The EGL version of ParaView (no need for X or DISPLAY)

import sys
from paraview.simple import *

if( len(sys.argv) != 2 ):
    print "Usage: pvpython " + sys.argv[0] + " <data_file>"
    exit(1)

reader = OpenDataFile(sys.argv[1])
reader.UpdatePipeline()

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [ 1920, 1080 ]

Show(reader, view)

displayData = Show(reader, view)

SaveScreenshot("img.png", view)
