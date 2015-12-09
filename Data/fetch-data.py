#!/usr/bin/env python

import pydas
import os
import sys
import shutil

def main():
    serverURL = "http://vivabrain.u-strasbg.fr/midas"

    print "Pydas AngioTK data downloader"
    print "Usage: " + sys.argv[0] + " <optinal:data_url_in_AngioTK_Community> ..."
    print "Example: " + sys.argv[0] + ""
    print " -> will download all the public and private datasets"
    print "Example: " + sys.argv[0] + " Public"
    print " -> will download all the public datasets"
    print ""

    # Build an array of data to download
    # defaults to public and private data
    baseURL = "/communities/AngioTK"
    dataURL = [ "Public", "Private" ]

    # If we have arguments, we replace the default directories to be downloaded
    if(len(sys.argv) > 1):
        dataURL = []
        for i in range(1, len(sys.argv)):
            dataURL.append(sys.argv[i])

    print "The following data will be downloaded:"
    for i in range(len(dataURL)):
        print "- " + baseURL + "/" + dataURL[i]

    print ""
    print "Using pydas " + pydas.__version__
    print "Using server URL: " + serverURL
    core_driver = pydas.drivers.CoreDriver(serverURL)
    print "Server version: " + core_driver.get_server_version()

    token = pydas.login(url=serverURL)
    print "API token = \"" + token + "\""

    # Eventually download data
    for i in range(len(dataURL)):
        remoteURL = baseURL + "/" + dataURL[i]
        localPath = "." #os.path.dirname("." + baseURL + "/" + dataURL[i])
        if(not os.path.exists(localPath)):
            os.makedirs(localPath)
        #else:
            #shutil.rmtree(localPath)
            #os.makedirs(localPath)
        pydas.download(remoteURL, local_path=localPath)

if __name__ == "__main__":
    main()
