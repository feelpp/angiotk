AngioTK Project
===============

This is the AngioTk project: a collaboration within the <a
href=http://icube-vivabrain.unistra.fr/index.php/Presentation>ANR MN
Vivabrain</a>

# The VIVABRAIN project
The Vivabrain project was a french publicly funded research project
involving*

- University of Reims
- University of Grenoble
- University of Strasbourg
- University Paris-Est ESIEE
- Kitware SAS.

The project is no longer funded but development continues. The main software outcome of this project is contained in
this repository. The software is under development. This document describes how to obtain it and compile it.

## Downloading the software


[![Join the chat at https://gitter.im/vivabrain/angiotk](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/vivabrain/angiotk?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


## Mailing lists

 - angiotk@cemosis.fr (discussions)
 - angiotk-commits@cemosis.fr (git commits)

## Build Status

AngioTk uses Travis-CI for continuous integration.
Travis-CI Build Status :

  - master branch : [![Build Status](https://magnum.travis-ci.com/feelpp/angiotk.svg?token=bYtmwPh98RexKnPLongV&branch=master)](https://magnum.travis-ci.com/feelpp/angiotk)

### Dashboard

  - [http://open.cdash.org/index.php?project=AngioTk](http://open.cdash.org/index.php?project=AngioTk)

### Cloning AngioTk

Then issue the following commands in order to clone the angiotk repository in feelpp/ research
```
git clone http://github.com/feelpp/angiotk.git
```

### Compiling AngioTk

You need the following dependencies for AngioTk to build:
* The vmtk library: [http://www.vmtk.org/download/](http://www.vmtk.org/download/)
* (Optional if using vmtk superbuild) ITK library at least at version 4 (http://www.itk.org)
* (Optional if using vmtk superbuild) VTK version 5.x (vmtk is not compatible with VTK 6.x for the moment)
* The Feel++ library: [https://github.com/feelpp/feelpp](https://github.com/feelpp/feelpp)

The following compilation order is strongly encouraged:
* First compile vmtk superbuild from the sources.
```
cd <vtmk_sources>
mkdir vmtk-build
cmake <vmtk_sources>
```
* Set the environment variables for vmtk, VTK, and ITK to be accessible in you shell. For example, you can use the `module` command or the `vmtk_env.sh` script located in `<vmtk_sources>/vmtk-build/Install`
* Build and install Feel++ (See the installation insctructions on the [Feel++ repository](https://github.com/feelpp/feelpp)). It is advised to let cmake detect the VTK library available, which should be the one built with vmtk. Enabling options such as In-Situ visualization will alter the VTK library chosen by cmake, you would end up with VTK 6.x which won't be compatible with AngioTk.
* Finally build AngioTk:
```
# If you didn't install Feel++ in a standard system installation dir
# you can help cmake in finding the library with the following command:
# export FEELPP_DIR=<feelpp_install_dir>
cd <angiotk_build_dir>
# Configure AngioTK with the modules you need
cmake <angiotk_sources> -DBUILD_MODULE_Meshing=ON
```
