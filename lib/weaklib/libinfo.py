import os, globals
from globals import *

#Create an executable build script which builds the weaklib library.
#The "make clean" ensures we always rebuild the library which is
#necessary to prevent the following example situation: we switch the FLASH
#compiler from absoft to lahey but the library is still built with absoft.
def create_build_script(absLibDir,buildFlag,args):
    buildScript = absLibDir + '/build.sh'
    if os.path.isfile(buildScript):
        os.remove(buildScript)

    incDir = os.path.join(absLibDir,'include')
    if not os.path.isdir(incDir): os.mkdir(incDir)

    setupVars = GVars.setupVars.getdict()

    if "weaklibACC" in setupVars:
        USE_OACC = str(setupVars["weaklibACC"]).upper()
    else:
        USE_OACC = "FALSE"

    if "weaklibOMP_OL" in setupVars:
        USE_OMP_OL = str(setupVars["weaklibOMP_OL"]).upper()
    else:
        USE_OMP_OL = "FALSE"

    if "weaklibOMP" in setupVars:
        USE_OMP = str(setupVars["weaklibOMP"]).upper()
    else:
        USE_OMP = "FALSE"

    fd = os.open(buildScript, os.O_WRONLY|os.O_CREAT|os.O_EXCL, 0x1e4) # 0x1e4 == 0o744
    fileObj = os.fdopen(fd, 'w')
    fileObj.write('#!/bin/sh\n')
    fileObj.write('#  This file should be executable!\n\n')
    fileObj.write('set -ex\n')  # set -e to fail when an error occurs, -x to trace commands
    fileObj.write('cd source/Distributions/ExternalLibrary\n')
    fileObj.write('make -f Makefile.Flash clean\n')
    fileObj.write('make -f Makefile.Flash -j8' +
                  ' BUILDFLAG=' + buildFlag +
                  ' USE_OACC=' + USE_OACC +
                  ' USE_OMP_OL=' + USE_OMP_OL +
                  ' USE_OMP=' + USE_OMP +
                  '\n')
    fileObj.write('cd ../../../\n')
    fileObj.close()


#We use this function to create a custom build file for our internal library.
#This is required because we want to use the compilation flags specified in
#the FFLAGS_[DEBUG|TEST|OPT] variable when building the library. This will
#avoid a prior problem where by the library was not compiled with reals
#promoted to 8 bytes, but the rest of the FLASH source used 8 byte reals.
def libinfo(relLibDir="",absLibDir="",buildFlag="",args="",macros=[]):

    args = args.lower()
    create_build_script(absLibDir,buildFlag,args)

    #Specify that we want to rebuild the library each time we resetup.
    #Also label the library as internal because we provide the source code.
    #Do not provide a path because we do not have several different
    #weaklib implementations that exist in different directories.
    return {"REBUILD":(not GVars.noClobber), "INTERNAL":""}
