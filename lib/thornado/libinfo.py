import os, globals
from globals import *

#Create an executable build script which builds the thornado library.
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

    USE_GPU = "FALSE"
    USE_CUDA = "FALSE"
    USE_CUBLAS = "FALSE"
    USE_HIP = "FALSE"
    USE_ROCM = "FALSE"
    USE_ONEMKL = "FALSE"
    USE_OACC = "FALSE"
    USE_OMP_OL = "FALSE"
    USE_OMP = "FALSE"

    momentClosure = "MINERBO"
    thornadoSolver = "EMAB"
    thornadoOrder = "ORDER_V"
    microphysics = ""

    if "thornadoGPU" in setupVars:
        thornadoGPU = str(setupVars["thornadoGPU"]).upper()
        if thornadoGPU in ["NVIDIA","AMD","INTEL"]:
            USE_GPU = "TRUE"
            if thornadoGPU == "AMD":
                USE_HIP = "TRUE"
                USE_ROCM = "TRUE"
            elif thornadoGPU == "INTEL":
                USE_ONEMKL = "TRUE"
            elif thornadoGPU == "NVIDIA":
                USE_CUDA = "TRUE"
                USE_CUBLAS = "TRUE"

    if "thornadoACC" in setupVars:
        USE_OACC = str(setupVars["thornadoACC"]).upper()

    if "thornadoOMP_OL" in setupVars:
        USE_OMP_OL = str(setupVars["thornadoOMP_OL"]).upper()

    if "thornadoOMP" in setupVars:
        USE_OMP = str(setupVars["thornadoOMP"]).upper()

    if "momentClosure" in setupVars:
        momentClosure = str(setupVars["momentClosure"]).upper()

    if "thornadoSolver" in setupVars:
        thornadoSolver = str(setupVars["thornadoSolver"]).upper()

    if "thornadoOrder" in setupVars:
        thornadoOrder = str(setupVars["thornadoOrder"]).upper()

    for unit in GVars.withUnits:
        if "WEAKLIB" == os.path.basename(unit).upper():
            microphysics = "WEAKLIB"

    fd = os.open(buildScript, os.O_WRONLY|os.O_CREAT|os.O_EXCL, 0x1e4) # 0x1e4 == 0o744
    fileObj = os.fdopen(fd, 'w')
    fileObj.write('#!/bin/sh\n')
    fileObj.write('#  This file should be executable!\n\n')
    fileObj.write('set -ex\n')  # set -e to fail when an error occurs, -x to trace commands
    fileObj.write('cd source/SandBox/Interface_FLASH\n')
    fileObj.write('make -f Makefile.Flash clean\n')
    fileObj.write('make -f Makefile.Flash -j8' +
                  ' BUILDFLAG=' + buildFlag +
                  ' USE_GPU=' + USE_GPU +
                  ' USE_CUDA=' + USE_CUDA +
                  ' USE_CUBLAS=' + USE_CUBLAS +
                  ' USE_HIP=' + USE_HIP +
                  ' USE_ROCM=' + USE_ROCM +
                  ' USE_ONEMKL=' + USE_ONEMKL +
                  ' USE_OACC=' + USE_OACC +
                  ' USE_OMP_OL=' + USE_OMP_OL +
                  ' USE_OMP=' + USE_OMP +
                  ' MOMENT_CLOSURE=' + momentClosure +
                  ' NEUTRINO_MATTER_SOLVER=' + thornadoSolver +
                  ' TWOMOMENT_ORDER=' + thornadoOrder +
                  ' MICROPHYSICS=' + microphysics +
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
    #thornado implementations that exist in different directories.
    return {"REBUILD":(not GVars.noClobber), "INTERNAL":""}
