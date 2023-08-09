import re, os, sys

#Create an executable build script which builds the library
#given by the "lib" argument.
def create_build_script(absLibDir,buildFlag):
    buildScript = absLibDir + '/build.sh'
    if os.path.isfile(buildScript):
      os.remove(buildScript)
 
    incDir = os.path.join(absLibDir,'include')
    if not os.path.isdir(incDir): os.mkdir(incDir)
   
    fileObj = open(buildScript, 'w')
    fileObj.write('#!/bin/bash\n')
    fileObj.write('#  Dont forget to make this file executable!\n\n')
    fileObj.write('cd ' + 'source\n')
    fileObj.write('make clean\n')
    fileObj.write('make BUILDFLAG=' + buildFlag + '\n') 
    fileObj.write('cd ../..\n')
    fileObj.close()
    os.chmod(buildScript, 0o774)

def libinfo(relLibDir="",absLibDir="",buildFlag="",args="",macros=[]):
    ans = {}
    args = args.lower()
    create_build_script(absLibDir,buildFlag)
    ans = {"REBUILD":1, "INTERNAL":""}
    return ans
