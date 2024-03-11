# from unitUtils import *
from macroProcessor import macroProcessor, makeVariantName
import subprocess
import shutil
import os

fortran_exts = [".F90", ".f90", ".F", ".f"]


def removeSuffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[: -len(suffix)]
    return input_string


def formatOutput(outpath):
    _, ext = os.path.splitext(outpath)
    if ext in fortran_exts:
        if hasattr(shutil, "which"):
            if shutil.which("fprettify"):
                subprocess.run("fprettify -w1 {}".format(outpath), shell=True, check=True)


# unitDir: path to unit directory with mc files
# objDir: path to object directory
# defsList: list of common defs and unit defs
# varList: list of variant names (unit has variants/varName.ini files)
def generateVariants(unitDir, objDir, defsList, varList, macroOnly=False):
    print("Generating variants {} for unit: {}".format(varList, unitDir))
    mcList = []
    mcListNoVariants = []
    baseList = []
    for f in os.listdir(unitDir):
        if os.path.splitext(f)[-1][-3:] == "-mc":
            mcPath = os.path.join(unitDir, f)
            with open(mcPath) as mcFile:
                lines = mcFile.read()
                if "!!NOVARIANTS" in lines:
                    mcListNoVariants.append(mcPath)
                else:
                    mcList.append(mcPath)
                    base, ext = os.path.splitext(f)
                    baseList.append(base + removeSuffix(ext, "-mc"))
    for var in varList:
        m = macroProcessor()
        if var != "":
            varDir = os.path.join(unitDir, var)
            if os.path.isdir(varDir):
                varDefs = [
                    os.path.join(varDir, f)
                    for f in os.listdir(varDir)
                    if ((os.path.splitext(f)[-1] == ".ini"))
                ]
                m.loadDefsList(defsList[0] + varDefs + defsList[1])
            else:
                m.loadDefsList(defsList[0] + defsList[1])
        else:
            m.loadDefsList(defsList[0] + defsList[1])
        for f in mcList:
            filebase, ext = os.path.splitext(os.path.basename(f))
            outfile = makeVariantName(filebase, var, removeSuffix(ext, "-mc"))
            outpath = os.path.join(objDir, outfile)
            processMcFile(m, f, outpath, macroOnly)

    # convert files with no variants
    m = macroProcessor()
    m.loadDefsList(defsList[0] + defsList[1])
    for f in mcListNoVariants:
        filebase, ext = os.path.splitext(os.path.basename(f))
        outfile = makeVariantName(filebase, "", removeSuffix(ext, "-mc"))
        outpath = os.path.join(objDir, outfile)
        processMcFile(m, f, outpath, macroOnly)

    if "null" in [v.lower() for v in varList]:
        baseList = []
    return baseList


# Checks specific conditions for processing a file.
# mp: Macroprocessor
# mc: file containing macros to process
# out: Output file
# macroOnly: flag for only processing files containing macros. 
def processMcFile(mp, mc, out, macroOnly):
    if os.path.islink(out):
        os.unlink(out)

    # mcIsNewer = not os.path.isfile(out) #or \
        #(os.path.isfile(out) and os.path.getmtime(mc) > os.path.getmtime(out))
    # macro file is newer or the flag is not set.
    # shouldProcess = (not macroOnly) or (macroOnly and mcIsNewer)

    # if shouldProcess:
    # print(f"Updating {mc} file")
    mp.convertFile(mc, out)
    formatOutput(out)


# unitDir: path to unit directory with mc files
# varList: list of variant names (unit has variants/varName.ini files)
def modifyMakefile(unitDir, makefile, varList):
    mcList = []
    mcListNoVariants = []
    for f in os.listdir(unitDir):
        if os.path.splitext(f)[-1][-3:] == "-mc":
            mcPath = os.path.join(unitDir, f)
            with open(mcPath) as mcFile:
                lines = mcFile.read()
                if "!!NOVARIANTS" in lines:
                    mcListNoVariants.append(mcPath)
                else:
                    mcList.append(mcPath)
    for f in mcList:
        filename = os.path.basename(f)
        filebase, ext = os.path.splitext(filename)
        baseObj = filebase + ".o"
        varObjs = []
        varFiles = []
        for var in varList:
            varObjs.append(makeVariantName(filebase, var, ".o"))
            varFiles.append(makeVariantName(filebase, var, removeSuffix(ext, "-mc")))

        with open(makefile, "r") as f:
            lines = f.read()
        lines = lines.replace(baseObj, " ".join(varObjs))
        lines = lines.replace(removeSuffix(filename, "-mc"), " ".join(varFiles))
        with open(makefile, "w") as f:
            f.write(lines)

    for f in mcListNoVariants:
        filename = os.path.basename(f)
        filebase, ext = os.path.splitext(os.path.basename(f))
        baseObj = filebase + ".o"
        varObjs = []
        varFiles = []
        varObjs.append(makeVariantName(filebase, "", ".o"))
        varFiles.append(makeVariantName(filebase, "", removeSuffix(ext, "-mc")))

        with open(makefile, "r") as f:
            lines = f.read()
        lines = lines.replace(baseObj, " ".join(varObjs))
        lines = lines.replace(removeSuffix(filename, "-mc"), " ".join(varFiles))
        with open(makefile, "w") as f:
            f.write(lines)


# DEPRECATED
# Add definitions from all files in passed list, then
# perform macro preprocessing on all files in current dir.
# def processFilesInCurrentDir( macroProc = None, defsList=[] ):
#    if macroProc is None:
#        m = macroProcessor()
#        for defs in defsList:
#            m.loadDefs(defsList)
#    else:
#        m = macroProc
#
#    # Get all .ini files from current directory
#    defs_in_dir = [f for f in os.listdir('.') if (os.path.isfile(f) and (os.path.splitext(f)[-1]==".ini"))]
#    # The following command will trigger a NameError exception if 'import ConfigParser' has failed:
#    m.loadDefs(defs_in_dir)
#
#    # Search current directory and convert all .F90-mc files
#    files = [f for f in os.listdir('.') if (os.path.isfile(f) and (".F90-mc" in f) and os.path.splitext(f)[-1]==".F90-mc")]
#    for f in files:
#        outfile = f.replace(".F90-mc",".F90")
#        print("Macro processor running on "+f)
#        m.convertFile(f,outfile)
