from unitUtils import *
from macroProcessor import *

# unitDir: path to unit directory with mc files
# objDir: path to object directory
# defsList: list of common defs and unit defs
# varList: list of variant names (unit has variants/varName.ini files)
def generateVariants(unitDir, objDir, defsList, varList):
    print("Generating variants {} for unit: {}".format(varList,unitDir))
    mcList = []
    mcListNoVariants = []
    baseList = []
    for f in os.listdir(unitDir):
        if os.path.splitext(f)[-1]==".F90-mc" :
            mcPath = os.path.join(unitDir,f)
            with open(mcPath) as mcFile:
                lines = mcFile.read()
                if '!!NOVARIANTS' in lines:
                    mcListNoVariants.append(mcPath)
                else:
                    mcList.append(mcPath)
                    baseList.append( os.path.splitext(f)[0]+".F90" )
    for var in varList:
        m = macroProcessor()
        if(var != ''):
            varDir = os.path.join(unitDir,var)
            if(os.path.isdir(varDir)):
                varDefs = [os.path.join(varDir,f) for f in os.listdir(varDir) if ((os.path.splitext(f)[-1]==".ini"))]
                m.loadDefsList(defsList[0] + varDefs + defsList[1])
            else:
                m.loadDefsList(defsList[0] + defsList[1])
        else:
            m.loadDefsList(defsList[0] + defsList[1])
        for f in mcList:
            filebase,_ = os.path.splitext( os.path.basename(f) )
            outfile = makeVariantName(filebase,var,"F90")
            outpath = os.path.join(objDir, outfile)
            if(os.path.islink(outpath)):
                os.unlink(outpath)

            m.convertFile(f,outpath)

    #convert files with no variants
    m = macroProcessor()
    m.loadDefsList(defsList[0] + defsList[1])
    for f in mcListNoVariants:
        filebase,_ = os.path.splitext( os.path.basename(f) )
        outfile = makeVariantName(filebase,'',"F90")
        outpath = os.path.join(objDir, outfile)
        if(os.path.islink(outpath)):
            os.unlink(outpath)

        m.convertFile(f,outpath)

    if 'null' in [ v.lower() for v in varList]:
        baseList = []
    return baseList

# unitDir: path to unit directory with mc files
# varList: list of variant names (unit has variants/varName.ini files)
def modifyMakefile(unitDir, makefile, varList):
    mcList = []
    mcListNoVariants = []
    for f in os.listdir(unitDir):
        if os.path.splitext(f)[-1]==".F90-mc" :
            mcPath = os.path.join(unitDir,f)
            with open(mcPath) as mcFile:
                lines = mcFile.read()
                if '!!Novariants' in lines:
                    mcListNoVariants.append(mcPath)
                else:
                    mcList.append(mcPath)
    for f in mcList:
        filebase,_ = os.path.splitext( os.path.basename(f) )
        baseObj = filebase + '.o'
        varObjs = []
        for var in varList:
            varObjs.append( makeVariantName(filebase,var,'o') )

        with open(makefile,'r') as f:
            lines = f.read()
        lines = lines.replace(baseObj,' '.join(varObjs))
        with open(makefile,'w') as f:
            f.write(lines)

    for f in mcListNoVariants:
        filebase,_ = os.path.splitext( os.path.basename(f) )
        baseObj = filebase + '.o'
        varObjs = []
        varObjs.append( makeVariantName(filebase,'','o') )

        with open(makefile,'r') as f:
            lines = f.read()
        lines = lines.replace(baseObj,' '.join(varObjs))
        with open(makefile,'w') as f:
            f.write(lines)


# Add definitions from all files in passed list, then
# perform macro preprocessing on all files in current dir.
def processFilesInCurrentDir( macroProc = None, defsList=[] ):
    if macroProc is None:
        m = macroProcessor()
        for defs in defsList:
            m.loadDefs(defsList)
    else:
        m = macroProc

    # Get all .ini files from current directory
    defs_in_dir = [f for f in os.listdir('.') if (os.path.isfile(f) and (os.path.splitext(f)[-1]==".ini"))]
    # The following command will trigger a NameError exception if 'import ConfigParser' has failed:
    m.loadDefs(defs_in_dir)

    # Search current directory and convert all .F90-mc files
    files = [f for f in os.listdir('.') if (os.path.isfile(f) and (".F90-mc" in f) and os.path.splitext(f)[-1]==".F90-mc")]
    for f in files:
        outfile = f.replace(".F90-mc",".F90")
        print("Macro processor running on "+f)
        m.convertFile(f,outfile)


