#!/usr/local/bin python
import os, re
import argparse
try:
  from configparser import ConfigParser    #Python 3
except ImportError: pass

#standardMacroDefLibrary = ['macroDefsTemplate.ini']
string_char = "\""
comment_char = "!"
macro_keyword = 'M'
macro_regex = r'\s*@\s*' + re.escape(macro_keyword) + r'\s*\w*(?:\s*\(\s*[\w.]*\s*(?:,\s*[\w.]*\s*)*\))?'
invocation_regex = r'(?P<indent>\s*)@\s*' + re.escape(macro_keyword) + \
                   r'\s*(?P<key>\w*)(?:\s*\((?P<arglist>\s*[\w.]*\s*(?:,\s*[\w.]*\s*)*)\))?'

class macroProcessor:
  def __init__(self):
    self.mdict = {}
    self.argdict = {}
    self.typedict = {}
    self.indentdict = {}
    self.keylist = []
    #self.loadDefsList(standardMacroDefLibrary)

  ######### LOADING DEFS ###########
  # Read a file and add contained macro definitions
  # to the macro dictionary.
  def loadDefs(self,filename):
    parser = ConfigParser(comment_prefixes=('#!','#!!'))
    parser.optionxform = lambda option: option # read keys case sensitively
    if filename is not None:
      parser.read(filename)
      for section in parser.sections():
        self.keylist.append(section)

        self.mdict[section] = parser.get(section,'definition').strip()

        try:
            args = parser.get(section,'args')
            self.argdict[section] = [item.strip() for item in args.strip().split(',')]
        except:
            self.argdict[section] = []

        try:
            self.typedict[section] = parser.get(section,'type')
        except:
            self.typedict[section] = ''
        try:
            indents = parser.get(section,'line_indents').strip().split(',')
            self.indentdict[section] = [int(i) for i in indents]
        except:
            self.indentdict[section] = [0]

  # Load a whole list of files.
  def loadDefsList(self,filenames):
    for filename in filenames:
      self.loadDefs(filename)

  ####### EXPAND INVOCATIONS OF MACROS #############

  def getMacroDef(self,key,args=[],indent=''):
    # get definition and replace args
    definition =  self.mdict[key]
    arglist = self.argdict[key]
    for i,arg in enumerate(args):
      if(i<len(arglist) ):
        arg_re = r'\b'+arglist[i]+r'\b' #don't substitute substrings
        definition = re.sub( arg_re, arg, definition)

    # add appropriate indent to all lines
    def_lines = definition.split('\n')
    if (len(def_lines[0].strip()) > 0):
      if (def_lines[0].strip()[0] == '#'):
        #ensure preprocessor directives are not inserted inline
        def_lines[0] = '\n' + def_lines[0]
    j = 0 # track place in per-line indent list
    n = len(self.indentdict[key])
    for i in range(len(def_lines)):
      lineStripped = def_lines[i].strip()
      lead = indent + ' '*self.indentdict[key][j]
      if (len(lineStripped) > 0):
        if (lineStripped[0] == '#'):
          lead = ''
      def_lines[i] = lead + def_lines[i]
      if(j<n-1): j = j+1
    definition = '\n'.join(def_lines)

    return definition

  def expandMacro(self, invocation, macroStack):
    expansion = invocation
    keymatch = False

    # use regex to get parts of invocation
    invocation_parts = re.match(invocation_regex, invocation)
    macroName = invocation_parts.group('key')
    indent = invocation_parts.group('indent')

    # check to make sure recursive calls don't cause infinite loops
    if macroName in macroStack:
      mlist = ', '.join(macroStack)
      msg = "Error: macro(s): (%s) causing a loop via macro recursion."%mlist
      raise SyntaxError(msg)

    for key in self.keylist:
      if (macroName == key):
        args=[]
        if len(self.argdict[key])>0:
          argtext = invocation_parts.group('arglist')
          if argtext is None:
            msg = "Error: argument list expected for macro %s"%macroName
            raise SyntaxError(msg)
          args = argtext.split(',')

        expansion = self.getMacroDef(macroName,args,indent)
        keymatch = True
        break

    # if expansion has a recursive macro, process lines again
    if (expansion.count(macro_keyword)>0 and keymatch):
      macroStack.append(macroName)
      expansionLines = expansion.split('\n')
      for i,line in enumerate(expansionLines):
        expansionLines[i] = self.processLine(line,macroStack)
      expansion = '\n'.join(expansionLines)

    return expansion


  ######### PROCESSING FILES  #################

  # Process a single line
  def processLine(self,lineIn,macroStack=None):
    if macroStack is None:
      macroStack = []

    lineOut = lineIn

    #scan line for comments and strip them
    lineStripped = lineIn
    inString = False
    for i,ch in enumerate(lineIn):
        if ch==string_char:
            inString = not inString
        if (ch==comment_char and not inString):
            lineStripped = lineIn[0:i]
            break

    if len(lineStripped.split())>=1:
      # find matches for macro invocation regex
      invocation_list = re.findall(macro_regex,lineStripped)

      # expand each invocation and replace
      for invocation in invocation_list:
        expansion = self.expandMacro(invocation,macroStack)
        lineOut = lineOut.replace(invocation, expansion, 1)

    return lineOut

  # Process a whole file
  def convertFile(self,filename,output):
    with open(output,'w') as f:
      lines = open(filename).readlines()
      for line in lines:
        f.write(self.processLine(line))

###########################################################
def makeVariantName(base, var, ext):
    if(var == '' or var.lower()=='null'):
        outfile = base + "." + ext
    else:
        outfile = base + "_" + var + "." + ext
    return outfile

# sourceDir - absolute path to source directory
# subDir - such that os.path.join(sourceDir,subDir) is absolute path to unit dir
def recursiveGetDefs(sourceDir,subDir,binDir,simDir):
    defsList = []
    unitList = subDir.split("/")
    currentUnit = sourceDir

    # get common defs from bin dir
    for f in os.listdir(binDir):
        fpath = os.path.join(binDir,f)
        if os.path.isfile(fpath) and os.path.splitext(fpath)[-1]==".ini":
            defsList.append(fpath )

    for unit in unitList:
        currentUnit = os.path.join(currentUnit,unit)
        for f in os.listdir(currentUnit):
            fpath = os.path.join(currentUnit,f)
            if os.path.isfile(fpath) and os.path.splitext(fpath)[-1]==".ini":
                defsList.append(fpath )

    # defs in simulation directory should overwrite others
    for f in os.listdir(simDir):
        fpath = os.path.join(simDir,f)
        if os.path.isfile(fpath) and os.path.splitext(fpath)[-1]==".ini":
            defsList.append(fpath )

    return defsList

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
                if '!!Novariants' in lines:
                    mcListNoVariants.append(mcPath)
                else:
                    mcList.append(mcPath)
                    baseList.append( os.path.splitext(f)[0]+".F90" )
    for var in varList:
        m = macroProcessor()
        if(var != ''):
            varDir = os.path.join(unitDir,var)
            varDefs = [os.path.join(varDir,f) for f in os.listdir(varDir) if ((os.path.splitext(f)[-1]==".ini"))]
            m.loadDefsList(defsList + varDefs)
        else:
            m.loadDefsList(defsList)
        for f in mcList:
            filebase,_ = os.path.splitext( os.path.basename(f) )
            outfile = makeVariantName(filebase,var,"F90")
            outpath = os.path.join(objDir, outfile)
            if(os.path.islink(outpath)):
                os.unlink(outpath)

            m.convertFile(f,outpath)

    #convert files with no variants
    m = macroProcessor()
    m.loadDefsList(defsList)
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


# Can run directly to process some file.
def main():
    parser=argparse.ArgumentParser(description='Macro Preprocessor for FLASH')
    parser.add_argument('--filename','-f',type=str,help='file to convert, if none is specified do all in current dir')
    parser.add_argument('--output','-o',type=str,help='output filename')
    parser.add_argument('--macroDefs','-m',type=str,help='file with list of extra macro definitions')
    args=parser.parse_args()

    m = macroProcessor()
    if(args.macroDefs is not None): m.loadDefs(args.macroDefs)

    if args.filename is None:
        processFilesInCurrentDir()

    else:
      if args.macroDefs is None:
        defs_in_dir = [f for f in os.listdir('.') if (os.path.isfile(f) and (".ini" in f))]
        m.loadDefs(defs_in_dir)
      # Convert just the given file
      if args.output is None:
        if (args.filename.count(".F90-mc")>=1):
          args.output = args.filename.replace(".F90-mc",".F90")
        else:
          args.output = args.filename+".out"
      m.convertFile(args.filename,args.output)


if __name__ == '__main__':
    main()
