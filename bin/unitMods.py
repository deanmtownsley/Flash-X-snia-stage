"""Module to call unit specific modules from bin/setup.py"""

import os, sys
import glob
import importlib
import inspect
from globals import SetupError, GVarsClass


def main(GVars, unitList):
    """
    Call external modules files present in the source directory
    """
    # Set absolute path for object directory
    objDir = os.path.join(GVars.flashHomeDir, GVars.objectDir)

    # Insert object directory path to sys path for loading python modules
    sys.path.insert(0, objDir)

    # Get a list of python modules from the object directory and strip
    # the directory names and suffix to get modules names
    pyModules = [module.replace(objDir + os.sep,"").replace(".py","")
                 for module in glob.glob(os.path.join(objDir, "*.py"))]

    # Build a list of units that have been included with their UnitMain
    # implementation. This allows for more control over which python
    # modules are loaded and executed.
    mainUnits = []
    for unit in unitList.getList():
        if "Main" in unit.split("/")[-1]:
            mainUnits.append(unit.split("/")[-1].replace("Main", "_setup"))

    # Loop over pyModules and check if they are present in mainUnits
    # If present, then import the modules and call their main function
    for unit in pyModules:
        if unit in mainUnits:
            module = importlib.import_module(unit)

            if len(inspect.signature(module.main).parameters) == 1:
                if (list(inspect.signature(module.main).parameters.values())[0].annotation == GVarsClass):
                    module.main(GVars)
                else:
                    raise SetupError(f"{unit}.main should only have one argument of type {GVarsClass}")
            else:
                raise SetupError(f"{unit}.main should only have one argument of type {GVarsClass}")
