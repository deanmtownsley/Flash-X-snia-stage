"""
A simple parser to split test.info into YAML files for individual
simulation directories and test.suite file for Private-Testing repo
"""
import FlashXTest
import yaml


def createSuiteDict(infoNode, suiteDict):
    """
    Create a list of node paths by recursively searching
    till the end of the tree

    Arguments
    ---------
    infoNode : FlashTest node object
    jobList  : Empty jobList
    """
    if infoNode.subNodes:
        for infoNode in infoNode.subNodes:
            createSuiteDict(infoNode, suiteDict)
    else:

        if infoNode.getPathBelowRoot() in suiteDict.keys():
            raise ValueError(f"Duplicate key {infoNode.getPathBelowRoot()}")

        else:

            xmlDict = {}
            for xmlEntry in infoNode.text:
                key, value = xmlEntry.split(":")
                xmlDict.update({key: value.strip()})

            xmlDict["nodeName"] = infoNode.getPathBelowRoot().replace("gce/", "")

            infoDict = {xmlDict["nodeName"]: xmlDict}
            suiteDict.update(infoDict)


if __name__ == "__main__":

    # Parse test.info from the working directory
    # using xml parser from FlashTest backend
    infoNode = FlashXTest.backend.FlashTest.lib.xmlNode.parseXml("test.info")

    # Create a dictionary for test nodes
    # using the recursive function createSuiteDict
    suiteDict = {}
    createSuiteDict(infoNode, suiteDict)

    # Create a list of setupNames using information
    # from suiteDict
    setupList = []
    for key in suiteDict.keys():
        setupList.append(suiteDict[key]["setupName"])

    # Remove duplicates from the list
    setupList = [*set(setupList)]
    print(f"set of setupList: {setupList}")

    # Create list of keys from existing test.info and
    # designated where they should belong
    # in "tests.yaml" or test suite
    yamlkeys = ["setupOptions", "parfiles", "restartParfiles", "transfers"]
    suitekeys = ["numProcs"]

    # Open *.suite file and start populating test suite
    with open("test.suite", "w") as suitefile:

        suitefile.write("# Test suite file for test gce\n")
        suitefile.write('# comments start with "#"\n')

        # Loop over setupList
        for setupName in setupList:

            suitefile.write(f"\n# {setupName}\n")

            # Open a tests.yaml file in setupName/tests folder and
            # start adding tests
            with open(f"SimulationMain/{setupName}/tests/tests.yaml", "w") as yamlfile:

                # Get list of test nodes matching
                # the setup name
                nodeList = [
                    node
                    for node in suiteDict.keys()
                    if suiteDict[node]["setupName"] == setupName
                ]

                yamlfile.write("# YAML file for test information\n")
                yamlfile.write('# comments start with "#"\n')

                # Loop over nodeList and start populating
                for nodeName in nodeList:

                    # Write information to yaml file
                    # This logic is a bit hacky to
                    # to manage paths for parfiles and restartParfiles
                    # TODO: ?? How will paths for transfers be handled ??
                    yamlfile.write(f"\n{nodeName}:\n")
                    for key in yamlkeys:
                        if key in suiteDict[nodeName].keys():
                            if key in ["parfiles", "restartParfiles"]:
                                if f"{setupName}/tests/" in suiteDict[nodeName][key]:
                                    suiteDict[nodeName][key] = suiteDict[nodeName][
                                        key
                                    ].replace(
                                        "<pathToSimulations>/" + setupName + "/tests/",
                                        "",
                                    )
                                elif f"<setupName>" in suiteDict[nodeName][key]:
                                    suiteDict[nodeName][key] = suiteDict[nodeName][
                                        key
                                    ].replace(
                                        "<pathToSimulations>/" + "<setupName>/",
                                        "",
                                    )
                                else:
                                    suiteDict[nodeName][key] = suiteDict[nodeName][
                                        key
                                    ].replace(
                                        "<pathToSimulations>/" + setupName + "/", ""
                                    )

                            yamlfile.write(f"  {key}: {suiteDict[nodeName][key]}\n")

                    # Write line to suitefile
                    suitefile.write(
                        f'{setupName} --test="{nodeName}" --nprocs={suiteDict[nodeName]["numProcs"]}\n'
                    )
