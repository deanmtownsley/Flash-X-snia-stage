"""
A simple parser to split test.info into YAML files for individual
simulation directories
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

    # Parse test.info
    infoNode = FlashXTest.backend.FlashTest.lib.xmlNode.parseXml("test.info")

    # Get suiteDict
    suiteDict = {}
    createSuiteDict(infoNode, suiteDict)

    setupList = []
    for key in suiteDict.keys():
        setupList.append(suiteDict[key]["setupName"])

    setupList = [*set(setupList)]

    yamlkeys = ["setupOptions", "parfiles", "restartParfiles", "transfers"]
    suitekeys = ["numProcs"]

    with open("test.suite", "w") as suitefile:

        suitefile.write("# Test suite file for test gce\n")
        suitefile.write('# comments start with "#"\n')

        for setupName in setupList:

            suitefile.write(f"\n# {setupName}\n")

            with open(f"SimulationMain/{setupName}/tests/tests.yaml", "w") as yamlfile:

                nodeList = [
                    node
                    for node in suiteDict.keys()
                    if suiteDict[node]["setupName"] == setupName
                ]

                yamlfile.write("# YAML file for test information\n")
                yamlfile.write('# comments start with "#"\n')

                for nodeName in nodeList:

                    yamlfile.write(f"\n{nodeName}:\n")
                    for key in yamlkeys:
                        if key in suiteDict[nodeName].keys():
                            if key == "parfiles":
                                if "tests" in suiteDict[nodeName][key]:
                                    suiteDict[nodeName][key] = suiteDict[nodeName][
                                        key
                                    ].replace(
                                        "<pathToSimulations>/" + setupName + "/tests/",
                                        "",
                                    )
                                else:
                                    suiteDict[nodeName][key] = suiteDict[nodeName][
                                        key
                                    ].replace(
                                        "<pathToSimulations>/" + setupName + "/", ""
                                    )

                            yamlfile.write(f"  {key}: {suiteDict[nodeName][key]}\n")

                    suitefile.write(
                        f'{setupName} --test={nodeName} --nprocs={suiteDict[nodeName]["numProcs"]}\n'
                    )
