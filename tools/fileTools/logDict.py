"""
    Parser for log files and preformance data

    Usage (works with log and log.csv files):

    # standard libraries
    import sys
    
    # insert path to flash-x
    sys.path.insert(0, f"/path/to/Flash-X/tools/fileTools")

    # load dictionary class for log files
    from logDict import LogDict

    # load dictionary from log file
    log_dict = LogDict("/path/to/LogFile.log")

    # Get number of calls for a timer
    log_dict["Grid_updateRefinement"]["num calls"]

    # Average time per proc
    log_dict["Grid_updateRefinement"]["avg/proc"]

    # Time per processor
    log_dict["Grid_updateRefinement"]["proc/0000"]
    log_dict["Grid_updateRefinement"]["proc/0001"]
    .
    .
    log_dict["Grid_updateRefinement"]["proc/0010"]
    .
    .

    # Get timer level
    log_dict["Grid_updateRefinement"]["level"]
"""

import os


class LogDict:
    """
    Class to handle dictionary of logfile entries
    """

    def __init__(self, logFile):
        """
        Initialization method

        Argument
        --------
        logFile : full path to `.log` file
        """
        # Create an empty dictionary
        # object for logfile entries
        self._logDict = {}

        # Call internal method to parse
        # logfile and populate dictionary
        self._setLogDict(logFile)

    def __getitem__(self, timerKey):
        """
        Internal getter, this done to enforce that
        contents of dictionary are read-only

        Arguments
        ---------
        timerKey : key for the timer entry
        """
        return self._logDict[timerKey]

    def __repr__(self):
        """
        Internal representation method to display
        dictionary contents
        """
        return f"{self._logDict}"

    @staticmethod
    def _fileToList(filePath):
        """
        Static method to convert a file to a list of lines

        Arguments
        ---------
        filePath: string (name of file - full path)

        Returns
        -------
        fileList: list of lines
        """
        # Create an empty list
        # to populate as the file is passed
        fileList = []

        # Open the input file in read-only mode
        with open(filePath, "r") as workingFile:

            # loop over lines
            # in working file
            for line in workingFile:

                # append to
                # file list
                fileList.append(line.strip())

        return fileList

    def _setLogDict(self, logFile):
        """
        Internal method to set logfile dictionary

        Arguments
        ---------
        logFile : path to logfile
        """
        # convert file to a list
        logList = self.__class__._fileToList(logFile)

        # Create empty lists for
        # parsing data from `logList`
        #
        # `accountingIndex` : list to store index of desired line
        #
        # `dashedIndex` : list to store indices of ------- lines
        #
        # `equalToIndex`  : list to store indices of ======= lines
        accountingIndex = []
        dashedIndex = []
        equalToIndex = []

        # list of reference keys to probe
        accountingList = ["max/proc", "min/proc", "avg/proc", "num calls"]

        # enumerate line and index and get locations
        # for markers defining the boundaries of
        # datasets that need to extracted
        for index, line in enumerate(logList):
            if all(keyword in line for keyword in accountingList):
                accountingIndex.append(index)
            if "-----------------------------" in line:
                dashedIndex.append(index)
            if "==============================" in line:
                equalToIndex.append(index)

        # create a bound list to store
        # bounds of desired datasets
        boundList = []

        # loop over accountingIndex and append
        # bounds
        for index in accountingIndex:
            boundList.append(
                [
                    min([i for i in dashedIndex if index < i]) + 1,
                    min([i for i in equalToIndex if index < i]) - 1,
                ]
            )

        # Check length of the line index
        # this to enforce only one entry exists
        # for the desired fields
        if len(boundList) > 1:
            raise ValueError("[logDict] Unrecognized Logfile")

        # loop over bounds and get desired data
        for bound in boundList:

            # extract key and values
            for index in range(bound[0], bound[1] + 1):

                # extract and update accounting_dict
                timerKey = " ".join(logList[index].split()[:-4])

                # Create a value list in proper data format
                valueList = [float(value) for value in logList[index].split()[-4:-1]]
                valueList.append(int(logList[index].split()[-1]))

                # Map valueList to a dictionary
                valueDict = {
                    key: value for key, value in zip(accountingList, valueList)
                }

                self._logDict.update({timerKey: valueDict})

        # Deal with `.log.csv` containing
        # data from multiple processors.
        # Start with an empty list to deal with exceptions
        csvList = []

        # check if csv path exists and
        # convert file to a list
        if os.path.exists(logFile + ".csv"):
            csvList = self.__class__._fileToList(logFile + ".csv")

        # iterate over each entry from csvList
        for csvEntry in csvList:

            # extract timer key
            timerKey = csvEntry.split(",")[0]

            # extract timer level
            timerLevel = int(csvEntry.split(",")[2])

            # extract values
            valueList = [float(value) for value in csvEntry.split(",")[3:]]

            # Blank dictionary
            valueDict = {}

            # Map valueList to a dictionary
            valueDict["all procs"] = valueList

            # Extract level of the timer
            valueDict["level"] = timerLevel

            # Update log dicitonary
            self._logDict[timerKey].update(valueDict)
