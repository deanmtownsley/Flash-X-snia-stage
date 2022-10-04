"""
    Parser library for log files and preformance data
"""

import os


class LogDict:
    """
    Class to handle dictionary of logfile entries
    """

    def __init__(self, log_file):
        """
        Initialization method

        Argument
        --------
        log_file : full path to `.log` file
        """
        # Create an empty dictionary
        # object for logfile entries
        self.log_dict = {}

        # Call internal method to parse
        # logfile and populate dictionary
        self._set_log_dict(log_file)

    def __getitem__(self, timer_key):
        """
        Internal getter, this done to enforce that
        contents of dictionary are read-only

        Arguments
        ---------
        timer_key : key for the timer entry
        """
        return self.log_dict[timer_key]

    def __repr__(self):
        """
        Internal representation method to display
        dictionary contents
        """
        return f"{self.log_dict}"

    @staticmethod
    def _file_to_list(file_path):
        """
        Static method to convert a file to a list of lines

        Arguments
        ---------
        file_path: string (name of file - full path)

        Returns
        -------
        file_list: list of lines
        """
        # Create an empty list
        # to populate as the file is passed
        file_list = []

        # Open the input file in read-only mode
        with open(file_path, "r") as working_file:

            # loop over lines
            # in working file
            for line in working_file:

                # append to
                # file list
                file_list.append(line.strip())

        return file_list

    def _set_log_dict(self, log_file):
        """
        Internal method to set logfile dictionary

        Arguments
        ---------
        log_file : path to logfile
        """
        # convert file to a list
        log_list = self.__class__._file_to_list(log_file)

        # Create empty lists for
        # parsing data from `log_list`
        #
        # `accounting_index` : list to store index of desired line
        #
        # `dashed_index` : list to store indices of ------- lines
        #
        # `equal_index`  : list to store indices of ======= lines
        accounting_index = []
        dashed_index = []
        equal_index = []

        # list of reference keys to probe
        accounting_list = ["max/proc", "min/proc", "avg/proc", "num calls"]

        # enumerate line and index and get locations
        # for markers defining the boundaries of
        # datasets that need to extracted
        for index, line in enumerate(log_list):
            if all(keyword in line for keyword in accounting_list):
                accounting_index.append(index)
            if "-----------------------------" in line:
                dashed_index.append(index)
            if "==============================" in line:
                equal_index.append(index)

        # create a bound list to store
        # bounds of desired datasets
        bound_list = []

        # loop over accounting_index and append
        # bounds
        for index in accounting_index:
            bound_list.append(
                [
                    min([i for i in dashed_index if index < i]) + 1,
                    min([i for i in equal_index if index < i]) - 1,
                ]
            )

        # Check length of the line index
        # this to enforce only one entry exists
        # for the desired fields
        if len(bound_list) > 1:
            raise ValueError("[flashparser] Unrecognized Logfile")

        # loop over bounds and get desired data
        for bound in bound_list:

            # extract key and values
            for index in range(bound[0], bound[1] + 1):

                # extract and update accounting_dict
                timer_key = " ".join(log_list[index].split()[:-4])

                # Create a value list in proper data format
                value_list = [float(value) for value in log_list[index].split()[-4:-1]]
                value_list.append(int(log_list[index].split()[-1]))

                # Map value_list to a dictionary
                value_dict = {
                    key: value for key, value in zip(accounting_list, value_list)
                }

                self.log_dict.update({timer_key: value_dict})

        # Deal with `.log.csv` containing
        # data from multiple processors.
        # Start with an empty list to deal with exceptions
        csv_list = []

        # check if csv path exists and
        # convert file to a list
        if os.path.exists(log_file + ".csv"):
            csv_list = self.__class__._file_to_list(log_file + ".csv")

        # iterate over each entry from csv_list
        for csv_entry in csv_list:

            # extract timer key
            timer_key = csv_entry.split(",")[0]

            # extract timer level
            timer_level = int(csv_entry.split(",")[2])

            # extract values
            value_list = [float(value) for value in csv_entry.split(",")[3:]]

            # Map value_list to a dictionary
            value_dict = {
                "proc/" + str(proc).zfill(4): value
                for proc, value in enumerate(value_list)
            }

            # Extract level of the timer
            value_dict["level"] = timer_level

            # Update log dicitonary
            self.log_dict[timer_key].update(value_dict)
