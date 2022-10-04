"""Script to parse log file will be converted to a small library"""


def fileToList(file_path):
    """
    Convert a file to a list of lines

    Arguments
    ---------
    file_path: string (name of file - full path)

    Returns
    -------
    file_list: list of lines
    """

    # Empty list for lines
    file_list = []

    # Open the input file in read-only mode
    with open(file_path, "r") as working_file:

        # loop over lines in working file
        for line in working_file:

            # file list
            file_list.append(line.strip())

    return file_list


def getPerformanceDict(file_path):
    """
    get a indices for accounting

    Arguments
    ---------
    file_path : string (name of file - full path)
    """

    # convert file to a list
    log_list = fileToList(file_path)

    # Create empty lists for parsing data
    accounting_index = []
    dashed_index = []
    equal_index = []

    # output dictionary
    accounting_dict = {}

    # key refrence
    accounting_list = ["max/proc", "min/proc", "avg/proc", "num calls"]

    # enumerate line and index and get locations
    # for markers defining the boundaries of
    # datasets that need to extracted
    for index, line in enumerate(log_list):
        if all(keyword in line for keyword in accounting_list):
            accounting_index.append(index)
        if "------------------------------" in line:
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

    if len(bound_list) > 1:
        raise ValueError("Unrecgonized logfile")

    # loop over bounds and get desired data
    for bound in bound_list:

        # extract key and values
        for index in range(bound[0], bound[1] + 1):

            # extract and update accounting_dict
            timer_key = " ".join(log_list[index].split()[:-4])

            # Create a value list in proper data format
            value_list = [float(value) for value in log_list[index].split()[-4:-1]]
            value_list.append(int(log_list[index].split()[-1]))

            # Map list to a dictionary
            value_dict = {key: value for key, value in zip(accounting_list, value_list)}

            accounting_dict.update({timer_key: value_dict})

    return accounting_dict
