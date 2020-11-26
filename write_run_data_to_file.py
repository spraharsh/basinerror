"""
Writes run data to file
"""

import yaml
import tempfile
import itertools as IT
import os




def write_run_data_to_file(optimizer_parameters, run_diagnostics, folder_location, name='run_parameters.yaml'):
    """
    Merges the optimizer paraters and run diagnostics and writes to a single file

    Parameters
    ----------
    optimizer_parameters: dict
        Dictionary of optimizer parameters for the run
    finaldata: ret container
        optimizer data that includes (function evaluations, hessian evaluations and steps)
    folder_location: path string
        location to write file to
    name: str
        name of the file
    """
    dictionary  = {
        "brand": "Ford",
        "model": "Mustang",
        "year": 1964}
    

    optimizer_parameters.update(run_diagnostics)

    print(optimizer_parameters)
    full_data = optimizer_parameters
    print(full_data)

    filename = uniquify(folder_location + '/' + name)
    
    with open(folder_location + '/' + name, 'w') as outfile:
        yaml.dump(full_data, outfile)

def uniquify(path, sep = ''):
    """
    makes the file name unique. From 
    https://stackoverflow.com/a/13852851/6788900

    """

    def name_sequence():
        count = IT.count()
        yield ''
        while True:
            yield '{s}{n:d}'.format(s = sep, n = next(count))
    orig = tempfile._name_sequence 
    with tempfile._once_lock:
        tempfile._name_sequence = name_sequence()
        path = os.path.normpath(path)
        dirname, basename = os.path.split(path)
        filename, ext = os.path.splitext(basename)
        fd, filename = tempfile.mkstemp(dir = dirname, prefix = filename, suffix = ext)
        tempfile._name_sequence = orig
    return filename

if __name__ == "__main__":
    # test data
    optimizer_parameters = {'c': 4, 'a': 10, 'b': 8, 'd': 6}
    run_diagnostics = {'nfev':0, 'nhev':0, 'nsev':0}
    location = os.getcwd()
    write_run_data_to_file(optimizer_parameters, run_diagnostics, location, name='write_test.yaml')
    