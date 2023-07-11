import os
import sys
import json

"""
def get_JSON(json_file):
    with open(json_file, "r") as in_d:
        json_data = json.load(in_d)
    return json_data
#end method

"""

def read_json(filename):

    if os.stat(filename).st_size == 0: 
        print("# -- Error -- file is empty:", filename)
        return []
    #end if

    with open(filename, "r") as fh:
        json_data = json.load(fh)
    fh.close()

    return json_data

#end method