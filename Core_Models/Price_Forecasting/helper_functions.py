# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import os

def generate_path(directory_list):
    dir = os.getcwd()
    dir = os.path.join(dir, "Core_Models", "HydroBoost")

    for d in directory_list:
        dir = os.path.join(dir, d)
        if os.path.exists(dir) == False:
            os.mkdir(dir)
