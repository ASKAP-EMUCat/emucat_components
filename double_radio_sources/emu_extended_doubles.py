### wrapper for dragn_hunter.py to integrate output into emucat
### needs to:
### 1) filter out components with a wise match
### 2) run dragn_hunter without writing to file
### 3) output minimal votable based on schema at https://aussrc.atlassian.net/browse/EMUCAT-41


import numpy as np, pandas as pd, argparse
from astropy.table import Table, hstack, vstack, unique
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from dragn_hunter import *


##############################################################################
###functions

def filter_matched_sources():
    "filter out components not to use because they've been matched to a host"
    return


def output_data_to_emucat():
    "subset required data from dragn_hunter based on emucat schema and write to votable"
    
    ###columns to output: id, components (as array), position (ra, dec),
    ### wise_id (once host finding done/ NULL if none), separation
    
    return
