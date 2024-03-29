**EMU extended doubles finding**

contact yjan.gordon@umanitoba.ca

Code for finding extended double radio sources in EMU. This is based on the code dragn_hunter.py (readme below) being developed for high resolution survey data, and this code has been specifically adapted for EMUcat.
Run on the command line as:

*python3 emu_doubles.py radio_catalogue --known_hosts=xid_filename*

where *radio_catalogue* is the file name of the SER component catalogue, and *xid_filename* is the filename of the SER simple cross ID file (e.g. nearest WISE match). The config file *config.txt* has been setup for EMUcat. The output is a votable file of the same name as the *radio_catalogue* file name with *_pairs* added before the extension, this is written to the same directory as the *radio_catalogue*. The output file currently contains the following columns

parameter | description [units]
----------|------------
id_1 | the identifier for the first component of the pair
id_2 | the identifier for the second component of the pair
pair_name | identifier for the pair of the format Jhhmmss.ss±ddmmss.s
cenRA | RA of the pair's flux weighted centroid [deg]
cenDEC | DEC of the pair's flux weighted centroid [deg]

note: for the columns id_1 and id_2, these column names are depended on the column name of the component identifier in the input data. For instance if the column 'component_id' is used to identify the components the output data will contain the columns 'component_id_1' and 'component_id_2'.




Additional information on the general dragn_hunter.py functionality is below (this is retained within the emu_doubles.py script).

**dragn_hunter.py readme**

command line code for searching for pairs of radio components that satisfy specified criteria, and flagging likely Double Radio sources associated with Active Galactic Nuclei (DRAGNs). Run as

*python3 dragn_hunter.py radio_catalogue*

where *radio_catalogue* is the filename of the radio catalogue data file. By default this code assumes that a config file named *config.txt* exists but a different config file can be called by *--config=filename*.
The output is a file with the pair data and (to be added) flagged likely DRAGNs. If importing the functionality from another code or using iPython in the command line, and outputting a file is not desired, then the option *--writepairs='False'* can be used to return an astropy table of the data instead.


**Config File Description**

parameter | description
----------|------------
name_col | name of the object name/id column in the radio catalogue file
ra_col | name of the RA column in the radio catalogue file (assumed to be in decimal degrees)
dec_col | name of the Dec column in the radio catalogue file (assumed to be in decimal degrees)
Speak_col | name of the peak flux column in the radio catalogue file
Stot_col | name of the total flux column in the radio catalogue file
maj_col | name of the major axis size column in the radio catalogue file
min_col | name of the minor axis size column in the radio catalogue file
pa_col | name of the position angle column in the radio catalogue file (assumed to be degrees)
min_brightness | minimum peak flux of radio sources to use in pair finding
min_size | minimum size of radio sources to use in pair finding
data_format | format of radio catalogue file


**Pair Data Description**

the output data is a list of pairs where each **pair** appears only once - i.e. component_1 of a pair does not appear later on as component_2 of the pair. For each component of the pair, the columns specified in the config file will be returned in the pair data with _1 or _2 appended to the column name. In addition several derived columns are given and are described below.

parameter | description [units]
----------|------------
Sep_pair | the angular separation of the pair from the catalogue positions of the constituent components [arcsec]
PA_pair | the position angle East of North from component_1 to component_2 in the pair [deg]
dPA_n | difference in position angle of component_n to PA_pair [deg]
abs_dPA | magnitude of dPA_n
Tflux_ratio | ratio of total flux of component_1 to component_2
Pflux_ratio | ratio of peak flux of component_1 to component_2
medianRA | median (not flux weighted) RA of the pair [deg]
medianDEC | median (not flux weighted) DEC of the pair [deg]
cenRA | RA of the pair's flux weighted centroid [deg]
cenDEC | DEC of the pair's flux weighted centroid [deg]
d_2NN | distance to the second nearest neighbour of the pair [arcsec]
pair_name | identifier for the pair of the format Jhhmmss.ss±ddmmss.s


**Code Dependencies**
The following python packages are required to run this code (version used in development):
* numpy (1.17.4)
* pandas (1.0.3)
* argparse (1.1)
* astropy (3.2)
