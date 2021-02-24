###code for finding Double Radio AGN (DRAGNs) from radio catalogue data
###finds pairs of radio components based on flux and size criteria
###developed by Yjan Gordon (yjan.gordon@umanitba.ca) Jan 2021

import numpy as np, pandas as pd, argparse
from astropy.table import Table, hstack, vstack, unique
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u

##############################################################################
###functions

def source_name(ra, dec, aprec=2, dp=5, prefix=''):
    ###create source name to nearest arcsec based on ra/dec
    ###truncated to dp so compatible with IAU standards
    ra, dec = np.array(ra), np.array(dec)
    
    cat = SkyCoord(ra=ra, dec=dec, unit='deg')
    
    astring = cat.ra.to_string(sep='', unit='hour', precision=dp, pad=True)
    dstring = cat.dec.to_string(sep='', precision=dp, pad=True, alwayssign=True)
    
    ###truncation index
    tinda = aprec+7
    tindd = tinda
    if aprec==1:
        tindd = tindd-1
    
    if prefix == '':
        sname = ['J' + astring[i][:tinda] + dstring[i][:tindd] for i in range(len(astring))]
    else:
        sname = [prefix + ' J' + astring[i][:tinda] + dstring[i][:tindd] for i in
                 range(len(astring))]


    return(sname)


def find_pairs(data, columndict = {'name': 'Component_name', 'ra': 'RA',
               'dec': 'DEC', 'peak': 'Peak_flux', 'total': 'Total_flux',
               'maj': 'DC_Maj', 'min': 'DC_Min', 'pa': 'DC_PA'},
               only_unique_pairs=True):
    'find pair candidates, include basic neighbour info, sepatation, PA, and sep to 3rdNN'
    
    ###limit to only necessary columns - make generic
    #    usecols = ['Component_name', 'RA', 'DEC', 'Peak_flux', 'Total_flux',
    #               'DC_Maj', 'DC_Min', 'DC_PA']
    usecols = list(columndict.values())
    
    ###setup paircols for hstack whilst at it
    p1cols = [i + '_1' for i in usecols]
    p2cols = [i + '_2' for i in usecols]
    
    pairdata = data[usecols]
    
    acol = columndict['ra']
    dcol = columndict['dec']
    spcol = columndict['peak']
    stcol = columndict['total']
    majcol = columndict['maj']
    mincol = columndict['min']
    pacol = columndict['pa']
    
    ###add units to positions and pos angles - must be decimal degs (makes work if not fits/votable)
    pairdata[acol].unit = 'deg'
    pairdata[dcol].unit = 'deg'
    pairdata[pacol].unit = 'deg'
    
    ###sky position and find neighbours
    poscat = SkyCoord(ra=pairdata[acol], dec=pairdata[dcol])
    
    nn1 = match_coordinates_sky(poscat, poscat, nthneighbor=2)
    nn2 = match_coordinates_sky(poscat, poscat, nthneighbor=3)
    
    
    ###create table of pairs
    p1, p2 = pairdata.copy(), pairdata.copy()
    p1.rename_columns(names=usecols, new_names=p1cols)
    p2.rename_columns(names=usecols, new_names=p2cols)
    
    pairs = hstack([p1, p2[nn1[0]]])
    
    c1 = SkyCoord(ra=pairs[acol+'_1'], dec=pairs[dcol+'_1'])
    c2 = SkyCoord(ra=pairs[acol+'_2'], dec=pairs[dcol+'_2'])
    
    
    ###add in pair separation, position angle, and dist to 2nd NN, flux ratios, centroid
    ###calculate position angles
    pa_pair = c1.position_angle(c2).to(u.deg)
    sep_pair = c1.separation(c2)
    centroid = c1.directional_offset_by(pa_pair, sep_pair/2) ###need to add flux weighted
    flux_weighting = np.array(pairs[stcol+'_1']/(pairs[stcol+'_1'] + pairs[stcol+'_2']))
    ###use 1-flux weighting to bring FW centroid closer to brightest component
    fwcentroid = c1.directional_offset_by(pa_pair, sep_pair*(1-flux_weighting))
    
    
    dpa1 = np.array(pairs[pacol+'_1'] - pa_pair)
    dpa2 = np.array(pairs[pacol+'_2'] - pa_pair)
    dpa = [dpa1, dpa2]
    adpa = []
    for pa in dpa:
        pa[pa<-180] = pa[pa<-180]+360
        pa[pa<-90] = pa[pa<-90]+180
        pa[pa>90] = pa[pa>90]-180
        adpa.append(np.abs(pa))
    #        pa = pa*u.deg

    
    pairs['Sep_pair'] = nn1[1].to(u.arcsec)
    pairs['PA_pair'] = pa_pair
    pairs['dPA_1'] = dpa[0]*u.deg
    pairs['dPA_2'] = dpa[1]*u.deg
    pairs['abs_dPA_1'] = adpa[0]*u.deg
    pairs['abs_dPA_2'] = adpa[1]*u.deg
    pairs['Tflux_ratio'] = np.array(pairs[stcol+'_1']/pairs[stcol+'_2'])
    pairs['Pflux_ratio'] = np.array(pairs[spcol+'_1']/pairs[spcol+'_2'])
    pairs['medianRA'] = centroid.ra
    pairs['medianDEC'] = centroid.dec
    pairs['cenRA'] = fwcentroid.ra
    pairs['cenDEC'] = fwcentroid.dec
    pairs['d_2NN'] = nn2[1].to(u.arcsec)
    pairs['pair_name'] = source_name(ra=pairs['cenRA'], dec=pairs['cenDEC'])

    if only_unique_pairs==True:
        pairs = unique(pairs, 'pair_name')
    
    return(pairs)


def preprocess_cat(catfile, flux_col, flux_min, size_col, size_min, fileformat='fits'):
    'preprocess radio catalogue to filter only components matching desired criteria'
    alldata = Table.read(catfile, format=fileformat)
    data = alldata[(alldata[flux_col]>=flux_min) & (alldata[size_col]>=size_min)]
    
    return(data)


def hunt_dragns(catfile, config_file, write_file=True):
    'find pairs and flag candidate dragns'
    
    ##extract config params
    catparams, coldict = config_parse(config_file=config_file)
    
    ##load and preprocesses cat data
    catdata = preprocess_cat(catfile=catfile, flux_col=coldict['peak'],
                             flux_min=catparams['flux_min'], size_col=coldict['maj'],
                             size_min=catparams['lobesize_min'],
                             fileformat=catparams['file_format'])
    
    ###find pairs and flag (need to add flagging)
    pairs = find_pairs(data=catdata, columndict=coldict)
    
    ###write to file or return data in code
    if write_file==True:
        ##create new file name based on old
        splitname = catfile.split('.')
        prename, extension = splitname[0], '.' + splitname[len(splitname)-1]
        newname = (prename + '_flux' + str(int(catparams['flux_min'])) + '_size'
                   + str(int(catparams['lobesize_min'])) + '_pairs' + extension)
        pairs.write(newname, format=catparams['file_format'])
        return
    else:
        return(pairs)


def config_parse(config_file):
    'extract required info from config file'
    cdata = pd.read_table(config_file, delim_whitespace=True).replace("'", "", regex=True)
    
    ##parameters for catalog preprocessing
    minflux = np.float(cdata['value'][np.where(cdata['parameter']=='min_brightness')[0][0]])
    minsizelobe = np.float(cdata['value'][np.where(cdata['parameter']=='min_size_lobe')[0][0]])
    maxsizecore = np.float(cdata['value'][np.where(cdata['parameter']=='max_size_core')[0][0]])
    dformat = cdata['value'][np.where(cdata['parameter']=='data_format')[0][0]]
    
    catparams = {'flux_min': minflux, 'lobesize_min': minsizelobe,
                 'coresize_max': maxsizecore, 'file_format':dformat}
    
    ##catalogue column names
    namecol = cdata['value'][np.where(cdata['parameter']=='name_col')[0][0]]
    racol = cdata['value'][np.where(cdata['parameter']=='ra_col')[0][0]]
    deccol = cdata['value'][np.where(cdata['parameter']=='dec_col')[0][0]]
    peakcol = cdata['value'][np.where(cdata['parameter']=='Speak_col')[0][0]]
    totcol = cdata['value'][np.where(cdata['parameter']=='Stot_col')[0][0]]
    majcol = cdata['value'][np.where(cdata['parameter']=='maj_col')[0][0]]
    mincol = cdata['value'][np.where(cdata['parameter']=='min_col')[0][0]]
    pacol = cdata['value'][np.where(cdata['parameter']=='pa_col')[0][0]]
    
    coldict = {'name': namecol, 'ra': racol, 'dec': deccol, 'peak': peakcol,
               'total': totcol, 'maj': majcol, 'min': mincol, 'pa': pacol}
    
    return(catparams, coldict)


def parse_args():
    "parse input args, i.e. target and config file names"
    parser = argparse.ArgumentParser(description="find radio pairs in catalogue data")
    parser.add_argument("radio_cat", help="radio component catalogue to find pairs in")
    parser.add_argument("--config", action='store', type=str, default='config.txt',
                        help="config file")
    parser.add_argument("--writepairs", action='store', type=str, default='True',
                        help="if True write pairs to file")
                                     
    args = parser.parse_args()
    
    ##make args.writepairs a bool
    if args.writepairs == 'True':
        args.writepairs = True
    
    return args


##############################################################################
###main code

if __name__ == '__main__':
    args = parse_args()
    if args.writepairs==True:
        hunt_dragns(catfile=args.radio_cat, config_file=args.config, write_file=True)
    else:
        pairs = hunt_dragns(catfile=args.radio_cat, config_file=args.config, write_file=False)
