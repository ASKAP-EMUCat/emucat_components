#!python3
import xml.etree.ElementTree as ET
import asyncio
import asyncpg
import os
import sys
import random
import string
import logging
import argparse
import configparser
import pyvo as vo

import numpy as np
from astropy.table import Table,vstack,hstack,join,setdiff,unique
from astropy.io.votable import parse_single_table, from_table, writeto
from astropy.coordinates import SkyCoord,ICRS
from astropy import units as u



logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s')


def convert(value, datatype):
    if datatype == 'float':
        return float(value)
    elif datatype == 'boolean':
        return bool(int(value))
    elif datatype == 'int':
        return int(value)
    elif datatype == 'long':
        return int(value)
    elif datatype == 'double':
        return float(value)
    elif datatype == 'string':
        return str(value)
    return value


def read_credentials(credentials):
    config = configparser.ConfigParser()
    config.read(credentials)
    user = config['emucat_database']['user']
    password = config['emucat_database']['password']
    database = config['emucat_database']['database']
    host = config['emucat_database']['host']
    port = config['emucat_database'].getint('port', 5432)
    return (user,password,database,host,port)



def import_properties(args):

    datatype = dict(asyncio.run(db_property_datatypes_select(args.credentials)))

    db_rows =  [['name', 'EMU IAU Name', '', datatype['string']],
                 'ra', 'J2000 Right Ascension', 'deg', datatype['double'],
                 'dec', 'J2000 Declination', 'deg', datatype['double']]

    properties = asyncio.run(db_properties_schema_upsert(db_rows,credentials))


    islands = asyncio.run(db_island_select(args.ser,args.credentials))


    # EMU IAU Name

    name_id = [p[0] for p in properties if p[1] == 'name'][0]

    for row in islands:
        if row['n_components'] > 1:
            name = position_to_EMU_name(row['i_ra_deg_cont'],row['i_dec_deg_cont'],source='I')

        elif row['n_components'] == 1:
            name = position_to_EMU_name(row['c_ra_deg_cont'],row['c_dec_deg_cont'],source='C')

        db_rows = [ name_id, row['sid'], name]
        asyncio.run(db_island_properties_upsert(db_rows,credentials))
 


async def db_property_datatypes_select(credentials):
    user,password,database,host,port = read_credentials(credentials)

    sql = "SELECT datatype,id from emucat.properties_schema_datatypes"

    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        async with conn.transaction():
            result = await conn.fetch(sql)
    finally:
        await conn.close()

    return result



async def db_properties_schema_upsert(rows,credentials):
    user,password,database,host,port = read_credentials(credentials)

    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        async with conn.transaction():
            await conn.executemany('INSERT INTO emucat.properties_schema ("name", "description", "unit", "datatype_id") '
                                'VALUES ($1, $2, $3, $4) '
                                'ON CONFLICT ("name") DO NOTHING',
                                rows)

            result = await conn.fetch("SELECT * from emucat.properties_schema")

    finally:
        await conn.close()
    return results



async def db_island_select(ser,credentials):
    user,password,database,host,port = read_credentials(credentials)

    sql = "SELECT s.id AS sid, \
                  i.ra_deg_cont AS i_ra_deg_cont, \
                  c.ra_deg_cont AS c_ra_deg_cont, \
                  i.dec_deg_cont AS i_dec_deg_cont, \
                  c.dec_deg_cont AS c_dec_deg_cont, \
                  i.n_components \
            FROM emucat.islands i, emucat.components c, emucat.sources_selavy_islands s, \
                emucat.source_extraction_regions ser, emucat.mosaics m  \
            WHERE i.mosaic_id = m.id AND c.mosaic_id = m.id AND s.mosaic_id = m.id \
            AND m.ser_id = ser.id AND c.component_id like '%a' AND s.island_id = i.island_id \
            AND c.island_id = i.island_id AND ser.name = '{0}' LIMIT 10".format(ser)

    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        async with conn.transaction():
            result = await conn.fetch(sql)
            logging.info(len(result))
    finally:
        await conn.close()

    return result



async def db_island_properties_upsert(rows,credentials):
    user,password,database,host,port = read_credentials(credentials)

    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        async with conn.transaction():
            await conn.executemany('INSERT INTO emucat.sources_selavy_islands_properties ("property_id", "source_id", "val") '
                           'VALUES ($1, $2, $3) '
                           'ON CONFLICT ("property_id", "source_id") DO NOTHING',
                           rows)
    finally:
        await conn.close()



# ** EMU IAU Name definition

# the HMS DMS portion of the EMU IAU name
def skycoord_to_string_position(c):
    return  '{0}{1}'.format(c.ra.to_string(unit=u.hourangle, sep='', precision=1, pad=True), 
                            c.dec.to_string(sep='', precision=0, alwayssign=True, pad=True))

# the full EMU IAU name
# version: 'E' early science, 'P' pilot or '1-9' main survey data relase
# type of source: 'C' component, 'I' island or 'S' source
# assume all EMUcat objects are sources in this catalog 
def position_to_EMU_name(ra,dec, version='P', source='S'):
    c = ICRS(ra*u.degree, dec*u.degree)
    radec_string = skycoord_to_string_position(c)
    return 'EMU {0}{1} J{2}'.format(version, source, radec_string)

# **********




def main():
    parser = argparse.ArgumentParser(prog='EMUCat', description='EMUCat catalog functions.')
    subparsers = parser.add_subparsers(help='sub-command help')
    parser.set_defaults(func=lambda args: parser.print_help())

    import_properties_parser = subparsers.add_parser('import_properties',
                                                       help='Import properties into EMUCat.')
    import_properties_parser.add_argument('-s', '--ser', help='Source extraction region.', type=str, required=True)
    import_properties_parser.add_argument('-c', '--credentials', help='Credentials file.', required=True)
    import_properties_parser.set_defaults(func=import_properties)

    args = parser.parse_args()
    args.func(args)



if __name__ == "__main__":
    try:
        main()
        exit(0)
    except Exception as e:
        logging.exception(e)
        exit(1)
