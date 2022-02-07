#!python3
import xml.etree.ElementTree as ET
import asyncio
import asyncpg
import os
import sys
import random
import string
import json
import logging
import argparse
import configparser
import pyvo as vo
import astropy


logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s')


async def db_moasic_select(conn, ser_name: str):
    mosaic_id = await conn.fetchrow('SELECT m.id '
                                    'FROM emucat.regions as se '
                                    'INNER JOIN emucat.mosaics as m '
                                    'ON se.id = m.ser_id '
                                    'where se.name = $1',
                                    ser_name)
    return mosaic_id


async def db_update_mosaic_header(conn, mosaic_id, json_hdr):
    await conn.fetchrow('UPDATE emucat.mosaics SET fits_header=$2 WHERE id=$1',
                        mosaic_id, json_hdr)


async def db_island_upsert_many(conn, row):
    await conn.executemany('INSERT INTO emucat.islands ('
                           '"mosaic_id", "island_id", "island_name", "n_components", "ra_hms_cont", '
                           '"dec_hms_cont", "ra_deg_cont", "dec_deg_cont", "freq", '
                           '"maj_axis", "min_axis", "pos_ang", "flux_int", "flux_int_err", '
                           '"flux_peak", "mean_background", "background_noise", "max_residual", '
                           '"min_residual", "mean_residual", "rms_residual", "stdev_residual", '
                           '"x_min", "x_max", "y_min", "y_max", "n_pix", '
                           '"solid_angle", "beam_area", "x_ave", "y_ave", "x_cen", "y_cen", '
                           '"x_peak", "y_peak", "flag_i1", "flag_i2", "flag_i3", "flag_i4", '
                           '"comment", "ra_dec") '
                           'VALUES($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15,'
                           '$16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29,'
                           '$30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41)',
                           row)


async def db_moasic_upsert(conn, row):
    mosaic_id = await conn.fetchrow('INSERT INTO emucat.mosaics ("ser_id", "table_version", "image_file", '
                                    '"flag_subsection",'
                                    '"subsection", "flag_statsec", "statsec", "search_type", "flag_negative", '
                                    '"flag_baseline", "flag_robuststats", "flag_fdr", "threshold", "flag_growth", '
                                    '"growth_threshold", "min_pix", "min_channels", "min_voxels", "flag_adjacent", '
                                    '"thresh_velocity", "flag_rejectbeforemerge", "flag_twostagemerging", '
                                    '"pixel_centre", "flag_smooth", "flag_atrous", "reference_frequency", '
                                    '"threshold_actual") '
                                    'VALUES($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15,'
                                    '$16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27) '
                                    'ON CONFLICT ("ser_id", "subsection") '
                                    'DO UPDATE SET subsection = EXCLUDED.subsection RETURNING id',
                                    *row)
    return mosaic_id[0]


async def db_mosaic_island_upsert_many(conn, rows):
    await conn.executemany('INSERT INTO emucat.mosaic_islands ("island_id") '
                           'VALUES ($1) '
                           'ON CONFLICT ("island_id") DO NOTHING',
                           rows)


async def db_sources_island_upsert_many(conn, rows):
    await conn.executemany('INSERT INTO emucat.sources_selavy_islands ("mosaic_id", "island_id") '
                           'VALUES ($1, $2) '
                           'ON CONFLICT ("mosaic_id", "island_id") DO NOTHING',
                           rows)


async def db_components_upsert_many(conn, rows):
    await conn.executemany('INSERT INTO emucat.components ("mosaic_id","island_id","component_id",'
                           '"component_name","ra_hms_cont","dec_hms_cont","ra_deg_cont",'
                           '"dec_deg_cont","ra_err","dec_err","freq","flux_peak","flux_peak_err",'
                           '"flux_int","flux_int_err","maj_axis","min_axis","pos_ang","maj_axis_err",'
                           '"min_axis_err","pos_ang_err","maj_axis_deconv","min_axis_deconv",'
                           '"pos_ang_deconv","maj_axis_deconv_err","min_axis_deconv_err",'
                           '"pos_ang_deconv_err","chi_squared_fit","rms_fit_gauss","spectral_index",'
                           '"spectral_curvature","spectral_index_err","spectral_curvature_err",'
                           '"rms_image","has_siblings","fit_is_estimate","spectral_index_from_tt",'
                           '"flag_c4","comment", "ra_dec")'
                           'VALUES($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15,'
                           '$16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29,'
                           '$30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40) '
                           'ON CONFLICT ("mosaic_id", "island_id", "component_name", '
                           '"ra_deg_cont", "dec_deg_cont") DO NOTHING',
                           rows)


async def db_lhr_upsert_many(conn, rows):
    await conn.executemany('INSERT INTO emucat.sources_lhr_allwise ("component_id", "wise_id", "w1_lr",'
                           '"w1_rel","w1_n_cont","w1_separation")'
                           'VALUES($1, $2, $3, $4, $5, $6) '
                           'ON CONFLICT ("component_id", "wise_id") '
                           'DO UPDATE SET '
                           '"w1_lr"=EXCLUDED."w1_lr",'
                           '"w1_rel"=EXCLUDED."w1_rel",'
                           '"w1_n_cont"=EXCLUDED."w1_n_cont",'
                           '"w1_separation"=EXCLUDED."w1_separation"',
                           rows)


async def db_match_nearest_neighbour_with_allwise(conn, ser_name: str, max_separation_rads: float):
    sql = "INSERT INTO emucat.sources_nearest_allwise (component_id, wise_id, separation) " \
          "select c.id, a.designation, a.distance " \
          "from emucat.components c, emucat.mosaics m, emucat.regions s,  " \
          "lateral " \
          "(select aw.designation, aw.ra_dec, " \
          "aw.ra_dec <-> spoint(c.ra_deg_cont*pi()/180.0, c.dec_deg_cont*pi()/180.0) distance " \
          "from emucat.allwise aw " \
          "where " \
          "aw.ra_dec @ scircle(spoint(c.ra_deg_cont*pi()/180.0, c.dec_deg_cont*pi()/180.0), $2)) a " \
          "where c.mosaic_id=m.id and m.ser_id=s.id and s.name=$1 " \
          "ON CONFLICT (component_id, wise_id) DO NOTHING RETURNING id"

    return await conn.fetch(sql, ser_name, max_separation_rads)


async def db_insert_extended_doubles(conn, rows):
    await conn.executemany('INSERT INTO emucat.sources_extended_doubles ("pair_name", "comp_id_1",'
                           '"comp_id_2", "cen_ra_dec")'
                           'VALUES($1, $2, $3, $4) '
                           'ON CONFLICT ("pair_name", "comp_id_1", "comp_id_2") '
                           'DO NOTHING',
                           rows)


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
    return value


def read_credentials(credentials: str):
    config = configparser.ConfigParser()
    config.read(credentials)
    user = config['emucat_database']['user']
    password = config['emucat_database']['password']
    database = config['emucat_database']['database']
    host = config['emucat_database']['host']
    port = config['emucat_database'].getint('port', 5432)
    return user, password, database, host, port


async def import_selavy_island_catalog(conn, ser_name: str, filename: str):
    mosaic_id = await db_moasic_select(conn, ser_name)
    if not mosaic_id:
        raise Exception(f'Mosaic not found for SER: {ser_name}')

    ns = {'ivoa': 'http://www.ivoa.net/xml/VOTable/v1.3'}

    root = await asyncio.get_event_loop().run_in_executor(None, ET.parse, filename)
    datatypes = []
    for field in root.findall('./ivoa:RESOURCE/ivoa:TABLE/ivoa:FIELD', ns):
        datatypes.append(field.get('datatype'))

    rows = []
    for i, tr in enumerate(root.findall('./ivoa:RESOURCE/ivoa:TABLE/ivoa:DATA/ivoa:TABLEDATA/ivoa:TR', ns)):
        cat = [convert(td.text.strip(), datatypes[j]) for j, td in enumerate(tr)]
        cat.insert(0, mosaic_id[0])
        rows.append(cat)
        ra_dec = f"({cat[6]}d, {cat[7]}d)"
        cat.append(ra_dec)

    await db_island_upsert_many(conn, rows)


async def import_selavy_catalog(conn, ser_name: str, filename: str):
    ser_id = await conn.fetchrow('SELECT id from emucat.regions where name=$1', ser_name)
    if not ser_id:
        raise Exception('Source extraction region not found')

    ser_id = ser_id[0]

    ns = {'ivoa': 'http://www.ivoa.net/xml/VOTable/v1.3'}

    root = await asyncio.get_event_loop().run_in_executor(None, ET.parse, filename)

    mosaic_map = {'ser_id': None,
                  'table_version': None,
                  'imageFile': None,
                  'flagSubsection': False,
                  'subsection': None,
                  'flagStatSec': False,
                  'StatSec': None,
                  'searchType': None,
                  'flagNegative': False,
                  'flagBaseline': False,
                  'flagRobustStats': False,
                  'flagFDR': False,
                  'threshold': None,
                  'flagGrowth': False,
                  'growthThreshold': 0.0,
                  'minPix': 0,
                  'minChannels': 0,
                  'minVoxels': 0,
                  'flagAdjacent': False,
                  'threshVelocity': None,
                  'flagRejectBeforeMerge': False,
                  'flagTwoStageMerging': False,
                  'pixelCentre': None,
                  'flagSmooth': False,
                  'flagATrous': False,
                  'Reference frequency': 0.0,
                  'thresholdActual': 0.0}

    for param in root.findall('./ivoa:RESOURCE/ivoa:TABLE/ivoa:PARAM', ns):
        key = param.get('name')
        value = convert(param.get('value'), param.get('datatype'))
        if key == 'imageFile':
            value = os.path.basename(value)
        mosaic_map[key] = value

    mosaic_map['ser_id'] = ser_id
    params = list(mosaic_map.values())
    mosaic_id = await db_moasic_upsert(conn, params)

    datatypes = []
    for field in root.findall('./ivoa:RESOURCE/ivoa:TABLE/ivoa:FIELD', ns):
        datatypes.append(field.get('datatype'))

    rows = []
    for i, tr in enumerate(root.findall('./ivoa:RESOURCE/ivoa:TABLE/ivoa:DATA/ivoa:TABLEDATA/ivoa:TR', ns)):
        cat = [convert(td.text.strip(), datatypes[j]) for j, td in enumerate(tr)]
        cat.insert(0, mosaic_id)
        ra_dec = f"({cat[6]}d, {cat[7]}d)"
        cat.append(ra_dec)
        rows.append(cat)

    island_rows = []
    source_island_rows = []
    for row in rows:
        # mosaic_id and island_id
        island_rows.append([row[1]])
        source_island_rows.append([row[0], row[1]])

    async with conn.transaction():
        await db_mosaic_island_upsert_many(conn, island_rows)
        await db_sources_island_upsert_many(conn, source_island_rows)
        await db_components_upsert_many(conn, rows)


async def import_lhr_catalog(conn, filename: str):

    ns = '{http://www.ivoa.net/xml/VOTable/v1.4}'

    data_types = [int, str, float, float, float, float]
    rows = []
    row = []
    type_count = 0

    for event, elem in ET.iterparse(filename, events=('start', 'end', 'start-ns', 'end-ns')):
        if event == 'end':
            if elem.tag == f"{ns}TR":
                # remove q_warning element
                #row.pop()
                rows.append(row)
                elem.clear()
            elif elem.tag == f"{ns}TD":
                value = data_types[type_count](elem.text)
                row.append(value)
                type_count += 1
        elif event == 'start':
            if elem.tag == f"{ns}TR":
                row = []
                type_count = 0

    async with conn.transaction():
        await db_lhr_upsert_many(conn, rows)


def import_lhr(args):
    asyncio.run(import_lhr_votable(args.input, args.credentials))


async def import_extended_double_catalog(conn, filename):
    ns = '{http://www.ivoa.net/xml/VOTable/v1.4}'

    data_types = [str, int, int, float, float]
    rows = []
    row = []
    type_count = 0

    for event, elem in ET.iterparse(filename, events=('start', 'end', 'start-ns', 'end-ns')):
        if event == 'end':
            if elem.tag == f"{ns}TR":
                dec = row.pop()
                ra = row.pop()
                point = f"({ra}d, {dec}d)"
                row.append(point)
                rows.append(row)
                elem.clear()
            elif elem.tag == f"{ns}TD":
                value = data_types[type_count](elem.text)
                row.append(value)
                type_count += 1
        elif event == 'start':
            if elem.tag == f"{ns}TR":
                row = []
                type_count = 0

    async with conn.transaction():
        await db_insert_extended_doubles(conn, rows)


async def _import_extended_doubles(filename: str, credentials: str):
    user, password, database, host, port = read_credentials(credentials)
    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        await import_extended_double_catalog(conn, filename)
    finally:
        await conn.close()


def import_extended_doubles(args):
    asyncio.run(_import_extended_doubles(args.input, args.credentials))


async def import_lhr_votable(filename: str, credentials: str):
    user, password, database, host, port = read_credentials(credentials)
    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        await import_lhr_catalog(conn, filename)
    finally:
        await conn.close()


def import_selavy_island(args):
    asyncio.run(import_selavy_island_votable(args.ser, args.input, args.credentials))


async def import_selavy_island_votable(ser_name: str, filename: str, credentials: str):
    user, password, database, host, port = read_credentials(credentials)
    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        await import_selavy_island_catalog(conn, ser_name, filename)
    finally:
        await conn.close()


def import_selavy(args):
    asyncio.run(import_selavy_votable(args.ser, args.input, args.credentials))


async def import_selavy_votable(ser_name: str, filename: str, credentials: str):
    user, password, database, host, port = read_credentials(credentials)
    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        await import_selavy_catalog(conn, ser_name, filename)
    finally:
        await conn.close()


def random_word(length):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(length))


async def import_allwise_catalog_from_csv(input_path: str):
    conn = await asyncpg.connect(user='admin', password='admin', database='emucat', host='localhost')
    async with conn.transaction():
        with open(input_path, 'rb') as f:
            table_name = random_word(6)
            temp_sql = f'CREATE TEMP TABLE {table_name} ON COMMIT DROP AS SELECT * FROM emucat.allwise WITH NO DATA;'
            await conn.fetchrow(temp_sql)
            copy_result = await conn.copy_to_table(table_name=table_name, source=f, format='csv')
            print(copy_result)
            copy_sql = f'INSERT INTO emucat.allwise SELECT * FROM {table_name} ORDER BY (designation) ' \
                       f'ON CONFLICT (designation) DO NOTHING'
            await conn.fetchrow(copy_sql)


async def _delete_components(ser: str, credentials: str):
    user, password, database, host, port = read_credentials(credentials)

    sql = 'DELETE FROM emucat.mosaics m WHERE m.ser_id IN ' \
          '(SELECT s.id FROM emucat.regions s WHERE s.name=$1)'

    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        async with conn.transaction():
            result = await conn.execute(sql, ser)
            logging.info(result)
    finally:
        await conn.close()


async def _match_nearest_neighbour_with_allwise(ser: str, credentials: str):
    user, password, database, host, port = read_credentials(credentials)

    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        async with conn.transaction():
            result = await db_match_nearest_neighbour_with_allwise(conn=conn,
                                                                   ser_name=ser,
                                                                   max_separation_rads=1.93925472249408e-5)
            logging.info(len(result))
    finally:
        await conn.close()


async def _check_preconditions(sbid: int, credentials: str):
    user, password, database, host, port = read_credentials(credentials)

    sql_update = 'UPDATE emucat.scheduling_blocks SET deposited = true WHERE sb_num = $1'
    sql_select = 'SELECT ser.name, ' \
                 'count(sb.sb_num) as sbid_count, ' \
                 'sum(case when sb.deposited = true then 1 else 0 end) as deposited_count ' \
                 'FROM emucat.regions as ser, ' \
                 'emucat.mosaic_prerequisites as mp, ' \
                 'emucat.scheduling_blocks as sb ' \
                 'WHERE ser.id = mp.ser_id and mp.sb_id = sb.id ' \
                 'and ser.submitted = false ' \
                 'GROUP BY ser.name'

    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        async with conn.transaction():
            await conn.execute(sql_update, sbid)

        async with conn.transaction():
            results = await conn.fetch(sql_select)

        ser = []
        for result in results:
            if result['sbid_count'] == result['deposited_count']:
                ser.append(result['name'])

        print(' '.join(ser), end='')

    finally:
        await conn.close()


async def _extract_insert_mosaic_header(ser: str, filename: str, credentials: str):
    json_hdr = {}
    with astropy.io.fits.open(filename) as hdu:
        for head in hdu[0].header:
            if isinstance(hdu[0].header[head], astropy.io.fits.header._HeaderCommentaryCards):
                json_hdr[head] = repr(hdu[0].header[head])
            else:
                json_hdr[head] = hdu[0].header[head]

    user, password, database, host, port = read_credentials(credentials)
    conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
    try:
        async with conn.transaction():
            mosaic_id = await db_moasic_select(conn, ser)
            if not mosaic_id:
                raise Exception(f'Mosaic not found for SER: {ser}')
            await db_update_mosaic_header(conn, mosaic_id[0], json.dumps(json_hdr))
    finally:
        await conn.close()


def delete_components(args):
    asyncio.run(_delete_components(args.ser, args.credentials))


def match_nearest_neighbour_with_allwise(args):
    asyncio.run(_match_nearest_neighbour_with_allwise(args.ser, args.credentials))


def check_preconditions(args):
    asyncio.run(_check_preconditions(int(args.sbid), args.credentials))


def extract_insert_mosaic_header(args):
    asyncio.run(_extract_insert_mosaic_header(args.ser, args.input, args.credentials))


def main():
    parser = argparse.ArgumentParser(prog='EMUCat', description='EMUCat catalog functions.')
    subparsers = parser.add_subparsers(help='sub-command help')
    parser.set_defaults(func=lambda args: parser.print_help())

    preconditions = subparsers.add_parser('emucat_preconditions',
                                          help='Check CASDA event against EMU preconditions and execute pipeline.')
    preconditions.add_argument('-s', '--sbid',
                               help='ASKAP scheduling block ID.',
                               type=str, required=True)
    preconditions.add_argument('-c', '--credentials',
                               help='Credentials file.', required=True)
    preconditions.set_defaults(func=check_preconditions)

    input_selavy_parser = subparsers.add_parser('import_selavy', help='Import selavy component catalog into EMUCat.')
    input_selavy_parser.add_argument('-s', '--ser', help='Source extraction region.', type=str, required=True)
    input_selavy_parser.add_argument('-i', '--input', help='Selavy votable.', type=str, required=True)
    input_selavy_parser.add_argument('-c', '--credentials', help='Credentials file.', required=True)
    input_selavy_parser.set_defaults(func=import_selavy)

    input_selavy_island_parser = subparsers.add_parser('import_selavy_island',
                                                       help='Import selavy islands into EMUCat.')
    input_selavy_island_parser.add_argument('-s', '--ser', help='Source extraction region.', type=str, required=True)
    input_selavy_island_parser.add_argument('-i', '--input', help='Selavy island votable.', type=str, required=True)
    input_selavy_island_parser.add_argument('-c', '--credentials', help='Credentials file.', required=True)
    input_selavy_island_parser.set_defaults(func=import_selavy_island)

    input_lhr_parser = subparsers.add_parser('import_lhr', help='Import lhr results into EMUCat.')
    input_lhr_parser.add_argument('-i', '--input', help='lhr votable.', type=str, required=True)
    input_lhr_parser.add_argument('-c', '--credentials', help='Credentials file.', required=True)
    input_lhr_parser.set_defaults(func=import_lhr)

    import_extended_doubles_parser = subparsers.add_parser('import_extended_doubles',
                                                       help='Import extended doubles into EMUCat.')
    import_extended_doubles_parser.add_argument('-i', '--input', help='extended doubles votable.', type=str,
                                                required=True)
    import_extended_doubles_parser.add_argument('-c', '--credentials', help='Credentials file.', required=True)
    import_extended_doubles_parser.set_defaults(func=import_extended_doubles)

    delete_components_parser = subparsers.add_parser('delete_components',
                                                     help='Delete all components of a SER.')
    delete_components_parser.add_argument('-s', '--ser', help='Source extraction region.', type=str, required=True)
    delete_components_parser.add_argument('-c', '--credentials', help='Credentials file.', required=True)
    delete_components_parser.set_defaults(func=delete_components)

    match_nearest_neighbour_with_allwise_parser = subparsers.add_parser('match_nearest_neighbour_with_allwise',
                                                                        help='Match nearest neighbour with allwise.')
    match_nearest_neighbour_with_allwise_parser.add_argument('-s', '--ser', help='Source extraction region.',
                                                             type=str, required=True)
    match_nearest_neighbour_with_allwise_parser.add_argument('-c', '--credentials',
                                                             help='Credentials file.', required=True)
    match_nearest_neighbour_with_allwise_parser.set_defaults(func=match_nearest_neighbour_with_allwise)

    extract_insert_mosaic_header_parser = subparsers.add_parser('extract_insert_mosaic_header',
                                                                help='Extract mosaic header and insert into Emucat.')
    extract_insert_mosaic_header_parser.add_argument('-s', '--ser', help='Source extraction region.',
                                                     type=str, required=True)
    extract_insert_mosaic_header_parser.add_argument('-i', '--input', help='Mosaic fits file.', type=str,
                                                     required=True)
    extract_insert_mosaic_header_parser.add_argument('-c', '--credentials',
                                                     help='Credentials file.', required=True)
    extract_insert_mosaic_header_parser.set_defaults(func=extract_insert_mosaic_header)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    try:
        main()
        exit(0)
    except Exception as e:
        logging.exception(e)
        exit(1)
