#!python3
import os
import sys
import csv
import argparse
import logging
import asyncpg
import asyncio
import configparser
import pyvo as vo

from pathlib import Path
from dl import queryClient as qc
from dl import authClient as ac
from dl import storeClient as sc


logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s')


def read_credentials(credentials: str):
    config = configparser.ConfigParser()
    config.read(credentials)
    user = config['emucat_database']['user']
    password = config['emucat_database']['password']
    database = config['emucat_database']['database']
    host = config['emucat_database']['host']
    port = config['emucat_database'].getint('port', 5432)
    return user, password, database, host, port

def read_noao_credentials(credentials: str):
    config = configparser.ConfigParser()
    config.read(credentials)
    user = config['noao']['user']
    password = config['noao']['password']
    return user, password

def connect_to_noao():
    return vo.dal.TAPService("https://datalab.noao.edu/tap")


async def _import_des_dr1_from_lhr(ser: str, output: str, credentials: str):
    loop = asyncio.get_event_loop()

    in_table = f"{ser.lower()}_allwise"
    out_table = f"{ser.lower()}_xmatched.csv"

    Path(output).mkdir(parents=True, exist_ok=True)

    try:
        os.remove(f"{output}/{out_table}")
    except OSError:
        pass

    user, password, database, host, port = read_credentials(credentials)

    # Get the lhr sources that dont already exist in des_dr1
    sql =   'SELECT aw.designation, aw.source_id, aw.ra, aw.dec ' \
            'FROM emucat.components c, emucat.mosaics m, emucat.regions s, ' \
            'emucat.sources_lhr_allwise lhr ' \
            'LEFT JOIN emucat.allwise as aw on lhr.wise_id = aw.designation ' \
            'LEFT JOIN emucat.des_dr1_allwise as des on aw.designation = des.wise_id ' \
            'WHERE c.mosaic_id=m.id ' \
            'AND m.ser_id=s.id ' \
            'AND lhr.component_id=c.id ' \
            'AND s.name=$1 ' \
            'AND des.wise_id is NULL '

    insert_conn = None
    conn = None
    try:
        noao = await loop.run_in_executor(None, connect_to_noao)
        insert_conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
        conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
        async with conn.transaction():
            records = await conn.fetch(sql, ser)

        if not records:
            return

        logging.info(f'Getting des_dr1 records from {len(records)} lhr matches.')

        sources = set()
        for row in records:
            sources.add((row['designation'], f"{row['source_id'].strip()}", row['ra'], row['dec']))

        with open(f"{output}/{in_table}.csv", 'w') as csvfile:
            filewriter = csv.writer(csvfile)
            filewriter.writerow(["designation", "source_id", "ra", "dec"])
            for s in sources:
                filewriter.writerow([s[0],s[1],s[2],s[3]])
        
        with open(f"{output}/{in_table}_schema.txt", 'w') as csvfile:
            filewriter = csv.writer(csvfile)
            filewriter.writerow(["designation", "text"])
            filewriter.writerow(["source_id", "text"])
            filewriter.writerow(["ra", "double precision"])
            filewriter.writerow(["dec", "double precision"])

        user, password = read_noao_credentials(credentials)
        token = ac.login(user, password)

        # create remote table at noao and import data
        resp = qc.mydb_create(in_table, f"{output}/{in_table}_schema.txt", drop=True)
        if resp != 'OK':
            raise Exception('mydb_create failed')

        resp = qc.mydb_import(in_table, f"{output}/{in_table}.csv", drop=True)
        if resp != 'OK':
            raise Exception('mydb_import failed')

        #qc.mydb_flush()

        query = "SELECT aw.designation, xm.id1, xm.distance, " \
                "des.ra, des.dec, des.mag_auto_g, des.mag_auto_r, " \
                "des.mag_auto_i, des.mag_auto_z, des.mag_auto_y " \
                f"FROM des_dr1.mag as des, des_dr1.x1p5__main__allwise__source as xm, mydb://{in_table} as aw " \
                "WHERE aw.source_id=xm.id2 AND xm.id1=des.coadd_object_id"

        logging.info(f'Cross matching with des_dr1.')
        print(query)
        a=qc.query(sql=query, fmt='csv', out=f"vos://{out_table}", drop=True,
        async_=True, wait=True, timeout=3000, poll=1)
        print(a)
        sc.get(f"vos://{out_table}", f"{output}/{out_table}")

        logging.info(f'Cross matching with des_dr1 complete.')

        des_rows = []
        xmatch_rows = []
        with open(f"{output}/{out_table}") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            # skip header
            next(csv_reader)
            for row in csv_reader:
                des_rows.append([int(row[1]), float(row[3]), float(row[4]), 
                                    float(row[5]), float(row[6]), float(row[7]),
                                    float(row[8]), float(row[9])])
                xmatch_rows.append([row[0], int(row[1]), float(row[2])])

        logging.info(f'Inserting {len(xmatch_rows)} into des_dr1.')

        async with insert_conn.transaction():
            await insert_conn.executemany('INSERT INTO emucat.des_dr1 ("coadd_object_id", '
                                            '"ra", "dec", "mag_auto_g", "mag_auto_r", '
                                            '"mag_auto_i", "mag_auto_z", "mag_auto_y") '
                                            'VALUES($1, $2, $3, $4, $5, $6, $7, $8) '
                                            'ON CONFLICT ("coadd_object_id") '
                                            'DO NOTHING',
                                            des_rows)

            await insert_conn.executemany('INSERT INTO emucat.des_dr1_allwise '
                                            '("wise_id", "coadd_object_id", "separation") '
                                            'VALUES($1, $2, $3) '
                                            'ON CONFLICT ("wise_id", "coadd_object_id") '
                                            'DO NOTHING',
                                            xmatch_rows)
    finally:
        if insert_conn:
            await insert_conn.close()
        if conn:
            await conn.close()

        try:
            qc.mydb_drop(in_table)
        except:
            pass


def import_des_dr1_from_lhr(args):
    asyncio.run(_import_des_dr1_from_lhr(args.ser, args.output, args.credentials))


def main():
    parser = argparse.ArgumentParser(prog='EMUCat NOAO Functions', description='EMUCat NOAO catalog functions.')
    subparsers = parser.add_subparsers(help='sub-command help')
    parser.set_defaults(func=lambda args: parser.print_help())
    
    import_des_dr1_from_lhr_parser = subparsers.add_parser('import_des_dr1_from_lhr',
                                                           help='Import des_dr1 catalog based on lhr matches.')
    import_des_dr1_from_lhr_parser.add_argument('-s', '--ser', help='Source extraction region.', type=str, required=True)
    import_des_dr1_from_lhr_parser.add_argument('-c', '--credentials', help='Credentials file.', required=True)
    import_des_dr1_from_lhr_parser.add_argument('-o', '--output', help='Output directory', type=str, required=True, default='./')
    import_des_dr1_from_lhr_parser.set_defaults(func=import_des_dr1_from_lhr)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    try:
        main()
        exit(0)
    except Exception as e:
        logging.exception(e)
        exit(1)