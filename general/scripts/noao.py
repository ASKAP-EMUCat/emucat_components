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
import copy

from astropy.io import ascii
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


async def _import_vhs_from_lhr(ser: str, output: str, credentials: str):

    ser_mod = copy.deepcopy(ser).replace('-', '_').replace('+', '_').lower()

    in_table = f"vhs_allwise_{ser_mod}"
    out_table = f"vhs_xmatched_{ser_mod}"

    Path(output).mkdir(parents=True, exist_ok=True)

    try:
        os.remove(f"{output}/{out_table}")
    except OSError:
        pass

    user, password, database, host, port = read_credentials(credentials)

    # Get the lhr sources that dont already exist in vhs
    sql =  '''
           SELECT aw.designation, trim(aw.source_id) as source_id 
           FROM emucat.components c, emucat.mosaics m, emucat.regions s, 
           emucat.sources_lhr_allwise lhr 
           LEFT JOIN emucat.allwise as aw on lhr.wise_id=aw.designation 
           LEFT JOIN emucat.vhs_dr5_allwise as vhs on aw.source_id=vhs.aw_source_id 
           WHERE c.mosaic_id=m.id 
           AND m.ser_id=s.id 
           AND lhr.component_id=c.id 
           AND s.name=$1 
           AND vhs.aw_source_id is NULL 
           GROUP BY aw.designation 
           ORDER BY aw.designation ASC
           '''

    insert_conn = None
    conn = None
    try:
        logging.info(f'Getting lhr sources that dont already exist in vhs.')

        conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
        records = await conn.fetch(sql, ser)
        if not records:
            return

        logging.info(f'Writing vhs records from {len(records)} lhr matches.')

        with open(f"{output}/{in_table}.csv", 'w') as csvfile:
            filewriter = csv.writer(csvfile)
            filewriter.writerow(["source_id"])
            for s in records:
                filewriter.writerow([f"{s[1]}"])

        with open(f"{output}/{in_table}_schema.txt", 'w') as csvfile:
            filewriter = csv.writer(csvfile)
            filewriter.writerow(["source_id", "text"])

        noao_user, noao_password = read_noao_credentials(credentials)
        token = ac.login(noao_user, noao_password)

        logging.info(f'Uploading table.')

        # create remote table at noao and import data
        resp = qc.mydb_create(in_table, f"{output}/{in_table}_schema.txt", drop=True)
        if resp != 'OK':
            raise Exception('mydb_create failed')

        resp = qc.mydb_import(in_table, f"{output}/{in_table}.csv", drop=True)
        if resp != 'OK':
            raise Exception('mydb_import failed')

        logging.info(f'Uploading complete.')
        
        query = "SELECT xm.*, vhs.* " \
                f"FROM vhs_dr5.vhs_cat_v3 as vhs, vhs_dr5.x1p5__vhs_cat_v3__allwise__source as xm, mydb://{in_table} as aw " \
                "WHERE xm.id1=vhs.sourceid AND aw.source_id=xm.id2"

        logging.info(f'Cross matching with vhs.')
    
        resp = qc.query(sql=query, fmt='table', out=f"vos://{out_table}", drop=True, timeout=600)
        if resp != 'OK':
            raise Exception(f'query: {resp}')
            
        sc.get(f"vos://{out_table}", f"{output}/{out_table}")

        logging.info(f'Cross matching with vhs complete.')

        vhs_rows = []
        xmatch_rows = []
        data = ascii.read( f"{output}/{out_table}")
        for r in data:
            vhs_rows.append([r['sourcename'], r['sourceid'],r['cueventid'],
            r['framesetid'],r['ra2000'],r['dec2000'],r['l'],
            r['b'],r['lambda'],r['eta'],r['priorsec'],
            r['ymjpnt'],r['ymjpnterr'],r['jmhpnt'],r['jmhpnterr'],
            r['hmkspnt'],r['hmkspnterr'],r['jmkspnt'],r['jmkspnterr'],
            r['ymjext'],r['ymjexterr'],r['jmhext'],r['jmhexterr'],
            r['hmksext'],r['hmksexterr'],r['jmksext'],r['jmksexterr'],
            r['mergedclassstat'],r['mergedclass'],r['pstar'],r['pgalaxy'],
            r['pnoise'],r['psaturated'],r['ebv'],r['ay'],
            r['aj'],r['ah'],r['aks'],r['ymjd'],
            r['ypetromag'],r['ypetromagerr'],r['yapermag3'],r['yapermag3err'],
            r['yapermag4'],r['yapermag4err'],r['yapermag6'],r['yapermag6err'],
            r['yapermagnoapercorr3'],r['yapermagnoapercorr4'],
            r['yapermagnoapercorr6'],r['yhlcorsmjradas'],
            r['ygausig'],r['yell'],r['ypa'],r['yerrbits'],
            r['yaverageconf'],r['yclass'],r['yclassstat'],r['ypperrbits'],
            r['yseqnum'],r['yxi'],r['yeta'],r['jmjd'],
            r['jpetromag'],r['jpetromagerr'],r['japermag3'],r['japermag3err'],
            r['japermag4'],r['japermag4err'],r['japermag6'],r['japermag6err'],
            r['japermagnoapercorr3'],r['japermagnoapercorr4'],
            r['japermagnoapercorr6'],r['jhlcorsmjradas'],
            r['jgausig'],r['jell'],r['jpa'],r['jerrbits'],
            r['javerageconf'],r['jclass'],r['jclassstat'],r['jpperrbits'],
            r['jseqnum'],r['jxi'],r['jeta'],r['hmjd'],
            r['hpetromag'],r['hpetromagerr'],r['hapermag3'],r['hapermag3err'],
            r['hapermag4'],r['hapermag4err'],r['hapermag6'],r['hapermag6err'],
            r['hapermagnoapercorr3'],r['hapermagnoapercorr4'],
            r['hapermagnoapercorr6'],r['hhlcorsmjradas'],
            r['hgausig'],r['hell'],r['hpa'],r['herrbits'],
            r['haverageconf'],r['hclass'],r['hclassstat'],r['hpperrbits'],
            r['hseqnum'],r['hxi'],r['heta'],r['ksmjd'],
            r['kspetromag'],r['kspetromagerr'],r['ksapermag3'],r['ksapermag3err'],
            r['ksapermag4'],r['ksapermag4err'],r['ksapermag6'],r['ksapermag6err'],
            r['ksapermagnoapercorr3'],r['ksapermagnoapercorr4'],
            r['ksapermagnoapercorr6'],r['kshlcorsmjradas'],
            r['ksgausig'],r['ksell'],r['kspa'],r['kserrbits'],
            r['ksaverageconf'],r['ksclass']])
            
            xmatch_rows.append([r['id1'], r['ra1'], r['ra2'], r['id2'], r['ra2'], r['dec2'], r['distance']])

        logging.info(f'Inserting {len(xmatch_rows)} into vhs.')

        INSERT = '''
                 INSERT INTO emucat.vhs_dr5_cat_v3 (
                    sourcename,sourceid,
                    cueventid,framesetid,
                    ra2000,dec2000,l,b,lambda,
                    eta,priorsec,ymjpnt,ymjpnterr,jmhpnt,
                    jmhpnterr,hmkspnt,hmkspnterr,jmkspnt,
                    jmkspnterr,ymjext,ymjexterr,jmhext,
                    jmhexterr,hmksext,hmksexterr,jmksext,
                    jmksexterr,mergedclassstat,mergedclass,pstar,
                    pgalaxy,pnoise,psaturated,ebv,
                    ay,aj,ah,aks,ymjd,ypetromag,
                    ypetromagerr,yapermag3,yapermag3err,yapermag4,
                    yapermag4err,yapermag6,yapermag6err,yapermagnoapercorr3,
                    yapermagnoapercorr4,yapermagnoapercorr6,
                    yhlcorsmjradas,ygausig,yell,ypa,
                    yerrbits,yaverageconf,yclass,yclassstat,
                    ypperrbits,yseqnum,yxi,yeta,jmjd,
                    jpetromag,jpetromagerr,japermag3,japermag3err,
                    japermag4,japermag4err,japermag6,japermag6err,
                    japermagnoapercorr3,japermagnoapercorr4,
                    japermagnoapercorr6,jhlcorsmjradas,
                    jgausig,jell,jpa,jerrbits,
                    javerageconf,jclass,jclassstat,jpperrbits,
                    jseqnum,jxi,jeta,hmjd,hpetromag,
                    hpetromagerr,hapermag3,hapermag3err,hapermag4,
                    hapermag4err,hapermag6,hapermag6err,hapermagnoapercorr3,
                    hapermagnoapercorr4,hapermagnoapercorr6,
                    hhlcorsmjradas,hgausig,hell,hpa,
                    herrbits,haverageconf,hclass,hclassstat,
                    hpperrbits,hseqnum,hxi,heta,
                    ksmjd,kspetromag,kspetromagerr,ksapermag3,
                    ksapermag3err,ksapermag4,ksapermag4err,ksapermag6,
                    ksapermag6err,ksapermagnoapercorr3,ksapermagnoapercorr4,ksapermagnoapercorr6,
                    kshlcorsmjradas,ksgausig,ksell,kspa,kserrbits,ksaverageconf,ksclass) 
                    VALUES ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,
                            $12,$13,$14,$15,$16,$17,$18,$19,
                            $20,$21,$22,$23,$24,$25,$26,$27,
                            $28,$29,$30,$31,$32,$33,$34,$35,
                            $36,$37,$38,$39,$40,$41,$42,$43,
                            $44,$45,$46,$47,$48,$49,$50,$51,
                            $52,$53,$54,$55,$56,$57,$58,$59,
                            $60,$61,$62,$63,$64,$65,$66,$67,
                            $68,$69,$70,$71,$72,$73,$74,$75,
                            $76,$77,$78,$79,$80,$81,$82,$83,
                            $84,$85,$86,$87,$88,$89,$90,$91,
                            $92,$93,$94,$95,$96,$97,$98,$99,
                            $100,$101,$102,$103,$104,$105,$106,$107,
                            $108,$109,$110,$111,$112,$113,$114,$115,
                            $116,$117,$118,$119,$120,$121,$122,$123,
                            $124,$125,$126,$127,$128,$129) 
                    ON CONFLICT (sourceid) 
                    DO NOTHING 
                 '''

        insert_conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
        async with insert_conn.transaction():
            await insert_conn.executemany(INSERT, vhs_rows)
            await insert_conn.executemany('''INSERT INTO emucat.vhs_dr5_allwise 
                                             (vhs_id, vhs_ra, vhs_dec, aw_source_id, aw_ra, aw_dec, separation) 
                                             VALUES($1, $2, $3, $4, $5, $6, $7) 
                                             ON CONFLICT (vhs_id, aw_source_id) 
                                             DO NOTHING''', xmatch_rows)
    finally:
        if insert_conn:
            await insert_conn.close()
        if conn:
            await conn.close()

        try:
            qc.mydb_drop(in_table)
        except:
            pass

        try:
            sc.rm(out_table)
        except:
            pass


async def _import_des_dr1_from_lhr(ser: str, output: str, credentials: str):
    loop = asyncio.get_event_loop()

    ser_mod = copy.deepcopy(ser).replace('-', '_').replace('+', '_').lower()

    in_table = f"{ser_mod}_des_dr1_allwise"
    out_table = f"{ser_mod}_des_dr1_xmatched.csv"

    Path(output).mkdir(parents=True, exist_ok=True)

    try:
        os.remove(f"{output}/{out_table}")
    except OSError:
        pass

    user, password, database, host, port = read_credentials(credentials)

    # Get the lhr sources that dont already exist in des_dr1
    sql =   """
            SELECT aw.designation, aw.source_id, aw.ra, aw.dec 
            FROM emucat.components c, emucat.mosaics m, emucat.regions s, 
            emucat.sources_lhr_allwise lhr 
            LEFT JOIN emucat.allwise as aw on lhr.wise_id = aw.designation 
            LEFT JOIN emucat.des_dr1_allwise as des on aw.designation = des.wise_id 
            WHERE c.mosaic_id=m.id 
            AND m.ser_id=s.id 
            AND lhr.component_id=c.id 
            AND s.name=$1 
            AND des.wise_id is NULL 
            GROUP BY aw.designation 
            ORDER BY aw.designation ASC
            """

    insert_conn = None
    conn = None
    try:
        logging.info(f'Getting lhr sources that dont already exist in des_dr1.')

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
        resp = qc.query(sql=query, fmt='csv', out=f"vos://{out_table}", drop=True, timeout=600) #async_=True, wait=True, timeout=3000, poll=1)
        #resp = qc.query(adql=query, fmt='csv', out=f"vos://{out_table}", drop=True, async_=True, wait=True, timeout=3000, poll=1)
        if resp != 'OK':
            raise Exception(f'query: {resp}')
            
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

        try:
            sc.rm(out_table)
        except:
            pass


async def _import_des_dr2_from_lhr(ser: str, output: str, credentials: str):
    loop = asyncio.get_event_loop()

    ser_mod = copy.deepcopy(ser).replace('-', '_').replace('+', '_').lower()

    in_table = f"{ser_mod}_des_dr2_allwise"
    out_table = f"{ser_mod}_des_dr2_xmatched.csv"

    Path(output).mkdir(parents=True, exist_ok=True)

    try:
        os.remove(f"{output}/{out_table}")
    except OSError:
        pass

    user, password, database, host, port = read_credentials(credentials)

    # Get the lhr sources that dont already exist in des_dr2
    sql =   """
            SELECT distinct(aw.designation), aw.source_id, aw.ra, aw.dec 
            FROM emucat.components c, emucat.mosaics m, emucat.regions s, 
            emucat.sources_lhr_allwise lhr 
            LEFT JOIN emucat.allwise as aw on lhr.wise_id = aw.designation 
            LEFT JOIN emucat.des_dr2_allwise as des on aw.designation = des.wise_id 
            WHERE c.mosaic_id=m.id 
            AND m.ser_id=s.id 
            AND lhr.component_id=c.id 
            AND s.name=$1 
            AND des.wise_id is NULL 
            GROUP BY aw.designation 
            ORDER BY aw.designation ASC
            """

    insert_conn = None
    conn = None
    try:
        logging.info(f'Getting lhr sources that dont already exist in des_dr2.')

        insert_conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
        conn = await asyncpg.connect(user=user, password=password, database=database, host=host, port=port)
        async with conn.transaction():
            records = await conn.fetch(sql, ser)

        if not records:
            return

        logging.info(f'Getting des_dr2 records from {len(records)} lhr matches.')

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

        query = "SELECT aw.designation, xm.id1, xm.distance, " \
                "des.ra, des.dec, des.alphawin_j2000, des.deltawin_j2000, "\
                "des.mag_auto_g_dered, des.mag_auto_i_dered, des.mag_auto_r_dered, "\
                "des.mag_auto_y_dered, des.mag_auto_z_dered, des.ebv_sfd98, "\
                "des.wavg_flux_psf_g, des.wavg_flux_psf_i, des.wavg_flux_psf_r, "\
                "des.wavg_flux_psf_y, des.wavg_flux_psf_z, des.wavg_fluxerr_psf_g, "\
                "des.wavg_fluxerr_psf_i, des.wavg_fluxerr_psf_r, des.wavg_fluxerr_psf_y, "\
                "des.wavg_fluxerr_psf_z, des.flux_auto_g, des.flux_auto_i, des.flux_auto_r, "\
                "des.flux_auto_y, des.flux_auto_z, des.fluxerr_auto_g, "\
                "des.fluxerr_auto_i, des.fluxerr_auto_r, des.fluxerr_auto_y, "\
                "des.fluxerr_auto_z, des.flags_g, des.flags_i, des.flags_r, "\
                "des.flags_y, des.flags_z, des.flux_radius_g, des.flux_radius_i, "\
                "des.flux_radius_r, des.flux_radius_y, des.flux_radius_z, "\
                "des.a_image, des.b_image, des.theta_j2000 "\
                f"FROM des_dr2.main as des, des_dr2.x1p5__main__allwise__source as xm, mydb://{in_table} as aw " \
                "WHERE aw.source_id=xm.id2 AND xm.id1=des.coadd_object_id"

        logging.info(f'Cross matching with des_dr2.')
        resp = qc.query(sql=query, fmt='csv', out=f"vos://{out_table}", drop=True, timeout=600) #async_=True, wait=True, timeout=3000, poll=1)
        if resp != 'OK':
            raise Exception(f'query: {resp}')

        sc.get(f"vos://{out_table}", f"{output}/{out_table}")

        logging.info(f'Cross matching with des_dr2 complete.')

        des_rows = []
        xmatch_rows = []
        with open(f"{output}/{out_table}") as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            # skip header
            next(csv_reader)
            for row in csv_reader:
                des_rows.append([int(row[1]), float(row[3]), float(row[4]), 
                                 float(row[5]), float(row[6]), float(row[7]),
                                 float(row[8]), float(row[9]), float(row[10]),
                                 float(row[11]), float(row[12]), float(row[13]),
                                 float(row[14]), float(row[15]), float(row[16]),
                                 float(row[17]), float(row[18]), float(row[19]),
                                 float(row[20]), float(row[21]), float(row[22]),
                                 float(row[23]), float(row[24]), float(row[25]),
                                 float(row[26]), float(row[27]), float(row[28]),
                                 float(row[29]), float(row[30]), float(row[31]),
                                 float(row[32]), int(row[33]), int(row[34]),
                                 int(row[35]), int(row[36]), int(row[37]),
                                 float(row[38]), float(row[39]), float(row[40]),
                                 float(row[41]), float(row[42]), float(row[43]),
                                 float(row[44]), float(row[45])
                                 ])
                xmatch_rows.append([row[0], int(row[1]), float(row[2])])

        logging.info(f'Inserting {len(xmatch_rows)} into des_dr2.')

        async with insert_conn.transaction():
            await insert_conn.executemany('INSERT INTO emucat.des_dr2 ("coadd_object_id", "ra", "dec", '
                                            'alphawin_j2000, deltawin_j2000, '
                                            'mag_auto_g_dered, mag_auto_i_dered, mag_auto_r_dered, '
                                            'mag_auto_y_dered, mag_auto_z_dered, ebv_sfd98, '
                                            'wavg_flux_psf_g, wavg_flux_psf_i, wavg_flux_psf_r, '
                                            'wavg_flux_psf_y, wavg_flux_psf_z, wavg_fluxerr_psf_g, '
                                            'wavg_fluxerr_psf_i, wavg_fluxerr_psf_r, wavg_fluxerr_psf_y, '
                                            'wavg_fluxerr_psf_z, flux_auto_g, flux_auto_i, flux_auto_r, '
                                            'flux_auto_y, flux_auto_z, fluxerr_auto_g, '
                                            'fluxerr_auto_i, fluxerr_auto_r, fluxerr_auto_y, '
                                            'fluxerr_auto_z, flags_g, flags_i, flags_r, '
                                            'flags_y, flags_z, flux_radius_g, flux_radius_i, '
                                            'flux_radius_r, flux_radius_y, flux_radius_z, '
                                            'a_image, b_image, theta_j2000) '
                                            'VALUES($1, $2, $3, $4, $5, $6, $7, $8, '
                                            '$9, $10, $11, $12, $13, $14, $15, $16, '
                                            '$17, $18, $19, $20, $21, $22, $23, $24, '
                                            '$25, $26, $27, $28, $29, $30, $31, $32, '
                                            '$33, $34, $35, $36, $37, $38, $39, $40, '
                                            '$41, $42, $43, $44) '
                                            'ON CONFLICT ("coadd_object_id") '
                                            'DO NOTHING',
                                            des_rows)

            await insert_conn.executemany('INSERT INTO emucat.des_dr2_allwise '
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

        try:
            sc.rm(out_table)
        except:
            pass


def import_des_dr1_from_lhr(args):
    asyncio.run(_import_des_dr1_from_lhr(args.ser, args.output, args.credentials))


def import_des_dr2_from_lhr(args):
    asyncio.run(_import_des_dr2_from_lhr(args.ser, args.output, args.credentials))


def import_vhs_from_lhr(args):
    asyncio.run(_import_vhs_from_lhr(args.ser, args.output, args.credentials))


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

    import_des_dr2_from_lhr_parser = subparsers.add_parser('import_des_dr2_from_lhr',
                                                           help='Import des_dr2 catalog based on lhr matches.')
    import_des_dr2_from_lhr_parser.add_argument('-s', '--ser', help='Source extraction region.', type=str, required=True)
    import_des_dr2_from_lhr_parser.add_argument('-c', '--credentials', help='Credentials file.', required=True)
    import_des_dr2_from_lhr_parser.add_argument('-o', '--output', help='Output directory', type=str, required=True, default='./')
    import_des_dr2_from_lhr_parser.set_defaults(func=import_des_dr2_from_lhr)

    import_vhs_from_lhr_parser = subparsers.add_parser('import_vhs_from_lhr',
                                                           help='Import vhs catalog based on lhr matches.')
    import_vhs_from_lhr_parser.add_argument('-s', '--ser', help='Source extraction region.', type=str, required=True)
    import_vhs_from_lhr_parser.add_argument('-c', '--credentials', help='Credentials file.', required=True)
    import_vhs_from_lhr_parser.add_argument('-o', '--output', help='Output directory', type=str, required=True, default='./')
    import_vhs_from_lhr_parser.set_defaults(func=import_vhs_from_lhr)


    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    try:
        main()
        exit(0)
    except Exception as e:
        logging.exception(e)
        exit(1)