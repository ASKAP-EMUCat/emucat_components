import asyncio
import asyncpg


async def match_with_allwise(conn, ser_name: str, max_separation_rads: float):
    sql = "INSERT INTO emucat.sources (component_id, wise_id, separation) " \
          "select c.id, a.designation, a.distance " \
          "from emucat.components c, emucat.mosaics m, emucat.source_extraction_regions s,  " \
          "lateral " \
          "(select aw.designation, aw.ra_dec, " \
          "aw.ra_dec <-> spoint(c.ra_deg_cont*pi()/180.0, c.dec_deg_cont*pi()/180.0) distance " \
          "from emucat.allwise aw " \
          "where " \
          "aw.ra_dec @ scircle(spoint(c.ra_deg_cont*pi()/180.0, c.dec_deg_cont*pi()/180.0), $2)) a " \
          "where c.mosaic_id=m.id and m.ser_id=s.id and s.name=$1 " \
          "ON CONFLICT (component_id, wise_id) DO NOTHING RETURNING id"

    return await conn.fetch(sql, ser_name, max_separation_rads)


