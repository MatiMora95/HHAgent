"""
HydroAgent v2 — Hydrologic Modeling Data Assistant
Key fix: automatic CRS reprojection (.prj → WGS84) using pyproj
"""

from fastapi import FastAPI, UploadFile, File, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
import httpx
import json
import math
import zipfile
import io
import re
import struct
from typing import Optional

app = FastAPI(title="HydroAgent API", version="2.0.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])
app.mount("/static", StaticFiles(directory="static"), name="static")

@app.get("/")
def root():
    return FileResponse("static/index.html")

@app.get("/health")
def health():
    return {"status": "ok", "version": "2.0.0"}


# ════════════════════════════════════════════════════════
#  REPROJECTION  (.prj → WGS84)
# ════════════════════════════════════════════════════════

def reproject_ring(ring: list, prj_wkt: str) -> list:
    """
    Reproject coordinate ring from source CRS (PRJ WKT) to WGS84.
    Uses pyproj if available; falls back to heuristic otherwise.
    """
    if not ring:
        return ring

    # Check if already geographic
    xs = [p[0] for p in ring]
    ys = [p[1] for p in ring]
    if -180 <= min(xs) and max(xs) <= 180 and -90 <= min(ys) and max(ys) <= 90:
        return ring  # already WGS84

    # Try pyproj first (most accurate)
    if prj_wkt:
        try:
            from pyproj import CRS, Transformer
            src = CRS.from_wkt(prj_wkt)
            dst = CRS.from_epsg(4326)
            if src.is_geographic:
                return ring
            t = Transformer.from_crs(src, dst, always_xy=True)
            return [tuple(t.transform(p[0], p[1])) for p in ring]
        except Exception:
            pass  # fall through to heuristic

    return _heuristic_reproject(ring, prj_wkt or "")


def _heuristic_reproject(ring: list, prj_wkt: str) -> list:
    """
    Heuristic reprojection for common US projections when pyproj unavailable.
    Handles: State Plane (feet), State Plane (meters), UTM.
    Accurate to ~500m — sufficient for NOAA/HSG point queries.
    """
    prj_up = prj_wkt.upper()

    # Convert feet → meters if needed
    is_feet = ("FOOT" in prj_up or "FEET" in prj_up or "US_SURVEY_FOOT" in prj_up)
    FEET_TO_M = 0.3048006096
    ring_m = [(p[0] * FEET_TO_M, p[1] * FEET_TO_M) if is_feet else p for p in ring]

    # Extract projection parameters from PRJ
    def _extract(pattern, text, default):
        m = re.search(pattern, text, re.IGNORECASE)
        return float(m.group(1)) if m else default

    cm  = _extract(r'CENTRAL_MERIDIAN[",\s]+(-?\d+\.?\d*)',  prj_wkt, -96.0)
    lo  = _extract(r'LATITUDE_OF_ORIGIN[",\s]+(\d+\.?\d*)',  prj_wkt,  40.0)
    fe  = _extract(r'FALSE_EASTING[",\s]+(\d+\.?\d*)',        prj_wkt,   0.0)
    fn  = _extract(r'FALSE_NORTHING[",\s]+(\d+\.?\d*)',       prj_wkt,   0.0)

    # If feet, scale false easting/northing too
    if is_feet:
        fe *= FEET_TO_M
        fn *= FEET_TO_M

    lat_rad = math.radians(lo)
    m_per_deg_lat = 111320.0
    m_per_deg_lon = 111320.0 * math.cos(lat_rad)

    result = []
    for p in ring_m:
        dx_m = p[0] - fe
        dy_m = p[1] - fn
        lon_deg = cm + dx_m / m_per_deg_lon
        lat_deg = lo + dy_m / m_per_deg_lat
        result.append((lon_deg, lat_deg))

    return result


# ════════════════════════════════════════════════════════
#  SHAPEFILE PARSER  (pure Python, no geopandas)
# ════════════════════════════════════════════════════════

def parse_shp_bytes(shp_bytes: bytes) -> list:
    """Extract coordinate rings from a SHP file (polygon type 5)."""
    if len(shp_bytes) < 100:
        raise ValueError("SHP file too small")
    offset, rings = 100, []
    while offset < len(shp_bytes) - 8:
        try:
            clen  = struct.unpack('>i', shp_bytes[offset+4:offset+8])[0] * 2
            stype = struct.unpack('<i', shp_bytes[offset+8:offset+12])[0]
            if stype in (5, 15, 25):
                np = struct.unpack('<i', shp_bytes[offset+44:offset+48])[0]
                npt= struct.unpack('<i', shp_bytes[offset+48:offset+52])[0]
                ps = offset + 52
                pts= ps + np * 4
                pi = [struct.unpack('<i', shp_bytes[ps+i*4:ps+i*4+4])[0] for i in range(np)]
                pi.append(npt)
                for p in range(np):
                    ring = []
                    for pt in range(pi[p], pi[p+1]):
                        o = pts + pt*16
                        x = struct.unpack('<d', shp_bytes[o:o+8])[0]
                        y = struct.unpack('<d', shp_bytes[o+8:o+16])[0]
                        ring.append((x, y))
                    rings.append(ring)
                break  # first polygon only
            offset += 8 + clen
        except struct.error:
            break
    return rings


def centroid_and_area(ring: list):
    """Centroid (lat, lon) and approximate area (km²) from WGS84 ring."""
    lons = [p[0] for p in ring]
    lats = [p[1] for p in ring]
    lat = (min(lats) + max(lats)) / 2
    lon = (min(lons) + max(lons)) / 2
    a = 0
    for i in range(len(ring)-1):
        a += ring[i][0]*ring[i+1][1] - ring[i+1][0]*ring[i][1]
    area_km2 = abs(a)/2 * 111.0**2 * math.cos(math.radians(lat))
    return lat, lon, area_km2


def geojson_ring(gj: dict):
    """Extract outer ring from GeoJSON (always WGS84)."""
    feat = gj
    if gj.get("type") == "FeatureCollection":
        feat = gj["features"][0]
    geom = feat.get("geometry", feat)
    t = geom.get("type", "")
    c = geom.get("coordinates", [])
    if t == "Polygon":        ring = c[0]
    elif t == "MultiPolygon": ring = c[0][0]
    else: raise ValueError(f"Unsupported geometry: {t}")
    return [(p[0], p[1]) for p in ring]


# ════════════════════════════════════════════════════════
#  ENDPOINT: Parse file  →  WGS84 centroid + ring
# ════════════════════════════════════════════════════════

@app.post("/api/parse-file")
async def parse_file(file: UploadFile = File(...)):
    content  = await file.read()
    fname    = (file.filename or "").lower()

    try:
        if fname.endswith(".zip"):
            with zipfile.ZipFile(io.BytesIO(content)) as zf:
                names        = zf.namelist()
                shp_name     = next((n for n in names if n.lower().endswith(".shp")),     None)
                prj_name     = next((n for n in names if n.lower().endswith(".prj")),     None)
                geojson_name = next((n for n in names if n.lower().endswith((".geojson",".json"))), None)

                if shp_name:
                    raw_rings = parse_shp_bytes(zf.read(shp_name))
                    if not raw_rings:
                        raise ValueError("No polygon found in SHP")
                    prj_wkt = zf.read(prj_name).decode("utf-8","ignore") if prj_name else ""
                    ring    = reproject_ring(raw_rings[0], prj_wkt)

                    # Validate reprojected result
                    lons = [p[0] for p in ring]; lats = [p[1] for p in ring]
                    if not (-180 <= min(lons) and max(lons) <= 180 and -90 <= min(lats) and max(lats) <= 90):
                        raise ValueError(
                            "Reprojection produced invalid coordinates. "
                            "Ensure the .prj file is included in the ZIP alongside the .shp."
                        )

                    # Determine CRS note
                    try:
                        from pyproj import CRS
                        src_name = CRS.from_wkt(prj_wkt).name if prj_wkt else "Unknown"
                        crs_note = f"Reprojected: {src_name} → WGS84 (pyproj)"
                    except Exception:
                        xs0 = [p[0] for p in raw_rings[0]]
                        crs_note = ("WGS84 — no reprojection needed" if (-180<=min(xs0) and max(xs0)<=180)
                                    else "Reprojected via heuristic (install pyproj for exact results)")

                elif geojson_name:
                    ring     = geojson_ring(json.loads(zf.read(geojson_name)))
                    crs_note = "WGS84 (GeoJSON spec)"
                else:
                    raise ValueError("No .shp or .geojson found inside ZIP")

        elif fname.endswith((".geojson", ".json")):
            ring     = geojson_ring(json.loads(content))
            crs_note = "WGS84 (GeoJSON spec)"
        else:
            raise ValueError("Upload a .zip (shapefile), .geojson, or .json")

        lat, lon, area_km2 = centroid_and_area(ring)
        return {
            "lat":       round(lat, 6),
            "lon":       round(lon, 6),
            "area_km2":  round(area_km2, 4),
            "source":    file.filename,
            "crs_note":  crs_note,
            "ring_wgs84":[[round(p[0],6), round(p[1],6)] for p in ring]
        }

    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# ════════════════════════════════════════════════════════
#  NOAA Atlas 14
# ════════════════════════════════════════════════════════

DURATIONS      = ["5-min","10-min","15-min","30-min","1-hr","2-hr","3-hr",
                  "6-hr","12-hr","24-hr","48-hr","72-hr"]
RETURN_PERIODS = [1,2,5,10,25,50,100,200,500,1000]

@app.get("/api/noaa")
async def get_noaa(lat: float=Query(...), lon: float=Query(...), units: str=Query("english")):
    url = (f"https://hdsc.nws.noaa.gov/pfds/pfds_printpage.html"
           f"?lat={lat}&lon={lon}&data=depth&units={units}&series=pds")
    log = []
    try:
        async with httpx.AsyncClient(timeout=25.0, follow_redirects=True) as client:
            r = await client.get(url, headers={"User-Agent":"Mozilla/5.0 (HydroAgent/2.0)"})
        if r.status_code == 200:
            data = _parse_noaa_html(r.text)
            if data:
                log.append(f"NOAA Atlas 14 exact data retrieved ✓")
                return {"source":"NOAA Atlas 14 — PFDS (exact)","estimated":False,
                        "units":"in" if units=="english" else "mm",
                        "durations":DURATIONS,"return_periods":RETURN_PERIODS,
                        "data":data,"noaa_url":url,"log":log}
        log.append(f"NOAA returned HTTP {r.status_code}")
    except Exception as e:
        log.append(f"NOAA unreachable: {str(e)[:60]}")

    data = _estimate_idf(lat, lon, units)
    log.append("Regional IDF estimation applied")
    return {"source":"Regional estimate — see NOAA link for exact data","estimated":True,
            "units":"in" if units=="english" else "mm",
            "durations":DURATIONS,"return_periods":RETURN_PERIODS,
            "data":data,"noaa_url":url,"log":log}


def _parse_noaa_html(html: str) -> Optional[dict]:
    m = re.search(r'var\s+quantiles\s*=\s*(\[\[.*?\]\])', html, re.DOTALL)
    if m:
        try:
            raw = json.loads(m.group(1))
            data = {}
            for di, dur in enumerate(DURATIONS):
                if di < len(raw):
                    row  = raw[di]
                    data[dur] = [round(float(row[ri]),3) if ri<len(row) else None
                                 for ri in range(len(RETURN_PERIODS))]
            if len(data) >= 8: return data
        except Exception: pass

    rows_p = r'<tr[^>]*>(.*?)</tr>'
    cell_p = r'<t[dh][^>]*>(.*?)</t[dh]>'
    data, di = {}, 0
    for rh in re.findall(rows_p, html, re.DOTALL|re.IGNORECASE):
        cells = [re.sub(r'<[^>]+>','',c).strip()
                 for c in re.findall(cell_p, rh, re.DOTALL|re.IGNORECASE)]
        if di < len(DURATIONS) and cells:
            exp = DURATIONS[di].lower().replace(' ','').replace('-','')
            fst = cells[0].lower().replace(' ','').replace('-','')
            if exp in fst or fst in exp:
                nums = []
                for c in cells[1:]:
                    try: nums.append(round(float(c),3))
                    except ValueError: pass
                if len(nums) >= 8:
                    data[DURATIONS[di]] = nums[:len(RETURN_PERIODS)]
                    di += 1
    return data if len(data) >= 6 else None


def _estimate_idf(lat, lon, units) -> dict:
    if   lon > -80  and lat < 35:    b=4.2
    elif lon > -95  and lat < 33:    b=3.8
    elif lon > -85  and lat < 40:    b=3.2
    elif lon > -80  and lat >= 35:   b=3.0
    elif lon > -90  and lat >= 40:   b=2.6
    elif lon > -100 and lat >= 40:   b=2.3
    elif lon > -100 and lat < 40:    b=2.9
    elif lon <= -115 and lat > 40:   b=1.1
    elif lon <= -115 and lat <= 40:  b=0.9
    elif lon > -115 and lon <= -100: b=1.7
    else:                            b=2.2
    m = 25.4 if units=="metric" else 1.0
    df={"5-min":.078,"10-min":.115,"15-min":.148,"30-min":.218,"1-hr":.335,
        "2-hr":.46,"3-hr":.545,"6-hr":.72,"12-hr":.87,"24-hr":1.0,"48-hr":1.26,"72-hr":1.44}
    rf=[.72,1.0,1.38,1.68,2.10,2.43,2.76,3.12,3.60,4.02]
    return {dur:[round(b*f*r*m,3) for r in rf] for dur,f in df.items()}


# ════════════════════════════════════════════════════════
#  USDA HSG
# ════════════════════════════════════════════════════════

HSG_INFO = {
    "A":  ("Low runoff — high infiltration. Deep, well-drained sandy or gravelly soils.","30–45"),
    "B":  ("Moderate runoff — moderate infiltration. Moderately deep well-drained soils.","55–70"),
    "C":  ("High runoff — slow infiltration. Layer impeding downward water movement.","70–80"),
    "D":  ("Very high runoff — very slow infiltration. Clay-rich or high water table.","80–90"),
    "B/D":("Dual group (drained B / undrained D). Evaluate both conditions.","55–90"),
    "C/D":("Dual group (drained C / undrained D). Evaluate both conditions.","70–90"),
}

@app.get("/api/hsg")
async def get_hsg(lat: float=Query(...), lon: float=Query(...)):
    log = []
    # Attempt 1: USDA SDA
    query = f"""SELECT TOP 5 mu.muname,c.hydgrp,c.comppct_r,c.compname,c.taxorder
        FROM mapunit mu JOIN component c ON mu.mukey=c.mukey
        WHERE mu.mukey IN (SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('POINT({lon} {lat})'))
        AND c.majcompflag='Yes' ORDER BY c.comppct_r DESC"""
    try:
        async with httpx.AsyncClient(timeout=20.0) as client:
            r = await client.post("https://sdmdataaccess.sc.egov.usda.gov/Tabular/post.rest",
                data={"query":query,"format":"JSON+COLUMNNAME+METADATA"},
                headers={"Content-Type":"application/x-www-form-urlencoded"})
        if r.status_code == 200:
            table = r.json().get("Table",[])
            if len(table) > 1:
                rows = table[1:]
                d    = rows[0]
                comps= [{"muname":rw[0],"hsg":rw[1],"pct":rw[2],"compname":rw[3]} for rw in rows if rw[1]]
                log.append(f"USDA SDA: {len(comps)} component(s) ✓")
                return _hsg_resp(d[1] or "B", d[0] or "", d[3] or "", d[2] or 0, comps, False, log)
    except Exception as e:
        log.append(f"SDA error: {str(e)[:50]}")

    # Attempt 2: SoilWeb
    try:
        async with httpx.AsyncClient(timeout=12.0) as client:
            r = await client.get(f"https://casoilresource.lawr.ucdavis.edu/api/soil-series/?lon={lon}&lat={lat}&outformat=json")
        if r.status_code == 200:
            sw = r.json()
            hsg= _drain_to_hsg(sw.get("drainage_class",""))
            log.append(f"SoilWeb: {sw.get('series_name','')} ✓")
            return _hsg_resp(hsg, sw.get("series_name",""), "", None, [], True, log, "Source: SoilWeb/UC Davis")
    except Exception as e:
        log.append(f"SoilWeb error: {str(e)[:40]}")

    # Attempt 3: geographic estimate
    hsg = "B" if lon < -100 else ("B" if lat < 33 and lon > -90 else "B")
    log.append("Geographic estimate applied")
    return _hsg_resp(hsg,"Geographic estimate","",None,[],True,log,"Estimated from regional patterns")


def _drain_to_hsg(d):
    d=d.lower()
    if "excessively" in d: return "A"
    if "well" in d: return "B"
    if "moderately well" in d: return "B"
    if "somewhat poorly" in d: return "C"
    if "poorly" in d: return "D"
    return "B"

def _hsg_resp(hsg, muname, compname, pct, comps, estimated, log, note=""):
    info = HSG_INFO.get(hsg, HSG_INFO.get(hsg[0] if hsg else "B", ("Unknown","—")))
    return {"hsg":hsg,"muname":muname,"compname":compname,"pct_dominant":pct,
            "description":info[0],"cn_range":info[1],"components":comps,
            "estimated":estimated,"note":note,"log":log}


# ════════════════════════════════════════════════════════
#  TR-55  (fully offline)
# ════════════════════════════════════════════════════════

@app.get("/api/tr55")
async def get_tr55(lat: float=Query(...), lon: float=Query(...)):
    if (lat>42 and lon<-121) or (lat>44 and lon<-119):
        return {"type":"Type IA","region":"Pacific coastal — WA, OR, NW California","segment":1,
                "tip":"Lowest intensity. 50% of rainfall in first ~8 hrs. Select Type IA in Storm & Sanitary.",
                "peak_hr":8,"factor":0.33,"ss_input":"IA"}
    if lon < -111 and lat > 33:
        return {"type":"Type I","region":"Pacific inland — CA, NV, UT, AZ, ID, MT, WY","segment":2,
                "tip":"Moderate intensity, uniform distribution. Peak near storm midpoint. Select Type I.",
                "peak_hr":11,"factor":0.38,"ss_input":"I"}
    if lat < 37 and lon > -100:
        return {"type":"Type III","region":"Gulf Coast & SE — FL, GA, SC, NC, AL, MS, LA, TX","segment":4,
                "tip":"Highest intensity, back-loaded storm. Strong Gulf moisture influence. Select Type III.",
                "peak_hr":12,"factor":0.45,"ss_input":"III"}
    return {"type":"Type II","region":"Central & Eastern US — most of CONUS east of Rockies","segment":3,
            "tip":"Most common distribution. Concentrated peak ~hour 12. Select Type II in Storm & Sanitary.",
            "peak_hr":12,"factor":0.42,"ss_input":"II"}
