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
    Handles: State Plane (feet/meters), UTM.
    Uses mid standard parallel for m_lon — accurate to ~200m for CONUS State Plane.
    """
    prj_up = prj_wkt.upper()

    is_feet = ("FOOT" in prj_up or "FEET" in prj_up or "US_SURVEY_FOOT" in prj_up)
    FEET_TO_M = 0.3048006096
    ring_m = [(p[0] * FEET_TO_M, p[1] * FEET_TO_M) if is_feet else p for p in ring]

    def _ex(pattern, text, default):
        m = re.search(pattern, text, re.IGNORECASE)
        return float(m.group(1)) if m else default

    cm  = _ex(r'CENTRAL_MERIDIAN[",\s]+(-?\d+\.?\d*)',      prj_wkt, -96.0)
    lo  = _ex(r'LATITUDE_OF_ORIGIN[",\s]+(\d+\.?\d*)',       prj_wkt,  40.0)
    fe  = _ex(r'FALSE_EASTING[",\s]+(\d+\.?\d*)',             prj_wkt,   0.0)
    fn  = _ex(r'FALSE_NORTHING[",\s]+(\d+\.?\d*)',            prj_wkt,   0.0)
    sp1 = _ex(r'STANDARD_PARALLEL_1[",\s]+(\d+\.?\d*)',       prj_wkt,   0.0)
    sp2 = _ex(r'STANDARD_PARALLEL_2[",\s]+(\d+\.?\d*)',       prj_wkt,   0.0)

    if is_feet:
        fe *= FEET_TO_M
        fn *= FEET_TO_M

    # Use mid standard parallel for lon scale (more accurate than lat_origin)
    ref_lat = (sp1 + sp2) / 2 if sp1 and sp2 else lo
    lat_rad = math.radians(ref_lat)
    m_lat   = 111320.0
    m_lon   = 111320.0 * math.cos(lat_rad)

    result = []
    for p in ring_m:
        lon_deg = cm + (p[0] - fe) / m_lon
        lat_deg = lo + (p[1] - fn) / m_lat
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

# Duration in hours — used to convert depth → intensity (depth / duration_hr)
DURATION_HR = {
    "5-min":  5/60,
    "10-min": 10/60,
    "15-min": 15/60,
    "30-min": 30/60,
    "1-hr":   1.0,
    "2-hr":   2.0,
    "3-hr":   3.0,
    "6-hr":   6.0,
    "12-hr":  12.0,
    "24-hr":  24.0,
    "48-hr":  48.0,
    "72-hr":  72.0,
}

def _depth_to_intensity(depth_data: dict) -> dict:
    """Convert depth table (in or mm) → intensity table (in/hr or mm/hr)."""
    intensity = {}
    for dur, values in depth_data.items():
        hr = DURATION_HR.get(dur, 1.0)
        intensity[dur] = [
            round(v / hr, 3) if v is not None else None
            for v in values
        ]
    return intensity


@app.get("/api/noaa")
async def get_noaa(lat: float=Query(...), lon: float=Query(...), units: str=Query("english")):
    url = (f"https://hdsc.nws.noaa.gov/pfds/pfds_printpage.html"
           f"?lat={lat}&lon={lon}&data=depth&units={units}&series=pds")
    log = []
    try:
        async with httpx.AsyncClient(timeout=25.0, follow_redirects=True) as client:
            r = await client.get(url, headers={"User-Agent":"Mozilla/5.0 (HydroAgent/2.0)"})
        if r.status_code == 200:
            depth_data = _parse_noaa_html(r.text)
            if depth_data:
                intensity_data = _depth_to_intensity(depth_data)
                log.append(f"NOAA Atlas 14 exact data retrieved ✓")
                return {"source":"NOAA Atlas 14 — PFDS (exact)","estimated":False,
                        "units":"in/hr" if units=="english" else "mm/hr",
                        "depth_units":"in" if units=="english" else "mm",
                        "durations":DURATIONS,"return_periods":RETURN_PERIODS,
                        "data":intensity_data,
                        "depth_data":depth_data,
                        "noaa_url":url,"log":log}
        log.append(f"NOAA returned HTTP {r.status_code}")
    except Exception as e:
        log.append(f"NOAA unreachable: {str(e)[:60]}")

    depth_data = _estimate_idf(lat, lon, units)
    intensity_data = _depth_to_intensity(depth_data)
    log.append("Regional IDF estimation applied")
    return {"source":"Regional estimate — see NOAA link for exact data","estimated":True,
            "units":"in/hr" if units=="english" else "mm/hr",
            "depth_units":"in" if units=="english" else "mm",
            "durations":DURATIONS,"return_periods":RETURN_PERIODS,
            "data":intensity_data,
            "depth_data":depth_data,
            "noaa_url":url,"log":log}


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
    "A":  ("Low runoff potential. High infiltration rate. Deep, well-drained sandy or gravelly soils. Transmission rate > 0.30 in/hr.", "30–45"),
    "B":  ("Moderately low runoff potential. Moderate infiltration rate. Moderately deep, well-drained soils with moderate to fine texture. Transmission 0.15–0.30 in/hr.", "55–70"),
    "C":  ("Moderately high runoff potential. Slow infiltration rate. Soils with a layer that impedes downward water movement or moderately fine texture. Transmission 0.05–0.15 in/hr.", "70–80"),
    "D":  ("High runoff potential. Very slow infiltration rate. Clays with high shrink-swell potential, soils with high water table, or shallow soils over impervious material. Transmission < 0.05 in/hr.", "80–90"),
    "B/D":("Dual hydrologic group. Drained condition = B. Undrained condition = D. Evaluate both for NRCS method.", "55–90"),
    "C/D":("Dual hydrologic group. Drained condition = C. Undrained condition = D. Evaluate both for NRCS method.", "70–90"),
}

@app.get("/api/hsg")
async def get_hsg(lat: float=Query(...), lon: float=Query(...)):
    log = []

    # Web Soil Survey PDF report URL (always include regardless of data source)
    wss_report_url = (
        f"https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx"
        f"?action=GetSoilReport&lat={lat}&lon={lon}"
    )
    # Simpler direct WSS URL that works as a link
    wss_url = "https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx"

    # Attempt 1: USDA SDA — extended query with more soil properties
    query = f"""SELECT TOP 10
        mu.muname, mu.mukey,
        c.hydgrp, c.comppct_r, c.compname, c.taxorder, c.taxsuborder,
        c.drainagecl, c.hydricrating,
        c.taxclname
    FROM mapunit mu
    JOIN component c ON mu.mukey = c.mukey
    WHERE mu.mukey IN (
        SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('POINT({lon} {lat})')
    )
    AND c.majcompflag = 'Yes'
    ORDER BY c.comppct_r DESC"""

    try:
        async with httpx.AsyncClient(timeout=20.0) as client:
            r = await client.post(
                "https://sdmdataaccess.sc.egov.usda.gov/Tabular/post.rest",
                data={"query": query, "format": "JSON+COLUMNNAME+METADATA"},
                headers={"Content-Type": "application/x-www-form-urlencoded"}
            )
        if r.status_code == 200:
            table = r.json().get("Table", [])
            if len(table) > 1:
                rows = table[1:]
                d = rows[0]
                # cols: muname, mukey, hydgrp, comppct_r, compname, taxorder,
                #       taxsuborder, drainagecl, hydricrating, taxclname
                hsg        = d[2] or "B"
                muname     = d[0] or ""
                mukey      = d[1] or ""
                compname   = d[4] or ""
                taxorder   = d[5] or ""
                taxsuborder= d[6] or ""
                drainage   = d[7] or ""
                hydric     = d[8] or ""
                taxclass   = d[9] or ""
                pct        = d[3] or 0

                comps = [{
                    "muname":   rw[0], "hsg": rw[2], "pct": rw[3],
                    "compname": rw[4], "taxorder": rw[5],
                    "drainage": rw[7], "hydric": rw[8]
                } for rw in rows if rw[2]]

                # Build WSS link with mukey for direct report
                wss_direct = f"https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx"
                log.append(f"USDA SDA: {len(comps)} component(s) found ✓")

                return _hsg_resp(
                    hsg, muname, compname, taxorder, taxsuborder,
                    drainage, hydric, taxclass, pct, comps,
                    False, log, "", wss_direct, lat, lon
                )
    except Exception as e:
        log.append(f"SDA error: {str(e)[:60]} — trying SoilWeb")

    # Attempt 2: SoilWeb (UC Davis)
    try:
        async with httpx.AsyncClient(timeout=12.0) as client:
            r = await client.get(
                f"https://casoilresource.lawr.ucdavis.edu/api/soil-series/"
                f"?lon={lon}&lat={lat}&outformat=json"
            )
        if r.status_code == 200:
            sw  = r.json()
            hsg = _drain_to_hsg(sw.get("drainage_class", ""))
            log.append(f"SoilWeb: series={sw.get('series_name','')} ✓")
            return _hsg_resp(
                hsg, sw.get("series_name",""), sw.get("series_name",""),
                "", "", sw.get("drainage_class",""), "", "", None,
                [], True, log, "Source: SoilWeb / UC Davis", wss_url, lat, lon
            )
    except Exception as e:
        log.append(f"SoilWeb error: {str(e)[:40]} — using geographic estimate")

    # Attempt 3: geographic estimate
    hsg = "B"
    log.append("Geographic estimate applied")
    return _hsg_resp(
        hsg, "Geographic estimate", "", "", "", "", "", "", None,
        [], True, log, "Estimated — open Web Soil Survey for exact data",
        wss_url, lat, lon
    )


def _drain_to_hsg(d):
    d = d.lower()
    if "excessively" in d: return "A"
    if "somewhat excessively" in d: return "A"
    if "well" in d: return "B"
    if "moderately well" in d: return "B"
    if "somewhat poorly" in d: return "C"
    if "poorly" in d or "very poorly" in d: return "D"
    return "B"

def _hsg_resp(hsg, muname, compname, taxorder, taxsuborder,
              drainage, hydric, taxclass, pct, comps,
              estimated, log, note, wss_url, lat, lon):
    info = HSG_INFO.get(hsg, HSG_INFO.get(hsg[0] if hsg else "B", ("Unknown","—")))
    # Build direct WSS PDF bookmark URL
    wss_pdf_url = (
        f"https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx"
    )
    return {
        "hsg":          hsg,
        "muname":       muname,
        "compname":     compname,
        "taxorder":     taxorder,
        "taxsuborder":  taxsuborder,
        "drainage_class": drainage,
        "hydric_rating":  hydric,
        "tax_class":    taxclass,
        "pct_dominant": pct,
        "description":  info[0],
        "cn_range":     info[1],
        "components":   comps,
        "estimated":    estimated,
        "note":         note,
        "wss_url":      wss_url,
        "wss_pdf_url":  wss_pdf_url,
        "lat":          lat,
        "lon":          lon,
        "log":          log
    }


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


# ════════════════════════════════════════════════════════
#  SOIL REPORT — full SSURGO data → HTML/PDF
# ════════════════════════════════════════════════════════

@app.get("/api/soil-report")
async def soil_report(lat: float=Query(...), lon: float=Query(...)):
    """
    Fetch comprehensive soil data from USDA SDA and return
    a full HTML report (printable as PDF from browser).
    """
    from fastapi.responses import HTMLResponse
    from datetime import datetime

    # ── Query 1: Component data ──────────────────────────────────────────────
    q_comp = f"""
    SELECT
        mu.muname, mu.mukey, mu.musym,
        c.compname, c.comppct_r, c.majcompflag,
        c.hydgrp, c.drainagecl, c.hydricrating,
        c.taxorder, c.taxsuborder, c.taxgrtgroup, c.taxsubgrp,
        c.taxclname, c.slope_r, c.slope_l, c.slope_h,
        c.elev_r, c.aspectrep,
        c.tfact, c.wei, c.weg,
        c.nirrcapcl, c.nirrcapscl,
        c.irrcapcl,  c.irrcapscl
    FROM mapunit mu
    JOIN component c ON mu.mukey = c.mukey
    WHERE mu.mukey IN (
        SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('POINT({lon} {lat})')
    )
    ORDER BY c.comppct_r DESC"""

    # ── Query 2: Horizon data for dominant component ─────────────────────────
    q_horiz = f"""
    SELECT TOP 20
        c.compname, c.comppct_r,
        h.hzname, h.hzdept_r, h.hzdepb_r,
        h.texture, h.texdesc,
        h.sandtotal_r, h.silttotal_r, h.claytotal_r,
        h.om_r, h.ph1to1h2o_r,
        h.ksat_r, h.awc_r,
        h.dbthirdbar_r,
        h.cec7_r
    FROM mapunit mu
    JOIN component c ON mu.mukey = c.mukey
    JOIN chorizon h  ON c.cokey  = h.cokey
    WHERE mu.mukey IN (
        SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('POINT({lon} {lat})')
    )
    AND c.majcompflag = 'Yes'
    ORDER BY c.comppct_r DESC, h.hzdept_r ASC"""

    # ── Query 3: Land capability and interpretations ─────────────────────────
    q_interp = f"""
    SELECT TOP 10
        c.compname, c.comppct_r,
        ci.rulename, ci.interphrc, ci.interphr
    FROM mapunit mu
    JOIN component c    ON mu.mukey  = c.mukey
    JOIN cointerp ci    ON c.cokey   = ci.cokey
    WHERE mu.mukey IN (
        SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('POINT({lon} {lat})')
    )
    AND c.majcompflag = 'Yes'
    AND ci.mrulename IN (
        'NRCS Irrigation Suitability',
        'Hydrologic Soil Group',
        'Flooding Frequency Class',
        'Ponding Frequency Class',
        'Depth to Restrictive Layer',
        'Shrink-Swell Potential',
        'Corrosion of Steel',
        'Corrosion of Concrete'
    )
    ORDER BY c.comppct_r DESC, ci.rulename"""

    SDA = "https://sdmdataaccess.sc.egov.usda.gov/Tabular/post.rest"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}

    comp_rows, horiz_rows, interp_rows = [], [], []

    async with httpx.AsyncClient(timeout=30.0) as client:
        try:
            r = await client.post(SDA, data={"query":q_comp,"format":"JSON+COLUMNNAME+METADATA"}, headers=headers)
            if r.status_code == 200:
                t = r.json().get("Table", [])
                comp_rows = t[1:] if len(t) > 1 else []
        except Exception: pass

        try:
            r = await client.post(SDA, data={"query":q_horiz,"format":"JSON+COLUMNNAME+METADATA"}, headers=headers)
            if r.status_code == 200:
                t = r.json().get("Table", [])
                horiz_rows = t[1:] if len(t) > 1 else []
        except Exception: pass

        try:
            r = await client.post(SDA, data={"query":q_interp,"format":"JSON+COLUMNNAME+METADATA"}, headers=headers)
            if r.status_code == 200:
                t = r.json().get("Table", [])
                interp_rows = t[1:] if len(t) > 1 else []
        except Exception: pass

    # ── Build HTML report ────────────────────────────────────────────────────
    now   = datetime.now().strftime("%B %d, %Y  %H:%M")
    muname= comp_rows[0][0] if comp_rows else "Unknown"
    musym = comp_rows[0][2] if comp_rows else "—"

    def safe(v, suffix="", decimals=None):
        if v is None or v == "": return "—"
        if decimals is not None:
            try: return f"{float(v):.{decimals}f}{suffix}"
            except: pass
        return str(v) + suffix

    def hdr(cols):
        return "<tr>" + "".join(f"<th>{c}</th>" for c in cols) + "</tr>"

    def row(cells):
        return "<tr>" + "".join(f"<td>{safe(c)}</td>" for c in cells) + "</tr>"

    # Component table
    comp_html = ""
    if comp_rows:
        comp_html = f"""
        <table>
          {hdr(["Component","Symbol","%","HSG","Drainage","Hydric","Tax. Order","Tax. Class","LCC (irr.)","LCC (non-irr.)","Slope (%)"])}
          {"".join(row([r[3],r[2],safe(r[4],"%"),safe(r[6]),safe(r[7]),safe(r[8]),safe(r[9]),safe(r[13]),
                        f"{safe(r[22])}{safe(r[23])}",f"{safe(r[20])}{safe(r[21])}",
                        f"{safe(r[14],decimals=1)}–{safe(r[16],decimals=1)}"]) for r in comp_rows[:8])}
        </table>"""

    # Horizon table
    horiz_html = ""
    if horiz_rows:
        horiz_html = f"""
        <table>
          {hdr(["Horizon","Depth top (cm)","Depth bot (cm)","Texture","Sand %","Silt %","Clay %",
                "OM %","pH","Ksat (µm/s)","AWC (in/in)","Bulk density","CEC"])}
          {"".join(row([r[2],safe(r[3]),safe(r[4]),safe(r[5]),safe(r[8],decimals=1),safe(r[9],decimals=1),
                        safe(r[10],decimals=1),safe(r[11],decimals=2),safe(r[12],decimals=1),
                        safe(r[13],decimals=3),safe(r[14],decimals=3),safe(r[15],decimals=2),safe(r[16],decimals=1)]) for r in horiz_rows)}
        </table>"""

    # Interpretations table
    interp_html = ""
    if interp_rows:
        interp_html = f"""
        <table>
          {hdr(["Interpretation","Rating class","Rating value"])}
          {"".join(row([r[2],safe(r[3]),safe(r[4],decimals=2)]) for r in interp_rows)}
        </table>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Soil Report — {muname}</title>
<style>
  @media print {{
    body {{ margin: 0; }}
    .no-print {{ display: none; }}
    h2 {{ page-break-before: always; }}
    h2:first-of-type {{ page-break-before: avoid; }}
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: Arial, sans-serif; font-size: 11px; color: #1a1a1a;
          background: #fff; padding: 28px 36px; max-width: 960px; margin: 0 auto; }}
  .cover {{ border-bottom: 3px solid #1a5c38; margin-bottom: 20px; padding-bottom: 14px; }}
  .cover h1 {{ font-size: 20px; color: #1a5c38; margin-bottom: 4px; }}
  .cover .sub {{ font-size: 12px; color: #555; }}
  .meta-grid {{ display: grid; grid-template-columns: repeat(4,1fr); gap: 10px;
                background: #f4f8f5; border: 1px solid #c8ddd0; border-radius: 4px;
                padding: 12px; margin: 16px 0; }}
  .meta-item .ml {{ font-size: 9px; color: #777; text-transform: uppercase;
                    letter-spacing: .06em; margin-bottom: 2px; }}
  .meta-item .mv {{ font-size: 13px; font-weight: 700; color: #1a5c38; }}
  h2 {{ font-size: 13px; color: #fff; background: #1a5c38; padding: 6px 12px;
        margin: 20px 0 8px; border-radius: 3px; }}
  h3 {{ font-size: 11px; color: #1a5c38; margin: 12px 0 5px; font-weight: 700; }}
  table {{ width: 100%; border-collapse: collapse; margin-bottom: 10px; font-size: 10px; }}
  th {{ background: #2e7d52; color: #fff; padding: 5px 7px; text-align: left;
        font-weight: 600; white-space: nowrap; }}
  td {{ padding: 4px 7px; border-bottom: 1px solid #e0e8e3; vertical-align: top; }}
  tr:nth-child(even) td {{ background: #f6faf7; }}
  .hsg-box {{ display: inline-block; width: 54px; height: 54px; border-radius: 50%;
              border: 3px solid #1a5c38; text-align: center; line-height: 48px;
              font-size: 26px; font-weight: 700; margin-right: 16px; vertical-align: middle; }}
  .hA {{ background:#d4edda;color:#145a32; }}
  .hB {{ background:#cce5ff;color:#0c407a; }}
  .hC {{ background:#fff3cd;color:#7a5300; }}
  .hD {{ background:#f8d7da;color:#721c24; }}
  .hsg-info {{ display: inline-block; vertical-align: middle; max-width: 75%; }}
  .hsg-info p {{ font-size: 11px; color: #444; margin-top: 4px; line-height: 1.5; }}
  .note {{ font-size: 10px; color: #777; border-left: 3px solid #2e7d52; padding-left: 8px;
           margin: 10px 0; line-height: 1.6; }}
  .print-btn {{ position: fixed; top: 18px; right: 18px;
                background: #1a5c38; color: #fff; border: none;
                padding: 10px 18px; font-size: 12px; border-radius: 4px;
                cursor: pointer; font-family: Arial; }}
  .print-btn:hover {{ background: #145030; }}
  footer {{ margin-top: 28px; padding-top: 10px; border-top: 1px solid #ccc;
            font-size: 9px; color: #999; }}
</style>
</head>
<body>

<button class="print-btn no-print" onclick="window.print()">🖨 Print / Save PDF</button>

<div class="cover">
  <h1>Soil Survey Report — {muname}</h1>
  <div class="sub">Generated by HydroAgent · USDA SSURGO Database · {now}</div>
</div>

<div class="meta-grid">
  <div class="meta-item"><div class="ml">Latitude</div><div class="mv">{lat:.5f}°N</div></div>
  <div class="meta-item"><div class="ml">Longitude</div><div class="mv">{abs(lon):.5f}°W</div></div>
  <div class="meta-item"><div class="ml">Map Unit</div><div class="mv">{musym}</div></div>
  <div class="meta-item"><div class="ml">Report Date</div><div class="mv">{datetime.now().strftime('%Y-%m-%d')}</div></div>
</div>

<h2>1 — Hydrologic Soil Group (HSG)</h2>
{"".join([f'''
<div style="margin-bottom:14px">
  <span class="hsg-box h{r[6][0] if r[6] else "B"}">{r[6] or "—"}</span>
  <div class="hsg-info">
    <strong>{r[3]} ({safe(r[4])}%)</strong>
    <p>Drainage: {safe(r[7])} &nbsp;|&nbsp; Hydric: {safe(r[8])} &nbsp;|&nbsp; Tax. order: {safe(r[9])}</p>
    <p>Taxonomic class: {safe(r[13])}</p>
  </div>
</div>
''' for r in comp_rows[:3]]) if comp_rows else "<p>No HSG data retrieved.</p>"}

<div class="note">
  <strong>HSG Reference (USDA-NRCS TR-55):</strong><br>
  Group A — High infiltration rate (&gt;0.30 in/hr). Sandy, gravelly. CN: 30–45<br>
  Group B — Moderate infiltration (0.15–0.30 in/hr). Moderately textured. CN: 55–70<br>
  Group C — Slow infiltration (0.05–0.15 in/hr). Layer impeding movement. CN: 70–80<br>
  Group D — Very slow infiltration (&lt;0.05 in/hr). Clay-rich or high WT. CN: 80–90
</div>

<h2>2 — Soil Components</h2>
{comp_html or "<p>No component data retrieved.</p>"}

<h2>3 — Soil Horizons (Dominant Component)</h2>
{horiz_html or "<p>No horizon data retrieved.</p>"}
<div class="note">Ksat = Saturated Hydraulic Conductivity (µm/s) &nbsp;|&nbsp; AWC = Available Water Capacity (in/in) &nbsp;|&nbsp; OM = Organic Matter &nbsp;|&nbsp; CEC = Cation Exchange Capacity (meq/100g)</div>

<h2>4 — Soil Interpretations</h2>
{interp_html or "<p>No interpretation data retrieved.</p>"}

<h2>5 — Data Sources &amp; References</h2>
<div class="note">
  <strong>Database:</strong> USDA-NRCS SSURGO (Soil Survey Geographic Database), accessed via Soil Data Access (SDA) REST API.<br>
  <strong>Web Soil Survey:</strong> <a href="https://websoilsurvey.sc.egov.usda.gov">websoilsurvey.sc.egov.usda.gov</a><br>
  <strong>Coordinates:</strong> {lat:.5f}°N, {abs(lon):.5f}°W (WGS84)<br>
  <strong>Note:</strong> Values shown are representative (median) values from SSURGO. Low and high values may differ.
  Always verify against the official Web Soil Survey for final design.
</div>

<footer>HydroAgent · USDA SSURGO · Report generated {now} · Coordinates: {lat:.5f}°N, {abs(lon):.5f}°W</footer>
</body>
</html>"""

    return HTMLResponse(content=html, media_type="text/html")

