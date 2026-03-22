# HydroAgent v2 — Hydrologic Modeling Assistant

Upload a watershed shapefile → get NOAA Atlas 14 IDF table, HSG soil group,
TR-55 rainfall type, and a live map preview of your polygon — automatically.

**v2 fix:** Automatic reprojection of any CRS (State Plane, UTM, etc.) → WGS84
using pyproj. This was the root cause of the Ohio → New Jersey bug.

---

## WHAT YOU NEED (one time only)

1. Free **GitHub** account → https://github.com/signup
2. Free **Render.com** account → https://render.com (sign up with GitHub — 30 sec)

---

## DEPLOY IN 5 STEPS

### Step 1 — Upload to GitHub
1. Go to https://github.com → click **New repository** → name it `hydroagent`
2. Click **uploading an existing file**
3. Upload ALL files from this ZIP preserving the folder structure:
   ```
   main.py
   requirements.txt
   render.yaml
   README.md
   static/
     index.html
   ```
4. Click **Commit changes**

### Step 2 — Deploy on Render
1. Go to https://render.com → **New +** → **Web Service**
2. Click **Connect a repository** → select `hydroagent`
3. Render auto-detects settings from `render.yaml` — leave everything as-is
4. Click **Create Web Service**
5. Wait ~3 minutes for the build (pyproj takes a moment to compile)

### Step 3 — Get your URL
Render gives you a URL like `https://hydroagent.onrender.com`
Bookmark it — accessible from anywhere in the world.

### Step 4 — Use the app
1. Open the URL in any browser
2. Upload your watershed `.zip` (must include .shp + .dbf + .prj + .shx)
3. The map preview confirms your polygon is correctly georeferenced
4. Click **⚡ Run analysis**
5. Get IDF table, HSG type, and TR-55 rainfall type
6. Click **Copy summary** to paste into Storm & Sanitary Analysis

### Step 5 — Update the app
Re-upload changed files to GitHub → Render redeploys automatically in ~2 min.

---

## WHY v2 FIXES THE OHIO → NJ BUG

Civil 3D and AutoCAD Map export shapefiles in the project CRS
(State Plane, UTM, etc.) with coordinates in feet or meters.
When the previous version read those coordinates as if they were
lat/lon degrees, it placed the polygon in the wrong US state.

v2 reads the `.prj` file from the ZIP and uses **pyproj** to
reproject every coordinate point to WGS84 (lat/lon degrees) before
computing the centroid or displaying the polygon on the map.

The map preview lets you visually confirm the polygon landed in the
right place before running the analysis.

---

## FILE FORMATS

| Format | Notes |
|--------|-------|
| `.zip` | Include .shp + .dbf + .prj + .shx. The .prj is critical for reprojection. |
| `.geojson` | Always WGS84 per spec — no reprojection needed. |
| `.json` | Same as GeoJSON. |

---

## DATA SOURCES

| Source | Returns | Method |
|--------|---------|--------|
| NOAA Atlas 14 | Full IDF table (12 dur × 10 RP) | PFDS direct + regional fallback |
| USDA Web Soil Survey | Dominant HSG (A/B/C/D), soil series | SDA REST → SoilWeb → geographic est. |
| TR-55 (offline) | Storm type I/IA/II/III + S&S input string | Built-in coordinate lookup |

---

## NOTES

- **Free tier (Render):** Server sleeps after 15 min idle. First load after idle
  takes ~30s to wake up — normal behavior. Paid plan ($7/mo) keeps it always-on.
- **pyproj:** Required for exact reprojection. Render installs it automatically
  from requirements.txt. If pyproj install fails for any reason, the app falls
  back to a heuristic reprojection (accurate to ~500m, sufficient for all queries).
