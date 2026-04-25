"""
fetch_metar.py -- Raw ASOS/METAR Data Pipeline

M : O -> D declaration:
    O = ASOS sensor instrument responses (thermometer, hygrometer,
        anemometer, barometer, visibility sensor)
    M = sensor calibration + METAR encoding -> numerical values
    D = { station_id, timestamp, tmpf, dwpf, relh, drct, sknt,
          alti, mslp, p01i, vsby, skyc1, skyl1, metar_raw }

    M does NOT declare atmospheric conditions.
    M declares what instruments reported in response to conditions.
    The C -> O -> D chain is preserved: we do not collapse C -> D.

Data source: Iowa Environmental Mesonet (IEM) ASOS archive
    - Raw observations, minimal QC
    - No API key required
    - CSV format with lat/lon/elevation
"""

import csv
import io
import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

# ===================================================================
# STATION GRAPH -- Natural topology of SB County region
# ===================================================================

STATIONS = {
    "SBA": {
        "name": "Santa Barbara Municipal",
        "lat": 34.4241, "lon": -119.8425, "elev_ft": 10,
        "context": "coastal, south-facing, lee of Santa Ynez Mtns"
    },
    "SMX": {
        "name": "Santa Maria Public",
        "lat": 34.8941, "lon": -120.4522, "elev_ft": 261,
        "context": "inland valley, north of mountains"
    },
    "VBG": {
        "name": "Vandenberg SFB",
        "lat": 34.7334, "lon": -120.5840, "elev_ft": 369,
        "context": "coastal headland, exposed to NW flow"
    },
    "OXR": {
        "name": "Oxnard / Ventura",
        "lat": 34.2010, "lon": -119.2070, "elev_ft": 45,
        "context": "coastal plain, SE of Santa Barbara"
    },
    "SBP": {
        "name": "San Luis Obispo",
        "lat": 35.2381, "lon": -120.6441, "elev_ft": 212,
        "context": "inland valley, far north, different air mass boundary"
    },
    "LPC": {
        "name": "Lompoc",
        "lat": 34.6656, "lon": -120.4685, "elev_ft": 88,
        "context": "Santa Ynez River valley, inland of Vandenberg"
    },
    "IZA": {
        "name": "Santa Ynez",
        "lat": 34.6068, "lon": -119.8813, "elev_ft": 671,
        "context": "mountain valley, north side of Santa Ynez Mtns"
    },
}

EDGES = [
    ("VBG", "LPC", "strong",
     "VBG coastal flow enters Lompoc valley directly"),
    ("VBG", "SMX", "strong",
     "NW flow path, same side of mountains"),
    ("LPC", "SMX", "strong",
     "valley-to-valley, no barrier"),
    ("LPC", "IZA", "moderate",
     "connected via Santa Ynez River valley but elevation change"),
    ("IZA", "SBA", "moderate",
     "mountain pass coupling, San Marcos / Cold Spring"),
    ("SMX", "IZA", "moderate",
     "indirect, through Lompoc valley or over ridge"),
    ("SBA", "OXR", "moderate",
     "along coast but 65km apart, Rincon point between"),
    ("SMX", "SBP", "weak",
     "100km, different synoptic regime boundary"),
    ("VBG", "SBP", "weak",
     "far north, connected only in strong synoptic events"),
]

COUPLING_WEIGHT = {"strong": 1.0, "moderate": 0.6, "weak": 0.3}


def build_adjacency():
    """Build adjacency list from edge declarations."""
    adj = {s: [] for s in STATIONS}
    for a, b, ctype, note in EDGES:
        w = COUPLING_WEIGHT[ctype]
        adj[a].append({"neighbor": b, "weight": w, "type": ctype, "note": note})
        adj[b].append({"neighbor": a, "weight": w, "type": ctype, "note": note})
    return adj


# ===================================================================
# IEM DATA FETCH
# ===================================================================

IEM_BASE = "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py"

IEM_VARS = [
    "tmpf", "dwpf", "relh", "drct", "sknt", "alti", "mslp",
    "p01i", "vsby", "skyc1", "skyl1", "metar",
]


def fetch_iem_data(hours=48):
    """Fetch raw METAR observations from IEM for all declared stations."""
    station_params = "&".join(f"station={s}" for s in STATIONS)
    data_params = "&".join(f"data={v}" for v in IEM_VARS)

    url = (
        f"{IEM_BASE}?"
        f"{station_params}&"
        f"{data_params}&"
        f"hours={hours}&"
        f"tz=UTC&"
        f"format=onlycomma&"
        f"latlon=yes&"
        f"elev=yes&"
        f"report_type=3"
    )

    print(f"Fetching {hours}h of METAR data for {len(STATIONS)} stations...")
    print(f"URL: {url[:120]}...")

    try:
        req = Request(url, headers={"User-Agent": "ABRCE-Weather/1.0"})
        with urlopen(req, timeout=30) as resp:
            raw = resp.read().decode("utf-8")
    except (URLError, HTTPError) as e:
        print(f"ERROR: Failed to fetch IEM data: {e}")
        sys.exit(1)

    reader = csv.DictReader(io.StringIO(raw))
    observations = []

    for row in reader:
        obs = {
            "station": row.get("station", ""),
            "timestamp": row.get("valid", ""),
            "lat": _safe_float(row.get("lat")),
            "lon": _safe_float(row.get("lon")),
            "elev_m": _safe_float(row.get("elevation")),
            "tmpf": _safe_float(row.get("tmpf")),
            "dwpf": _safe_float(row.get("dwpf")),
            "relh": _safe_float(row.get("relh")),
            "drct": _safe_float(row.get("drct")),
            "sknt": _safe_float(row.get("sknt")),
            "alti": _safe_float(row.get("alti")),
            "mslp": _safe_float(row.get("mslp")),
            "p01i": _safe_float(row.get("p01i")),
            "vsby": _safe_float(row.get("vsby")),
            "skyc1": row.get("skyc1", "M"),
            "skyl1": _safe_float(row.get("skyl1")),
            "metar": row.get("metar", ""),
        }
        observations.append(obs)

    print(f"Retrieved {len(observations)} observations")

    by_station = {}
    for obs in observations:
        s = obs["station"]
        by_station.setdefault(s, []).append(obs)

    for s in STATIONS:
        count = len(by_station.get(s, []))
        missing_vars = []
        if count > 0:
            sample = by_station[s][0]
            for v in ["tmpf", "dwpf", "relh", "drct", "sknt", "alti"]:
                if sample[v] is None:
                    missing_vars.append(v)
        status = f"{count} obs" + (f", missing: {missing_vars}" if missing_vars else "")
        print(f"  {s} ({STATIONS[s]['name']}): {status}")

    return observations


def _safe_float(val):
    """Convert to float, return None for missing."""
    if val is None or val == "" or val == "M":
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


# ===================================================================
# NDBC BUOY DATA
# ===================================================================

BUOYS = {
    "46054": {
        "name": "W. Santa Barbara Channel",
        "lat": 34.274, "lon": -120.459,
        "context": "offshore, upwind of SB in NW flow"
    },
    "46011": {
        "name": "Santa Maria (offshore)",
        "lat": 34.868, "lon": -120.857,
        "context": "offshore, upwind of VBG/SMX"
    },
}

NDBC_BASE = "https://www.ndbc.noaa.gov/data/realtime2"


def fetch_buoy_data():
    """Fetch recent buoy observations from NDBC."""
    buoy_data = {}

    for bid, info in BUOYS.items():
        url = f"{NDBC_BASE}/{bid}.txt"
        print(f"Fetching buoy {bid} ({info['name']})...")

        try:
            req = Request(url, headers={"User-Agent": "ABRCE-Weather/1.0"})
            with urlopen(req, timeout=15) as resp:
                raw = resp.read().decode("utf-8")
        except (URLError, HTTPError) as e:
            print(f"  WARNING: Failed to fetch buoy {bid}: {e}")
            buoy_data[bid] = []
            continue

        lines = raw.strip().split("\n")
        if len(lines) < 3:
            print(f"  WARNING: No data for buoy {bid}")
            buoy_data[bid] = []
            continue

        observations = []
        for line in lines[2:]:
            parts = line.split()
            if len(parts) < 15:
                continue
            try:
                obs = {
                    "buoy_id": bid,
                    "timestamp": f"{parts[0]}-{parts[1]}-{parts[2]}T{parts[3]}:{parts[4]}:00Z",
                    "lat": info["lat"],
                    "lon": info["lon"],
                    "wdir": _safe_float_ndbc(parts[5]),
                    "wspd": _safe_float_ndbc(parts[6]),
                    "gust": _safe_float_ndbc(parts[7]),
                    "pres": _safe_float_ndbc(parts[12]),
                    "atmp": _safe_float_ndbc(parts[13]),
                    "wtmp": _safe_float_ndbc(parts[14]),
                    "dewp": _safe_float_ndbc(parts[15]) if len(parts) > 15 else None,
                }
                observations.append(obs)
            except (IndexError, ValueError):
                continue

        buoy_data[bid] = observations
        print(f"  Retrieved {len(observations)} observations")

    return buoy_data


def _safe_float_ndbc(val):
    """NDBC uses 99, 999, 9999 etc. for missing values."""
    if val is None:
        return None
    try:
        v = float(val)
        if v in (99.0, 99.00, 999.0, 9999.0):
            return None
        return v
    except (ValueError, TypeError):
        return None


# ===================================================================
# SNAPSHOT BUILDER
# ===================================================================

def build_snapshots(observations, interval_minutes=60):
    """Align observations to regular time intervals."""
    if not observations:
        return []

    timestamps = []
    for obs in observations:
        try:
            t = datetime.fromisoformat(obs["timestamp"].replace("Z", "+00:00"))
            timestamps.append(t)
        except (ValueError, TypeError):
            continue

    if not timestamps:
        return []

    t_min = min(timestamps)
    t_max = max(timestamps)

    dt = timedelta(minutes=interval_minutes)
    t_start = t_min.replace(minute=0, second=0, microsecond=0)
    half = timedelta(minutes=interval_minutes / 2)

    by_station = {}
    for obs in observations:
        s = obs["station"]
        try:
            t = datetime.fromisoformat(obs["timestamp"].replace("Z", "+00:00"))
            by_station.setdefault(s, []).append((t, obs))
        except (ValueError, TypeError):
            continue

    for s in by_station:
        by_station[s].sort(key=lambda x: x[0])

    snapshots = []
    t_curr = t_start
    while t_curr <= t_max:
        snap = {"timestamp": t_curr.isoformat(), "stations": {}}

        for s in STATIONS:
            if s not in by_station:
                snap["stations"][s] = None
                continue

            best = None
            best_dt = half + timedelta(seconds=1)
            for t_obs, obs in by_station[s]:
                d = abs(t_obs - t_curr)
                if d < best_dt:
                    best_dt = d
                    best = obs

            if best_dt <= half:
                snap["stations"][s] = best
            else:
                snap["stations"][s] = None

        snapshots.append(snap)
        t_curr += dt

    total_slots = len(snapshots) * len(STATIONS)
    filled = sum(
        1 for snap in snapshots
        for s in snap["stations"]
        if snap["stations"][s] is not None
    )
    print(f"\nBuilt {len(snapshots)} snapshots ({interval_minutes}min intervals)")
    print(f"Coverage: {filled}/{total_slots} station-slots filled "
          f"({100*filled/total_slots:.1f}%)")

    return snapshots


# ===================================================================
# OUTPUT
# ===================================================================

def save_data(observations, buoy_data, snapshots, output_dir="data"):
    """Save raw data and snapshots to files."""
    outdir = Path(output_dir)
    outdir.mkdir(exist_ok=True)

    obs_path = outdir / "metar_raw.json"
    with open(obs_path, "w") as f:
        json.dump(observations, f, indent=2, default=str)
    print(f"\nSaved {len(observations)} raw observations to {obs_path}")

    buoy_path = outdir / "buoy_raw.json"
    with open(buoy_path, "w") as f:
        json.dump(buoy_data, f, indent=2, default=str)
    total_buoy = sum(len(v) for v in buoy_data.values())
    print(f"Saved {total_buoy} buoy observations to {buoy_path}")

    snap_path = outdir / "snapshots.json"
    with open(snap_path, "w") as f:
        json.dump(snapshots, f, indent=2, default=str)
    print(f"Saved {len(snapshots)} snapshots to {snap_path}")

    topo = {
        "stations": STATIONS,
        "edges": [
            {"from": a, "to": b, "type": t, "weight": COUPLING_WEIGHT[t], "note": n}
            for a, b, t, n in EDGES
        ],
        "adjacency": {
            s: [{"neighbor": e["neighbor"], "weight": e["weight"]}
                for e in adj_list]
            for s, adj_list in build_adjacency().items()
        }
    }
    topo_path = outdir / "graph_topology.json"
    with open(topo_path, "w") as f:
        json.dump(topo, f, indent=2)
    print(f"Saved graph topology to {topo_path}")


# ===================================================================
# MAIN
# ===================================================================

if __name__ == "__main__":
    hours = int(sys.argv[1]) if len(sys.argv) > 1 else 48

    print("=" * 60)
    print("ABRCE Weather Cell Analysis -- Data Pipeline")
    print("=" * 60)
    print(f"\nM : O -> D declaration:")
    print(f"  O = ASOS sensor responses at {len(STATIONS)} stations")
    print(f"  M = instrument calibration -> numerical METAR values")
    print(f"  D = scalar values on graph nodes (no interpolation)")
    print(f"  C -> O -> D chain preserved: no conditions assumed")
    print(f"\nFetching {hours} hours of data...\n")

    observations = fetch_iem_data(hours=hours)
    buoy_data = fetch_buoy_data()
    snapshots = build_snapshots(observations, interval_minutes=60)
    save_data(observations, buoy_data, snapshots)

    print("\n" + "=" * 60)
    print("Data pipeline complete. Next: abrce_graph_ops.py")
    print("=" * 60)
