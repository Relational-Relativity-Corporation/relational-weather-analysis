"""
run_analysis.py -- Execute ABRCE Weather Cell Analysis

Simulation build protocol:
    1. Fetch data and verify it looks physically reasonable
    2. Build graph and verify topology
    3. Run operators on a single snapshot, inspect raw output
    4. Run full temporal evolution
    5. Extract signals
    6. Report
"""

import json
import sys
from pathlib import Path
from datetime import datetime

from fetch_metar import (
    STATIONS, EDGES, COUPLING_WEIGHT,
    fetch_iem_data, fetch_buoy_data, build_snapshots, save_data
)
from abrce_graph_ops import (
    StationGraph, operator_a, operator_b, operator_r, operator_c,
    compute_rho, operator_e, compute_derived_fields,
    run_full_analysis, track_evolution, extract_precip_signals,
    ANALYSIS_VARS, DERIVED_FIELDS
)


def verify_data(observations):
    """Step 1: Verify raw data looks physically reasonable."""
    print("\n" + "=" * 60)
    print("STEP 1: Data Verification")
    print("=" * 60)

    if not observations:
        print("ERROR: No observations retrieved. Cannot proceed.")
        return False

    range_checks = {
        "tmpf": (-20, 130, "Temperature F"),
        "dwpf": (-40, 100, "Dewpoint F"),
        "relh": (0, 100, "Relative humidity %"),
        "drct": (0, 360, "Wind direction deg"),
        "sknt": (0, 100, "Wind speed kt"),
        "alti": (28.0, 31.0, "Altimeter inHg"),
    }

    issues = []
    for obs in observations:
        for var, (lo, hi, desc) in range_checks.items():
            val = obs.get(var)
            if val is not None and (val < lo or val > hi):
                issues.append(f"  {obs['station']} {obs['timestamp']}: "
                             f"{desc} = {val} (expected {lo}-{hi})")

    if issues:
        print(f"WARNING: {len(issues)} values outside expected ranges:")
        for iss in issues[:10]:
            print(iss)
        if len(issues) > 10:
            print(f"  ... and {len(issues)-10} more")
    else:
        print("All values within expected physical ranges.")

    by_station = {}
    for obs in observations:
        by_station.setdefault(obs["station"], []).append(obs)

    print("\nStation coverage:")
    for s in sorted(STATIONS.keys()):
        obs_list = by_station.get(s, [])
        if obs_list:
            times = sorted(o["timestamp"] for o in obs_list)
            print(f"  {s}: {len(obs_list)} obs, "
                  f"{times[0][:16]} to {times[-1][:16]}")
        else:
            print(f"  {s}: NO DATA")

    return True


def verify_graph(graph):
    """Step 2: Verify graph topology."""
    print("\n" + "=" * 60)
    print("STEP 2: Graph Topology Verification")
    print("=" * 60)

    print(f"Nodes: {graph.n}")
    print(f"Directed edges: {graph.edge_count()}")
    print(f"Undirected edges: {graph.edge_count() // 2}")

    print("\nAdjacency:")
    for s in graph.node_ids:
        nbs = graph.neighbors(s)
        nb_str = ", ".join(f"{nb}({w:.1f})" for nb, w in nbs)
        print(f"  {s} [{STATIONS[s]['name']}]: {nb_str}")

    print("\nEdge distances (km):")
    seen = set()
    for (a, b), d in sorted(graph.distances.items()):
        key = tuple(sorted([a, b]))
        if key not in seen:
            seen.add(key)
            print(f"  {a} <-> {b}: {d:.1f} km")

    return True


def verify_single_snapshot(graph, snapshots):
    """Step 3: Run operators on ONE snapshot, inspect raw output."""
    print("\n" + "=" * 60)
    print("STEP 3: Single Snapshot Operator Verification")
    print("=" * 60)

    best_snap = None
    best_count = 0
    for snap in snapshots:
        count = sum(1 for s in snap["stations"] if snap["stations"][s] is not None)
        if count > best_count:
            best_count = count
            best_snap = snap

    if best_snap is None or best_count < 3:
        print("ERROR: No snapshot with sufficient station coverage.")
        return False

    print(f"Using snapshot: {best_snap['timestamp']}")
    print(f"Stations reporting: {best_count}/{graph.n}")

    field = best_snap["stations"]

    print("\nRaw observations:")
    for s in graph.node_ids:
        obs = field.get(s)
        if obs:
            print(f"  {s}: T={obs.get('tmpf','M')}F  "
                  f"Td={obs.get('dwpf','M')}F  "
                  f"RH={obs.get('relh','M')}%  "
                  f"Wind={obs.get('drct','M')}deg/{obs.get('sknt','M')}kt  "
                  f"Alt={obs.get('alti','M')}inHg")
        else:
            print(f"  {s}: NO DATA")

    print("\n--- Operator A (temperature gradients) ---")
    a_tmpf = operator_a(graph, field, "tmpf")
    for (src, dst), g in sorted(a_tmpf.items()):
        if g is not None:
            direction = "warm->cool" if g > 0 else "cool->warm"
            print(f"  {src} -> {dst}: {g:+.2f}F ({direction})")

    print("\n--- Rho (circulation strength from T gradients) ---")
    rho = compute_rho(graph, a_tmpf)
    for s, r in sorted(rho.items()):
        print(f"  {s}: rho = {r:.4f}" if r is not None else f"  {s}: None")

    print("\n--- Operator B (accumulated T gradients) ---")
    b_tmpf = operator_b(graph, a_tmpf)
    for (src, dst), g in sorted(b_tmpf.items()):
        if g is not None:
            print(f"  {src} -> {dst}: {g:+.3f}")

    print("\n--- Operator R (circulation from T field) ---")
    r_tmpf = operator_r(graph, b_tmpf, rho)
    for s, v in sorted(r_tmpf.items()):
        if v is not None:
            indicator = ""
            if abs(v) > 0.5:
                indicator = " <- STRONG"
            elif abs(v) > 0.1:
                indicator = " <- moderate"
            print(f"  {s}: {v:+.4f}{indicator}")

    print("\n--- Operator E (full evolution, bounded) ---")
    e_tmpf = operator_e(graph, field, "tmpf")
    for s, v in sorted(e_tmpf.items()):
        if v is not None:
            print(f"  {s}: {v:+.4f}")

    print("\n--- Derived: Td depression ---")
    derived = compute_derived_fields(graph, field)
    for s in graph.node_ids:
        d = derived.get(s, {})
        td_dep = d.get("td_depression")
        theta = d.get("theta_approx")
        if td_dep is not None:
            saturation = "NEAR SATURATION" if td_dep < 5 else (
                "moist" if td_dep < 15 else "dry")
            theta_str = f"  theta_approx = {theta:.1f}F" if theta else ""
            print(f"  {s}: Td_dep = {td_dep:.1f}F ({saturation}){theta_str}")

    return True


def run_full_evolution(graph, snapshots, rho_base=0.3):
    """Step 4: Run full temporal evolution."""
    print("\n" + "=" * 60)
    print("STEP 4: Full Temporal Evolution")
    print("=" * 60)

    evolution = track_evolution(graph, snapshots, rho_base)
    print(f"Processed {len(evolution)} time steps")

    target = "KSBA"
    print(f"\nEvolution at {target} ({STATIONS[target]['name']}):")
    print(f"{'Time':>20s}  {'T_E':>8s}  {'Td_dep_E':>8s}  "
          f"{'P_E':>8s}  {'rho(T)':>8s}")
    print("-" * 60)

    for result in evolution:
        ts = result["timestamp"]
        try:
            t = datetime.fromisoformat(ts)
            time_str = t.strftime("%Y-%m-%d %H:%M")
        except (ValueError, TypeError):
            time_str = ts[:16]

        t_e = result["variables"].get("tmpf", {}).get("e_evolved", {}).get(target)
        td_e = result["variables"].get("td_depression", {}).get("e_evolved", {}).get(target)
        p_e = result["variables"].get("alti", {}).get("e_evolved", {}).get(target)
        rho_t = result["variables"].get("tmpf", {}).get("rho", {}).get(target)

        t_str = f"{t_e:>+8.4f}" if t_e is not None else f"{'M':>8s}"
        td_str = f"{td_e:>+8.4f}" if td_e is not None else f"{'M':>8s}"
        p_str = f"{p_e:>+8.4f}" if p_e is not None else f"{'M':>8s}"
        r_str = f"{rho_t:>8.4f}" if rho_t is not None else f"{'M':>8s}"
        print(f"{time_str:>20s}  {t_str}  {td_str}  {p_str}  {r_str}")

    return evolution


def run_signal_extraction(graph, evolution):
    """Step 5: Extract precipitation signals."""
    print("\n" + "=" * 60)
    print("STEP 5: Precipitation Signal Extraction")
    print("=" * 60)

    signals = extract_precip_signals(graph, evolution)

    if signals:
        latest_ts = sorted(signals.keys())[-1]
        latest = signals[latest_ts]

        print(f"\nLatest signals ({latest_ts[:16]}):")
        print(f"{'Station':>8s}  {'Composite':>10s}  {'Moisture':>10s}  "
              f"{'Convergence':>12s}  {'Circulation':>12s}")
        print("-" * 60)

        for s in sorted(latest.keys()):
            sig = latest[s]
            comp = sig["composite"]
            moist = sig["components"].get("moisture")
            conv = sig["components"].get("moisture_convergence")
            circ = sig["components"].get("moisture_circulation")

            m_str = f"{moist:>10.4f}" if moist is not None else f"{'M':>10s}"
            c_str = f"{conv:>12.4f}" if conv is not None else f"{'M':>12s}"
            ci_str = f"{circ:>12.4f}" if circ is not None else f"{'M':>12s}"
            print(f"{s:>8s}  {comp:>10.4f}  {m_str}  {c_str}  {ci_str}")

        print("\nPrecipitation assessment (ABRCE structural signal):")
        for s in sorted(latest.keys()):
            sig = latest[s]
            comp = sig["composite"]
            if comp > 0.5:
                level = "HIGH -- multiple convergent indicators"
            elif comp > 0.3:
                level = "MODERATE -- some structural support"
            elif comp > 0.1:
                level = "LOW -- weak or isolated indicators"
            else:
                level = "MINIMAL -- no structural signal"
            print(f"  {s} ({STATIONS[s]['name']}): {level}")

    return signals


def save_results(evolution, signals, output_dir="data"):
    """Save analysis results."""
    outdir = Path(output_dir)
    outdir.mkdir(exist_ok=True)

    summary = []
    for result in evolution:
        step = {"timestamp": result["timestamp"], "variables": {}}
        for var, data in result["variables"].items():
            step["variables"][var] = {
                "rho": data.get("rho", {}),
                "r_circulation": data.get("r_circulation", {}),
                "e_evolved": data.get("e_evolved", {}),
            }
        if "delta_e" in result:
            step["delta_e"] = result["delta_e"]
        summary.append(step)

    evo_path = outdir / "evolution_summary.json"
    with open(evo_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSaved evolution summary to {evo_path}")

    sig_path = outdir / "signals.json"
    with open(sig_path, "w") as f:
        json.dump(signals, f, indent=2, default=str)
    print(f"Saved signals to {sig_path}")


# ===================================================================
# MAIN
# ===================================================================

if __name__ == "__main__":
    hours = int(sys.argv[1]) if len(sys.argv) > 1 else 48
    rho_base = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3

    print("=" * 60)
    print("ABRCE Weather Cell Analysis")
    print("Graph Topology -- No Interpolation")
    print("=" * 60)
    print(f"\nDomain D: {len(STATIONS)} station nodes, "
          f"{len(EDGES)} coupling edges")
    print(f"M : O -> D : ASOS sensor responses -> numerical values")
    print(f"Operators: A -> B -> R -> C -> E on graph topology")
    print(f"rho_base: {rho_base}")
    print(f"Fetch window: {hours} hours")

    graph = StationGraph(STATIONS, EDGES)

    observations = fetch_iem_data(hours=hours)
    buoy_data = fetch_buoy_data()

    if not verify_data(observations):
        print("\nData verification failed. Stopping.")
        sys.exit(1)

    snapshots = build_snapshots(observations, interval_minutes=60)
    save_data(observations, buoy_data, snapshots)

    verify_graph(graph)

    if not verify_single_snapshot(graph, snapshots):
        print("\nSnapshot verification failed. Stopping.")
        sys.exit(1)

    evolution = run_full_evolution(graph, snapshots, rho_base)
    signals = run_signal_extraction(graph, evolution)
    save_results(evolution, signals)

    print("\n" + "=" * 60)
    print("Analysis complete.")
    print("=" * 60)
    print("\nDeclared open conditions:")
    print("  1. Graph edge weights (coupling strength) are structural")
    print("     declarations, not empirically calibrated. Sensitivity")
    print("     testing across weight values is needed.")
    print(f"  2. rho_base = {rho_base} -- not tuned. Try 0.1-0.5.")
    print("  3. Temporal resolution (1hr) may miss sub-hourly dynamics.")
    print("  4. B path alignment uses geographic bearing as axis proxy.")
    print("     Topographic channeling may make effective flow paths")
    print("     differ from straight-line bearings.")
    print("  5. No A-equivalent preprocessing exists at the sensor level.")
    print("     ASOS values contain absolute reference (station elevation,")
    print("     calibration offsets). The residual structural gap from")
    print("     section 3.10.6 of the Triad document applies.")
