"""
abrce_graph_ops.py -- ABRCE Operators on Graph Topology

This module implements the ABRCE operator kernel on an irregular
graph rather than a regular lattice.

The structural principles are identical to operators.rs:
    E(x, rho) = C(R(B(A(x)), rho))
    A -> B -> R -> C -> E

What changes is the topology:
    - Lattice: cells on a periodic ring/torus, fixed neighbors
    - Graph: stations as nodes, edges as declared coupling paths

TOPOLOGY DECLARATION (per Verifier note):
    The graph topology is part of the M extension -- it encodes
    the spatial and atmospheric coupling structure through which
    observables relate. It is NOT part of operator invariance.

    Operators A, B, R, C are invariant: they compute differences,
    accumulate along paths, cross-couple, and bound -- regardless
    of the topology they operate on.

    The graph (node set, edge set, edge weights, bearing geometry)
    is declared infrastructure. Changing the graph changes what
    M maps onto D, not how operators process D.

    Consequence: results depend on graph declaration. This is
    analogous to how lattice results depend on grid resolution.
    Neither is an operator property -- both are M properties.

Operator A:
    Lattice: x[i] - x[i+1]  (forward neighbor on ring)
    Graph:   x[i] - x[j]    for each edge (i,j) in adjacency

    Output is per-EDGE, not per-cell. Each edge carries a
    directed difference. NO coupling weight applied at A --
    A is the pure primitive, extracting what M(o) reported
    at one station minus what M(o) reported at another.

Operator B:
    Lattice: gradient at i accumulates with gradient at i+1
             along the same axis
    Graph:   gradient on edge (i,j) accumulates with gradients
             on edges (j,k) that continue through j -- path
             accumulation through shared vertices

    Coupling weights applied HERE -- they encode how strongly
    gradients propagate through each vertex, which is an
    accumulation property (B's domain), not an extraction
    property (A's domain).

    Bearing alignment (cosine of bearing change) serves as
    the graph analog of "same axis." This is declared as part
    of topology geometry, not operator primitive.

Operator R:
    Lattice: cross-couples N-S gradients with E-W gradients
    Graph:   at each node, cross-couples gradients arriving
             from different edges -- the circulation is the
             asymmetry between incoming relational paths

Operator C:
    Unchanged. Pointwise: x / (1 + |x|)

Domain D is unchanged: { x in R^n | n < inf, |x[i]| < inf }
All quantifiers bounded over D. No claim beyond D.
"""

import math
from typing import Optional


# ===================================================================
# GRAPH STRUCTURE
# ===================================================================

class StationGraph:
    """
    Irregular graph of weather stations with declared edges.

    Nodes: station IDs (strings)
    Edges: (node_a, node_b, weight) -- weight from coupling declaration
    """

    def __init__(self, stations: dict, edges: list):
        self.stations = stations
        self.node_ids = sorted(stations.keys())
        self.node_index = {s: i for i, s in enumerate(self.node_ids)}
        self.n = len(self.node_ids)

        coupling_w = {"strong": 1.0, "moderate": 0.6, "weak": 0.3}
        self.adj = {s: [] for s in self.node_ids}
        self.edge_list = []

        for a, b, ctype, note in edges:
            w = coupling_w.get(ctype, 0.5)
            self.adj[a].append((b, w))
            self.adj[b].append((a, w))
            self.edge_list.append((a, b, w))
            self.edge_list.append((b, a, w))

        self.distances = {}
        for a, b, w in self.edge_list:
            if (a, b) not in self.distances:
                d = _haversine(
                    stations[a]["lat"], stations[a]["lon"],
                    stations[b]["lat"], stations[b]["lon"]
                )
                self.distances[(a, b)] = d
                self.distances[(b, a)] = d

    def neighbors(self, node_id):
        return self.adj.get(node_id, [])

    def edge_count(self):
        return len(self.edge_list)


def _haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    rlat1, rlon1 = math.radians(lat1), math.radians(lon1)
    rlat2, rlon2 = math.radians(lat2), math.radians(lon2)
    dlat = rlat2 - rlat1
    dlon = rlon2 - rlon1
    a = math.sin(dlat/2)**2 + math.cos(rlat1)*math.cos(rlat2)*math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(a))


# ===================================================================
# OPERATOR A -- Relational Gradient Extraction (Graph)
# ===================================================================

def operator_a(graph: StationGraph,
               field: dict,
               variable: str) -> dict:
    """
    A : field -> gradient structure (on edges)

    For each directed edge (i -> j), compute:
        a[(i,j)] = field[i] - field[j]

    Pure relational difference. No coupling weight, no scaling.
    A is the primitive -- it extracts exactly what M(o) reported
    at station i minus what M(o) reported at station j.

    Coupling geometry (edge weights) is applied downstream in B,
    where it belongs: weights encode how gradients propagate
    through the graph, which is an accumulation property,
    not an extraction property.

    Returns: dict of {(from, to): gradient_value}

    If either station has a None value (missing observation),
    that edge gradient is None. Missing data propagates honestly.
    """
    gradients = {}

    for a, b, w in graph.edge_list:
        val_a = _extract_value(field, a, variable)
        val_b = _extract_value(field, b, variable)

        if val_a is not None and val_b is not None:
            # Pure relational difference. No coupling weight here.
            # A is the primitive -- it extracts what M(o) reported
            # at station a minus what M(o) reported at station b.
            # Coupling geometry is B and R's domain.
            gradients[(a, b)] = val_a - val_b
        else:
            gradients[(a, b)] = None

    return gradients


def _extract_value(field: dict, station_id: str, variable: str):
    station_data = field.get(station_id)
    if station_data is None:
        return None
    val = station_data.get(variable)
    if val is None:
        return None
    return float(val)


# ===================================================================
# OPERATOR B -- Local Relational Accumulation (Graph)
# ===================================================================

def operator_b(graph: StationGraph,
               gradients: dict) -> dict:
    """
    B : gradient structure -> accumulated gradient structure

    On a graph, gradient on edge (i,j) accumulates with gradients
    on edges (j,k) that continue through j. Alignment is computed
    from geographic bearing change. Coupling weight applied here.

    alignment(i,j,k) = cos(bearing_change)
        +1 = same direction (accumulate)
         0 = perpendicular (R's job)
        -1 = reversal (opposing gradient)
    """
    accumulated = {}

    for (i, j), g_ij in gradients.items():
        if g_ij is None:
            accumulated[(i, j)] = None
            continue

        acc = g_ij

        # Find continuing edges from j (excluding back to i)
        for k, w_jk in graph.neighbors(j):
            if k == i:
                continue

            g_jk = gradients.get((j, k))
            if g_jk is None:
                continue

            align = _path_alignment(graph, i, j, k)

            # Coupling weight applied HERE, not in A.
            # w_jk encodes how strongly the atmospheric path
            # j->k couples to the path i->j. This is accumulation
            # geometry -- how far the gradient's relational reach
            # extends through this vertex -- not extraction.
            acc += g_jk * align * w_jk

        accumulated[(i, j)] = acc

    return accumulated


def _path_alignment(graph: StationGraph, i: str, j: str, k: str) -> float:
    """
    Compute alignment of path i->j->k based on geographic bearing.
    Returns cosine of bearing change.
    """
    si, sj, sk = graph.stations[i], graph.stations[j], graph.stations[k]

    bear_ij = _bearing(si["lat"], si["lon"], sj["lat"], sj["lon"])
    bear_jk = _bearing(sj["lat"], sj["lon"], sk["lat"], sk["lon"])

    delta = bear_jk - bear_ij
    while delta > 180:
        delta -= 360
    while delta < -180:
        delta += 360

    return math.cos(math.radians(delta))


def _bearing(lat1, lon1, lat2, lon2):
    rlat1 = math.radians(lat1)
    rlat2 = math.radians(lat2)
    dlon = math.radians(lon2 - lon1)
    x = math.sin(dlon) * math.cos(rlat2)
    y = (math.cos(rlat1) * math.sin(rlat2) -
         math.sin(rlat1) * math.cos(rlat2) * math.cos(dlon))
    return (math.degrees(math.atan2(x, y)) + 360) % 360


# ===================================================================
# COMPUTE RHO
# ===================================================================

def compute_rho(graph: StationGraph,
                gradients: dict,
                rho_base: float = 0.3) -> dict:
    """
    Per-node rho from A's gradient output.
    rho[i] = rho_base * max_grad / (1 + max_grad)
    """
    rho = {}

    for s in graph.node_ids:
        max_grad = 0.0
        has_data = False

        for neighbor, w in graph.neighbors(s):
            g = gradients.get((s, neighbor))
            if g is not None:
                has_data = True
                max_grad = max(max_grad, abs(g))

        if has_data:
            rho[s] = rho_base * max_grad / (1.0 + max_grad)
        else:
            rho[s] = None

    return rho


# ===================================================================
# OPERATOR R -- Antisymmetric Circulation (Graph)
# ===================================================================

def operator_r(graph: StationGraph,
               accumulated: dict,
               rho: dict) -> dict:
    """
    R : accumulated gradients x rho -> evolved field (on nodes)

    VERIFIER NOTE (structural tension, accepted as generalization):
        On the lattice, R has exactly 2 axes and 1 cross-term.
        On the graph, R computes cross-coupling across all edge
        pairs, weighted by sin(angle). This is "all-pair
        antisymmetry" rather than "axis-pair antisymmetry."

        This is accepted as a generalization because:
        (a) The graph has no natural axes -- imposing axis clusters
            would project lattice structure onto a topology that
            doesn't have it.
        (b) sin(angle) naturally suppresses near-parallel pairs
            and amplifies perpendicular pairs.
        (c) The sum (not average) preserves magnitude of strong
            perpendicular cross-terms.

    Returns: dict of {node_id: evolved_value}
    """
    evolved = {}

    for s in graph.node_ids:
        r = rho.get(s)
        if r is None:
            evolved[s] = None
            continue

        neighbors = graph.neighbors(s)
        if len(neighbors) < 2:
            evolved[s] = 0.0
            continue

        incoming = []
        for nb, w in neighbors:
            g = accumulated.get((nb, s))
            if g is not None:
                bear = _bearing(
                    graph.stations[nb]["lat"], graph.stations[nb]["lon"],
                    graph.stations[s]["lat"], graph.stations[s]["lon"]
                )
                incoming.append((bear, g, nb))

        if len(incoming) < 2:
            evolved[s] = 0.0
            continue

        circ_total = 0.0

        for idx_a in range(len(incoming)):
            for idx_b in range(idx_a + 1, len(incoming)):
                bear_a, g_a, _ = incoming[idx_a]
                bear_b, g_b, _ = incoming[idx_b]

                angle = bear_b - bear_a
                while angle > 180:
                    angle -= 360
                while angle < -180:
                    angle += 360

                cross = math.sin(math.radians(angle))
                circ_total += cross * (g_a - g_b)

        # circ_total is a SUM across pairs, not an average.
        # Rationale: averaging dilutes strong perpendicular
        # cross-terms with weak near-parallel pairs. sin(angle)
        # already suppresses near-parallel contributions, so
        # summing lets perpendicular terms dominate naturally.
        # On the lattice, R sums one cross-term (dz - dm) --
        # it does not average. Summing is the closer analog.
        evolved[s] = r * circ_total

    return evolved


# ===================================================================
# OPERATOR C -- Bounded Coherence
# ===================================================================

def operator_c(field: dict) -> dict:
    """C(x) = x / (1 + |x|). Output in (-1, 1). Pointwise."""
    bounded = {}
    for k, v in field.items():
        if v is not None:
            bounded[k] = v / (1.0 + abs(v))
        else:
            bounded[k] = None
    return bounded


# ===================================================================
# OPERATOR E -- Composite Evolution
# ===================================================================

def operator_e(graph: StationGraph,
               field: dict,
               variable: str,
               rho_base: float = 0.3) -> dict:
    """E(x, rho) = C(R(B(A(x)), rho(A(x))))"""
    a = operator_a(graph, field, variable)
    rho = compute_rho(graph, a, rho_base)
    b = operator_b(graph, a)
    r = operator_r(graph, b, rho)
    c = operator_c(r)
    return c


# ===================================================================
# MULTI-VARIABLE ANALYSIS
# ===================================================================

ANALYSIS_VARS = {
    "tmpf": "Temperature (F) -- thermometer response",
    "dwpf": "Dewpoint (F) -- moisture content indicator",
    "alti": "Altimeter setting (inHg) -- pressure at station",
    "sknt": "Wind speed (kt) -- anemometer response",
    "relh": "Relative humidity (%) -- derived T/Td ratio",
    "vsby": "Visibility (mi) -- transmissometer response",
}

DERIVED_FIELDS = {
    "td_depression": "T - Td: saturation deficit, key precipitation indicator",
    "theta_approx": "Potential temperature proxy: T + 0.0054 * elev_ft",
}


def compute_derived_fields(graph: StationGraph, field: dict) -> dict:
    """Compute derived observables from raw sensor data."""
    derived = {}

    for s in graph.node_ids:
        obs = field.get(s)
        if obs is None:
            derived[s] = {}
            continue

        d = {}
        tmpf = obs.get("tmpf")
        dwpf = obs.get("dwpf")
        elev = graph.stations[s].get("elev_ft", 0)

        if tmpf is not None and dwpf is not None:
            d["td_depression"] = tmpf - dwpf

        if tmpf is not None:
            d["theta_approx"] = tmpf + 0.0054 * elev

        derived[s] = d

    return derived


def run_full_analysis(graph: StationGraph,
                      snapshot: dict,
                      rho_base: float = 0.3) -> dict:
    """Run ABRCE operators on all analysis variables for one snapshot."""
    field = snapshot.get("stations", {})
    results = {"timestamp": snapshot.get("timestamp", ""), "variables": {}}

    for var, desc in ANALYSIS_VARS.items():
        a = operator_a(graph, field, var)
        rho = compute_rho(graph, a, rho_base)
        b = operator_b(graph, a)
        r = operator_r(graph, b, rho)
        e = operator_c(r)

        results["variables"][var] = {
            "description": desc,
            "a_gradients": {f"{k[0]}->{k[1]}": v for k, v in a.items()},
            "rho": rho,
            "b_accumulated": {f"{k[0]}->{k[1]}": v for k, v in b.items()},
            "r_circulation": r,
            "e_evolved": e,
        }

    derived = compute_derived_fields(graph, field)

    for dvar in DERIVED_FIELDS:
        derived_field = {}
        for s in graph.node_ids:
            dvals = derived.get(s, {})
            if dvar in dvals:
                derived_field[s] = {dvar: dvals[dvar]}
            else:
                derived_field[s] = None

        a = operator_a(graph, derived_field, dvar)
        rho = compute_rho(graph, a, rho_base)
        b = operator_b(graph, a)
        r = operator_r(graph, b, rho)
        e = operator_c(r)

        results["variables"][dvar] = {
            "description": DERIVED_FIELDS[dvar],
            "raw_values": {s: derived.get(s, {}).get(dvar) for s in graph.node_ids},
            "a_gradients": {f"{k[0]}->{k[1]}": v for k, v in a.items()},
            "rho": rho,
            "b_accumulated": {f"{k[0]}->{k[1]}": v for k, v in b.items()},
            "r_circulation": r,
            "e_evolved": e,
        }

    return results


# ===================================================================
# TEMPORAL EVOLUTION TRACKING
# ===================================================================

def track_evolution(graph: StationGraph,
                    snapshots: list,
                    rho_base: float = 0.3) -> list:
    """Run ABRCE analysis across all snapshots."""
    evolution = []

    for i, snap in enumerate(snapshots):
        result = run_full_analysis(graph, snap, rho_base)

        if i > 0 and evolution:
            prev = evolution[-1]
            deltas = {}

            for var in result["variables"]:
                curr_e = result["variables"][var].get("e_evolved", {})
                prev_e = prev["variables"].get(var, {}).get("e_evolved", {})

                delta = {}
                for s in graph.node_ids:
                    c = curr_e.get(s)
                    p = prev_e.get(s)
                    if c is not None and p is not None:
                        delta[s] = c - p
                    else:
                        delta[s] = None
                deltas[var] = delta

            result["delta_e"] = deltas

        evolution.append(result)

    return evolution


# ===================================================================
# SIGNAL EXTRACTION
# ===================================================================

def extract_precip_signals(graph: StationGraph,
                           evolution: list) -> dict:
    """
    Extract structurally grounded precipitation signals
    from ABRCE evolution.

    PROJECTION LAYER DISCIPLINE (per Verifier note):
        Everything below this line is POST-OPERATOR interpretation.
        These are observational summaries of what the operators
        produced -- not operator-layer constructs.

        The operators computed relational structure.
        This function READS that structure and summarizes it.
        It does NOT feed back into operators or modify their output.

        "Precipitation signal" is a projection-layer label applied
        to a convergence of operator outputs. The operators do not
        know what precipitation is. They detected relational
        gradients, accumulation patterns, and circulation structure.
        WE interpret those patterns as precipitation-relevant.

        If this interpretation is wrong, the operators are not wrong.
        The projection is wrong.

    Returns: dict with per-station, per-timestep signal strength
    """
    signals = {}

    for t_idx, result in enumerate(evolution):
        ts = result.get("timestamp", f"t{t_idx}")
        signals[ts] = {}

        for s in graph.node_ids:
            sig = {"components": {}, "composite": 0.0}

            td_dep = result["variables"].get("td_depression", {})
            td_raw = td_dep.get("raw_values", {}).get(s)
            if td_raw is not None:
                sig["components"]["moisture"] = max(0, 1.0 - td_raw / 30.0)
            else:
                sig["components"]["moisture"] = None

            td_a = td_dep.get("a_gradients", {})
            convergence = 0.0
            edge_count = 0
            for nb, w in graph.neighbors(s):
                g = td_a.get(f"{nb}->{s}")
                if g is not None:
                    convergence += max(0, g)
                    edge_count += 1
            if edge_count > 0:
                sig["components"]["moisture_convergence"] = convergence / edge_count
            else:
                sig["components"]["moisture_convergence"] = None

            td_circ = td_dep.get("r_circulation", {}).get(s)
            if td_circ is not None:
                sig["components"]["moisture_circulation"] = abs(td_circ)
            else:
                sig["components"]["moisture_circulation"] = None

            if t_idx > 0 and "delta_e" in result:
                p_delta = result["delta_e"].get("alti", {}).get(s)
                if p_delta is not None:
                    sig["components"]["pressure_forcing"] = abs(p_delta)
                else:
                    sig["components"]["pressure_forcing"] = None

            components = [v for v in sig["components"].values()
                         if v is not None and v > 0]
            if components:
                product = 1.0
                for c in components:
                    product *= c
                sig["composite"] = product ** (1.0 / len(components))
            else:
                sig["composite"] = 0.0

            signals[ts][s] = sig

    return signals


# ===================================================================
# CONVENIENCE
# ===================================================================

def build_default_graph():
    """Build the SB County station graph from fetch_metar declarations."""
    from fetch_metar import STATIONS, EDGES
    return StationGraph(STATIONS, EDGES)
