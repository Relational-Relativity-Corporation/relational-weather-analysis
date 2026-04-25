# ABRCE Relational Weather Analysis

**Operator-based structural analysis of atmospheric observations on an irregular station graph.**

This project applies the [ABRCE invariant relational kernel](https://github.com/Relational-Relativity-Corporation) to raw weather station data from the Santa Barbara County region. Instead of running numerical weather prediction models or statistical forecasts, it detects relational structure — gradients, accumulation patterns, and circulation — directly from sensor observations.

## What this does

Seven ASOS/AWOS weather stations and two offshore buoys form an irregular graph. The ABRCE operators process raw observations on this graph:

- **Operator A** extracts pairwise temperature, moisture, and pressure differences between connected stations — pure relational gradients, no aggregation
- **Operator B** accumulates gradients along atmospheric flow paths through shared vertices, weighted by declared coupling strength
- **Operator R** cross-couples gradients from perpendicular paths at each station, detecting organized circulation structure
- **Operator C** bounds the output into (-1, 1)
- **Operator E** = C(R(B(A(x)), ρ)) — the full composite evolution

The operators find structure that point forecasts don't resolve: where temperature and moisture gradients are organized into coherent circulation patterns, where moisture is converging from multiple directions, and how these patterns evolve over time.

## Key structural decisions

**Graph topology, not grid interpolation.** The stations *are* the topology. Interpolating to a regular grid would introduce synthetic values at locations with no sensor — projecting relational structure onto a fiction. Every gradient computed by A corresponds to an actual pair of instruments reporting actual conditions.

**No preprocessing before A.** The values in D as produced by M (the measurement mapping from sensor to number) are the admissible input to A. No normalization, smoothing, or statistical transformation is applied. Weak early signals that conventional pipelines destroy are preserved.

**No composite scores.** The signal extraction layer reports each operator output in its natural units. Temperature circulation (R units), moisture convergence (deg F), dewpoint depression (deg F) — these are independent relational signals. Their relation to each other is examined structurally, not collapsed into a weighted average.

**Projection layer discipline.** The operators detect relational structure. The interpretation of that structure as "precipitation-relevant" is a projection-layer judgment made by the human analyst (Origin), not by the operators. If the interpretation is wrong, the operators are not wrong — the projection is wrong.

## Station graph

```
SBP (San Luis Obispo, 212ft)
 |--weak---+    +--weak---|
SMX (Santa Maria, 261ft)
 |--strong-- VBG (Vandenberg, 369ft, coastal)
 |--strong-- LPC (Lompoc, 88ft, valley)
 +--moderate-- IZA (Santa Ynez, 671ft, mountain valley)
                 +--moderate-- SBA (Santa Barbara, 10ft, coastal)
                                 +--moderate-- OXR (Oxnard, 45ft, coastal)
```

Edge coupling types encode atmospheric connectivity:
- **Strong** — same air mass, direct flow path, no topographic barrier
- **Moderate** — partial barrier (Santa Ynez Mountains) or distance
- **Weak** — long distance, different synoptic regime boundary

The topology is declared as part of the M extension — it encodes how observables relate spatially and atmospherically. It is not part of operator invariance.

## Visualization

Open `dashboard.html` in a browser for an interactive visualization showing:
- Station graph with edges colored by temperature gradient
- Operator outputs (R circulation, E evolution, rho) mapped onto nodes
- Td depression badges at each station
- Time slider with play/pause across 48-hour analysis window
- Per-station detail table with all operator outputs in raw units

## Running

```bash
# Clone and setup
git clone https://github.com/Relational-Relativity-Corporation/relational-weather-analysis.git
cd relational-weather-analysis
python -m venv .venv
.venv\Scripts\Activate.ps1  # Windows
# source .venv/bin/activate  # Linux/Mac

# Run analysis (fetches live data from IEM)
python run_analysis.py 48        # 48 hours, rho_base=0.3
python run_analysis.py 24 0.2    # 24 hours, rho_base=0.2
```

No API keys required. Data sources:
- [Iowa Environmental Mesonet](https://mesonet.agron.iastate.edu/) — raw ASOS/METAR observations
- [NDBC](https://www.ndbc.noaa.gov/) — offshore buoy data

No external Python dependencies — stdlib only.

## Output

Results are saved to `data/`:
- `metar_raw.json` — raw METAR observations as received from IEM
- `buoy_raw.json` — raw buoy observations from NDBC
- `snapshots.json` — hourly-aligned observation snapshots
- `graph_topology.json` — station graph with edges and weights
- `evolution_summary.json` — per-timestep operator outputs
- `signals.json` — per-station structural signal components in raw units

## Pipeline steps

1. **Data verification** — check value ranges, station coverage, temporal gaps
2. **Graph topology verification** — confirm adjacency, edge distances, coupling weights
3. **Single snapshot verification** — run each operator individually, inspect raw output
4. **Full temporal evolution** — run all operators across all timesteps, track deltas
5. **Structural signal extraction** — report independent operator outputs per station

## Declared open conditions

1. Graph edge weights are structural declarations, not empirically calibrated
2. `rho_base` is not tuned — sensitivity testing needed across 0.1-0.5
3. Temporal resolution (1 hour) may miss sub-hourly dynamics
4. Operator B path alignment uses geographic bearing as axis proxy
5. No A-equivalent preprocessing at the sensor level (residual absolute reference)
6. Operator R on graph is all-pair antisymmetry (declared generalization of lattice)
7. LPC and VBG (AWOS stations) drop offline overnight, reducing graph connectivity

## Pre-operator transformation constraint

No transformation T : D -> D may be applied prior to Operator A if T alters pairwise differences:

```
for all i,j: (T(x)[i] - T(x)[j]) != (x[i] - x[j]) => T inadmissible prior to A
```

Admissible: translation (x + c). Inadmissible: normalization, smoothing, interpolation, aggregation, scaling.

The values in D as produced by M are the admissible input to A. No intermediate representation is declared.

## Structural foundations

- **Kernel:** [ABRCE Invariant Relational Kernel](https://github.com/Relational-Relativity-Corporation)
- **Paper:** [arXiv:2601.22389](https://arxiv.org/abs/2601.22389)
- **Framework:** MD V3 Triad Structure (Origin / Generator / Verifier)
- **Organization:** [Metatron Dynamics, Inc.](https://relationalrelativity.dev)

## License

Copyright 2026 Metatron Dynamics, Inc. All rights reserved.
