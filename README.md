# snowball-limit-cycle

[![DOI](https://zenodo.org/badge/1046432247.svg)](https://doi.org/10.5281/zenodo.16988254)

This repository contains the code and notebooks used to generate the figures for **“Repeated Snowball-hothouse cycles within the Neoproterozoic Sturtian Glaciation.”** The model couples the global carbon, oxygen, and phosphorus cycles to a 1D climate module to explore Snowball limit cycling. As a matter of courtesy, we request that people using this code please cite Minsky et al. (in review).

## Installation

Ensure you have Python 3.8+ and Jupyter installed. Install required dependencies using:

```bash
pip install -r requirements.txt
```

## Code structure

- `model.py` — Core ICEBOX box model
- `ice_albedo.py` — Climate module: albedo parameterization and CO2-temperature interpolation from PCM-LBL data
- `helpers.py` — Plotting utilities and diagnostic calculations (durations, LIP volume, etc.)

## Notebooks and figure outputs

The numeric prefixes match figure numbers. “generate” notebooks run the model and save intermediate data for slower operations; “show” notebooks read saved data and render figures.

### Main text

- `1_canonical_hardball.ipynb` — Canonical hard Snowball trajectories; used to assemble Fig. 1.
- `2_limit_cycle_example.ipynb` — Example limit-cycle run; produces Fig. 2 and supplementary figures Fig. S4 and Fig. S5
- `3_heatmap_generate.ipynb` — Runs the 2-parameter sweep used for Fig. 3
- `3_heatmap_show.ipynb` — Produces Fig. 3

### Supplementary figures

- `SI_bifurcation.ipynb` — Produces Fig. S3: climate bifurcation
- `SI_Cimb_mins.ipynb` — Produces Fig. S6: minimum carbon imbalance threshold for limit cycling
- `SI_LIP_weathering.ipynb` — Produces Fig. S12: LIP weathering parameter space
- `SI_LIP_volume.ipynb` — Produces Fig. S13: LIP volume constraints for 56 Myr duration
- `SI_duration_overall.ipynb` — Produces Fig. S8: overall limit cycle duration across $C_{imb}$-$\tau$ parameter space
- `SI_duration_snowball.ipynb` — Produces Fig. S9: Snowball duration in $V_C$-$W_{sea}$ parameter space
- `SI_duration_interglacial.ipynb` — Produces Fig. S11: interglacial duration sensitivity to $\tau$
- `SI_sweep_generate.ipynb` — Runs Latin hypercube sweep for Fig. S7
- `SI_sweep_show.ipynb` — Produces Fig. S7: parameter distribution comparisons
- `SI_O2_contours.ipynb` — Produces Fig. S14–S15: oxygenation outcomes across $V_C$ and $V_{red}$
- `SI_sensitivities.ipynb` — Produces Fig. S10 and Fig. S16: sensitivity plots for oxygen and interglacial durations
- `SI_Marinoan.ipynb` — Produces Fig. S17: combined Sturtian–Marinoan scenario example

## Contact

Contact: cminsky@g.harvard.edu
