# M31 Bulge Microlensing Detection via Difference Imaging Analysis (DIA)

## Project Overview
This repository contains the core Python data analysis and optimization scripts developed as part of a Third Year Physics Laboratory experiment at the University of Manchester. The research and experiment were conducted in collaboration with Ashwin Knight. 

Please note that the code provided in this repository represents specific computational modules (hyperparameter optimization and automated transient detection) which form the backbone of the larger experimental pipeline.

## Experiment Summary (TL;DR)
For those who wish to understand the core findings without reading the full report:

* Goal: To identify gravitational microlensing candidates in the M31 (Andromeda) galactic bulge and estimate the masses of the lens objects (such as dark matter candidates or compact objects) using observational data from the Angstrom Project.
* Methodology: 
    * Difference Imaging Analysis (DIA): Utilized the legacy 'Isis' package to subtract static background light from images, isolating true flux variations caused by astrophysical events.
    * Hybrid Optimization: Automated the highly sensitive DIA kernel and background parameter tuning process. This was achieved using a two-stage approach: a global stochastic search (SciPy's Dual Annealing) to avoid local minima, followed by precise local refinement (Powell's method).
    * Anomaly Detection: Extracted light curves and ranked them using a smoothed geometric mean of z-scores to systematically find statistically significant transient events.

For the theoretical physics background, mathematical modelling of the Point-Source Point-Lens (PSPL) events, and data analysis, please refer to the attached Report and Presentation files.

## Repository Structure
```text
.
|-- hyper_parameter_optimization/
|   |-- optimize_global.py       # Global parameter search using Dual Annealing
|   `-- optimize_local_powell.py # Local parameter refinement using Powell's method
|-- search_microlensing/
|   |-- search_microlensing.py   # Light curve ranking and transient detection script
|   `-- example_images/          # Folder containing sample output plots and light curve data
|-- Hong_MinKi_DIA_Report.pdf    # Full experimental research report
|-- Presentation_dia.pdf         # Summary presentation slides
└-- README.md
```

## Code Module Descriptions

### 1. Hyperparameter Optimization (`hyper_parameter_optimization/`)
Difference Imaging Analysis requires highly precise tuning of Gaussian components to properly match the Point Spread Function (PSF) and background between reference and target images. 
Because the Isis C-shell macros are computationally expensive to run and the parameters are strongly correlated, the parameter tuning was strategically decoupled into two scripts:
* optimize_global.py: Wraps the legacy C-shell scripts and utilizes Generalized Simulated Annealing (Dual Annealing with local search disabled) to explore the broad parameter space and identify promising global minima.
* optimize_local_powell.py: Takes the regions found by the global search and applies Powell's direction-set method to precisely pinpoint the optimal continuous parameters (sigmas).

### 2. Microlensing Search (`search_microlensing/`)
Once the optimal difference images were generated, a variability map was used to extract raw light curves. The `search_microlensing.py` script automates the filtering and identification of true transient events among thousands of data points:
* Calculates the absolute z-score of the flux variations to standardize deviations.
* Applies a rolling geometric mean (window smoothing) to suppress observational noise and outliers.
* Ranks all light curves based on their maximum anomaly score and automatically generates visual plots for the top candidates for further scientific inspection.
