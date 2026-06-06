# Stability Region Mapper for EFTCAMB / Cobaya Models

**Author:** Lorenzo Baldazzi  
**Affiliation:** University of Rome Tor Vergata  

A pair of Python scripts to map and visualise the **stability region** of any modified-gravity or dark-energy model implemented in [EFTCAMB](https://eftcamb.org) and sampled via [Cobaya](https://cobaya.readthedocs.io). The approach is completely general: although the worked example targets the Shift Symmetric Horndeski model, adapting to any other model requires only editing a configuration block.

---

## Table of Contents

- [Overview](#overview)
- [How It Works](#how-it-works)
  - [1. Sobol Quasi-Monte Carlo Sampling](#1-sobol-quasi-monte-carlo-sampling)
  - [2. Stability Evaluation via Cobaya](#2-stability-evaluation-via-cobaya)
  - [3. Adaptive Boundary Refinement](#3-adaptive-boundary-refinement)
  - [4. Visualisation](#4-visualisation)
- [Scripts](#scripts)
  - [`stability_map.py`](#stability_mappy)
  - [`stability_viz.py`](#stability_vizpy)
- [Requirements](#requirements)
- [Usage](#usage)
  - [Generating the Stability Map](#generating-the-stability-map)
  - [Visualising the Results](#visualising-the-results)
- [Adapting to Your Own Model](#adapting-to-your-own-model)
- [Output Format](#output-format)
- [Example: Shift Symmetric Horndeski](#example-shift-symmetric-horndeski)

---

## Overview

When running a Bayesian parameter estimation with EFTCAMB, a large fraction of the prior volume is often physically *unstable* — the model fails ghost or gradient stability conditions and EFTCAMB returns a log-likelihood of $-\infty$. Understanding *where* this happens in parameter space is useful for:

- diagnosing MCMC convergence issues caused by large unstable regions;
- designing better-informed priors;
- producing publication-quality stability maps as figures.

This toolkit answers the question: **given a set of free parameters, which combinations lead to a stable model?**

---

## How It Works

### 1. Sobol Quasi-Monte Carlo Sampling

The free parameters are sampled over a user-defined bounding box using a **scrambled Sobol sequence** (via `scipy.qmc.Sobol`). A Sobol sequence is a low-discrepancy QMC sequence: it fills the parameter space more uniformly than pseudo-random draws, achieving an integration error that scales as $\mathcal{O}((\log N)^d / N)$ rather than the Monte Carlo $\mathcal{O}(N^{-1/2})$. The sample size is internally rounded up to the next power of 2 to exploit Sobol's optimal balance properties.

### 2. Stability Evaluation via Cobaya

For each sampled point $\theta$, a minimal Cobaya model is instantiated (once per worker process) with:

- the **theory block** copied verbatim from the user's YAML (preserving all EFTCAMB flags and `extra_args`);
- the **`one` likelihood** — a no-op that returns $\log\mathcal{L} = 0$ whenever the theory initialises successfully, and $-\infty$ otherwise;
- all nuisance parameters **pinned** at fixed reference values.

A point is labelled **stable** if `model.loglike(θ)` returns a finite value, meaning EFTCAMB completed successfully under the chosen ghost/gradient stability flags. Evaluation is parallelised over a `multiprocessing.Pool` using the `spawn` start method, which is safe with EFTCAMB's Fortran/C shared state.

### 3. Adaptive Boundary Refinement

After the base sample is evaluated, an optional boundary-refinement loop increases point density near the stability boundary, where classification uncertainty is highest.

**Uncertainty score.** For each point $i$, let $f_i$ be the fraction of stable points among its $k$ nearest neighbours (in coordinates normalised to $[0,1]^d$). The uncertainty score is:

$$u_i = 1 - |2f_i - 1|$$

This equals 0 in a homogeneous neighbourhood (all stable or all unstable) and 1 when exactly half the neighbours are of each class — the strongest signal of proximity to the boundary.

**Refinement.** Points in the top decile of $u_i$ are selected. A new Sobol batch is drawn inside the axis-aligned bounding box of those points, padded slightly and clipped to the physical prior. This concentrates new evaluations where they are most informative.

### 4. Visualisation

The labelled point cloud is saved to a pickle file and passed to the visualisation script, which produces two corner-style figures:

- **Scatter corner**: each 2-D panel shows stable (green) and unstable (red) points projected onto a pair of parameters. Diagonal panels show 1-D marginal histograms.
- **Marginalised stability-fraction heatmap**: each 2-D panel shows the estimated probability of stability $\hat{p}(\text{stable}\,|\,\theta_i, \theta_j)$, obtained by binning points and averaging the binary label over all other parameters. Bins with fewer than a minimum count are masked.

---

## Scripts

### `stability_map.py`

Generates and saves the labelled point cloud.

| Section | Responsibility |
|---|---|
| **Configuration** | Parameter ranges, fixed values, YAML path |
| `build_cobaya_info` | Constructs the minimal Cobaya `info` dict |
| `make_model` | Lazy-imports and instantiates the Cobaya model |
| `is_stable` | Single-point stability check |
| `sobol_points` | Sobol QMC sampling in an arbitrary box |
| `evaluate_points` | Parallel evaluation driver |
| `boundary_uncertainty` | K-NN uncertainty score |
| `refine_box` | Draws the next refinement batch |
| `main` | CLI entry point, orchestrates all steps |

### `stability_viz.py`

Reads the pickle produced by `stability_map.py` and writes two PNG figures.

| Function | Output |
|---|---|
| `scatter_corner` | Corner scatter plot (stable/unstable coloured points) |
| `marginal_fraction_corner` | Corner heatmap of marginalised $P(\text{stable})$ |

---

## Requirements

```
cobaya
eftcamb          # compiled and on PYTHONPATH / LD_LIBRARY_PATH
camb
scipy            # >= 1.7 for scipy.qmc.Sobol
scikit-learn     # for NearestNeighbors in the refinement step
numpy
matplotlib
pyyaml
```

---

## Usage

### Generating the Stability Map

```bash
# Quick smoke test (~1 min on 4 cores)
python stability_map.py --base 64 --refine-iters 0

# Full run
python stability_map.py \
    --yaml  /path/to/your_model.yaml \
    --base  4096 \
    --refine-iters 3 \
    --refine-batch 1024 \
    --output stability_map.pkl
```

| Argument | Default | Description |
|---|---|---|
| `--yaml` | (set in script) | Path to the Cobaya YAML |
| `--base` | 4096 | Base Sobol sample size |
| `--refine-iters` | 3 | Number of boundary-refinement passes |
| `--refine-batch` | 1024 | Points added per refinement pass |
| `--workers` | `ncpu - 1` | Parallel worker processes |
| `--seed` | 42 | Sobol scramble seed |
| `--output` | (set in script) | Output pickle path |
| `--serial` | False | Disable parallelism (useful for debugging) |

### Visualising the Results

```bash
python stability_viz.py stability_map.pkl --bins 30 --outdir ./figs
```

| Argument | Default | Description |
|---|---|---|
| `pickle_path` | (positional) | Pickle file from `stability_map.py` |
| `--bins` | 25 | Bins per axis for the heatmap |
| `--min-count` | 3 | Minimum points per bin before masking |
| `--max-scatter` | 20000 | Subsample cap for the scatter plot |
| `--outdir` | `.` | Output directory |

---

## Adapting to Your Own Model

The scripts are model-agnostic. To use them with a different EFTCAMB model, edit the configuration block at the top of `stability_map.py`:

```python
# 1. Point to your Cobaya YAML
DEFAULT_YAML = "/path/to/your_model.yaml"

# 2. Define the parameters you want to scan and their ranges
SS_RANGES = {
    "your_param_1": (lo_1, hi_1),
    "your_param_2": (lo_2, hi_2),
    # add as many as needed
}

# 3. Pin all other parameters at fixed reference values
SS_FIXED = {
    "ombh2": 0.0224,
    "omch2": 0.118,
    # ... LCDM parameters, EFT fixed values, etc.
}
```

No other changes are required. The LATEX display names used in figures can be updated in `stability_viz.py`:

```python
LATEX = {
    "your_param_1": r"$\alpha_1$",
    "your_param_2": r"$m$",
}
```

The sampling, evaluation, refinement, and plotting pipeline then runs unchanged for any number of free parameters $d \geq 1$.

---

## Output Format

`stability_map.py` writes a pickle file containing a single dictionary:

```python
{
    "param_names":  list[str],          # ordered parameter names
    "param_ranges": dict[str, (lo, hi)],# parameter bounding boxes
    "ss_fixed":     dict[str, float],   # pinned parameter values
    "points":       np.ndarray,         # shape (N, d), sampled points
    "stable":       np.ndarray,         # shape (N,), bool stability labels
    "iter_info":    list[dict],         # per-phase statistics
    "yaml":         str,                # absolute path to the YAML used
}
```

This format is self-contained and can be read independently of the scripts that produced it.

---

## Example: Shift Symmetric Horndeski

The default configuration explores the 2-D parameter space of the **Shift Symmetric Horndeski** model as implemented in EFTCAMB, scanning:

| Parameter | Symbol | Range |
|---|---|---|
| `Shift_Symmetric_alphaB0` | $\alpha_{B,0}$ | $[-0.1,\ 3.0]$ |
| `Shift_Symmetric_m` | $m$ | $[-1.5,\ 10.0]$ |

The six $\Lambda$CDM parameters are held fixed at fiducial values, and the remaining model parameters (`alphaK0`, `EFTw0`, `EFTwa`) are declared as constants directly in the YAML. This example is included purely as a concrete illustration; the methodology applies equally to any model whose stability region one wishes to characterise.
