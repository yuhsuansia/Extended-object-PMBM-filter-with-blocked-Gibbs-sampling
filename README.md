# Extended-Object PMBM Filtering with Gibbs Sampling

MATLAB implementation of a Poisson multi-Bernoulli mixture (PMBM) filter for
tracking multiple extended objects. Each object's kinematic state, measurement
rate, and extent are represented by a gamma-Gaussian-inverse-Wishart (GGIW)
density. Object measurements follow a Poisson point process, and clutter is
uniform Poisson clutter. There is no separate detection-probability parameter.

This standalone folder contains standard blocked Gibbs and collapsed Gibbs
association sampling, with either simple or DBSCAN initialization. It runs one
measurement sequence at a time. Monte Carlo experiment drivers, saved experiment
results, cross-method plotting scripts, and PMB approximations are excluded.
The random sampling inside the Gibbs filter and scenario generator is retained.

## Related Publication

This repository accompanies the paper:

**[Efficient Implementations of Extended Object PMBM Filters with Blocked Gibbs Sampling](https://arxiv.org/abs/2604.24010)**  
Yuxuan Xia, &Aacute;ngel F. Garc&iacute;a-Fern&aacute;ndez, and Lennart Svensson.

The paper has been accepted for publication in *IEEE Transactions on Aerospace
and Electronic Systems*. A preprint is available on arXiv as
[arXiv:2604.24010](https://arxiv.org/abs/2604.24010).

## Requirements

- MATLAB, tested with R2026a.
- Statistics and Machine Learning Toolbox for `poissrnd`, `gamrnd`, `iwishrnd`,
  and `dbscan`.

No MEX compilation, downloaded dataset, or other folder from the original
project is required. Start with a clean MATLAB session so that functions with
the same names in other filter implementations do not shadow this code.

## Run the Demonstration

Set MATLAB's Current Folder to this folder, then run:

```matlab
main
```

Or, from a terminal in this folder, with MATLAB on your executable search path:

```sh
matlab -batch "main"
```

`main.m` sets up the paths, generates the truth and measurements, plots the
ground truth, runs the filter, and prints each time step followed by timing and
GOSPA/PGOSPA summaries. Filtering time excludes metric evaluation and plotting.
Nothing is saved automatically.

The MATLAB workspace contains:

- `cfg` and `par`: scenario and filter settings.
- `data`: trajectories, extents, object rates, and measurements.
- `results.est{k}`: estimates at time step `k` from the highest-weight global
  hypothesis, with existence greater than `par.existThreshold`.
- `results.gospa`, `results.gospa_decomp`: GOSPA and its localization, missed,
  and false components.
- `results.pgospa`, `results.pgospa_decomp`: the weighted average of the
  per-global-hypothesis MB PGOSPA values and their components.
- `results.tStep`, `results.totalFilterTime`: per-step and total filter runtime.
- `results.state`: the final filter state, predicted to the next time step.

## Configuration

Edit `scenario/defaultConfig.m` for the main experiment settings:

| Setting | Default | Meaning |
| --- | --- | --- |
| `numGibbsIter` | `100` | Stored hypotheses before deduplication: one initializer and 99 Gibbs sweeps |
| `collapsedGibbs` | `true` | Marginalize Bernoulli existence in the sampler |
| `enableClustering` | `false` | Apply DBSCAN to measurements unassigned by the greedy initializer |
| `meanMeasPerObject` | `10` | Mean of the generated object-rate prior and the birth-rate prior |
| `meanClutter` | `10` | Expected clutter measurements per scan; must be positive |
| `truthSeed` | `0` | Seed for true trajectories, extents, and object rates |
| `measurementSeed` | `0` | Seed for measurement generation |
| `numSteps` | `100` | Sequence length |
| `plotGroundTruth` | `true` | Show the ground-truth plot |

The default ten objects move approximately toward the center of a 300 m by
300 m surveillance region, with staggered births and deaths. Object rates and
extents are sampled once per object and remain constant along each trajectory.
When changing `numObjects`, update `appearance` to have one column per object.

The four supported configurations correspond to:

| Variant | `collapsedGibbs` | `enableClustering` |
| --- | --- | --- |
| PMBM-Gibbs1 | `false` | `false` |
| PMBM-Gibbs2 | `false` | `true` |
| PMBM-C-Gibbs1 | `true` | `false` |
| PMBM-C-Gibbs2 | `true` | `true` |

`scenario/buildParameters.m` contains the motion, GGIW prediction, birth,
pruning, and DBSCAN parameters. The inherited DBSCAN defaults are
`dbscanEps = sqrt(trace(extentNominal)/2)` and `dbscanMinPts = 2`. Noise points
become singleton new-object candidates. The simple initializer uses singleton
candidates without clustering. Distance-partition initialization is excluded.

There is no separate maximum-global-hypothesis setting: `numGibbsIter` bounds
the number of candidate global hypotheses per scan; deduplication and
`globalHypoPruneThreshold` reduce this number further. GOSPA uses `p = 1`,
`c = 20`, and `alpha = 2`, set in `runFilter.m`.

For example, run collapsed Gibbs with DBSCAN and 500 stored hypotheses:

```matlab
setupPath;
cfg = defaultConfig();
cfg.collapsedGibbs = true;
cfg.enableClustering = true;
cfg.numGibbsIter = 500;
par = buildParameters(cfg);
data = generateScenario(cfg, par, cfg.measurementSeed);
results = runFilter(data, par);
```

## Reproducibility and Existing Measurements

Scenario generation preserves the source implementation's equations and order
of random draws. It resets the truth stream using `truthSeed`, then resets the
measurement stream using `measurementSeed`. Defaults are both zero, as in the
original demonstration. Within the same MATLAB environment, equal seeds and
scenario parameters reproduce the same truth and measurements regardless of
the Gibbs mode, iteration count, or initialization choice.

`runFilter` does not regenerate measurements or reset the random stream. For a
repeatable filter run on existing `data`, set the RNG immediately before it:

```matlab
rng(42, 'twister');
results = runFilter(data, par, false);  % false disables progress printing
```

Each `data.measCell{k}` must be a `2`-by-`M_k` measurement matrix. Use
`zeros(2,0)` for a scan without measurements. Evaluation also requires
`data.truthTracks` (4-by-T-by-N) and `data.truthExtents` (2-by-2-by-T-by-N),
with `NaN` entries when an object is absent. The state order is
`[x; y; vx; vy]`; extent matrices are 2-by-2 covariances.

For measurements without ground truth, call the filtering API directly:

```matlab
setupPath;
cfg = defaultConfig();
par = buildParameters(cfg);
cache = buildCache(par);
birth = buildPoissonBirth(par);
state = initFilterState(birth);
% Z is the 2-by-M measurement matrix for the current scan.
[state, estimates, posterior] = pmbmStep(state, birth, Z, par, cache);
```

Repeat the final call for successive scans. `posterior` is the updated MBM
before prediction; `state` is ready for the next scan. Match the model and
clutter parameters to the supplied measurements before building the cache.

## Tests

From this folder:

```matlab
setupPath;
testResults = runtests('tests');
assertSuccess(testResults);
```

The tests cover reproducibility, both Gibbs modes and initialization choices,
empty and singleton scans, probability normalization, GGIW updates, and
categorical sampling edge cases. They do not run Monte Carlo performance
experiments.

## Layout

- `main.m`: runnable demonstration.
- `runFilter.m`: filtering and evaluation for one supplied sequence.
- `setupPath.m`: runtime path setup.
- `core/`: PMBM lifecycle and hypothesis construction, caching, and pruning.
- `association/`: initializers and standard/collapsed Gibbs association steps.
- `ggiw/`: prediction, update, likelihoods, moment matching, and state sampling.
- `scenario/`: configuration and synthetic truth/measurement generation.
- `metrics/`: GOSPA, PGOSPA, and assignment utilities.
- `math/`: numerical and categorical-sampling helpers.
- `scripts/`: ground-truth plotting.
- `tests/`: deterministic regression and smoke tests.

Shared initial and sampled hypothesis construction, state drawing, Gaussian
scoring, and categorical sampling are consolidated. The singleton PPP GGIW
update retains its original numerical regularization, and the standard and
collapsed association kernels remain separate because their probabilities differ.

See `THIRD_PARTY_NOTICES.md` for the bundled assignment solver's attribution.
