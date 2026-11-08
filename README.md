# Spatial Catch Allocation MSE

Code, analysis data, figures, and manuscript for:

**Can empirical catch allocation strategies based on single-region stock assessments provide adequate advice for managing spatially structured fish populations?**

The study compares five catch-allocation strategies in a closed-loop management
strategy evaluation conditioned on western Atlantic black sea bass
(*Centropristis striata*). Sixteen operating-model scenarios vary population
structure, movement, and recruitment. Results summarize 100 paired retained
realizations per scenario and strategy over a 50-year management period.

## Repository contents

| Directory | Contents |
| --- | --- |
| `manuscript/` | Main manuscript, editable main tables, and supplementary materials |
| `figures/main/` | Figures 1-7 as PNGs; numerical figures also have vector PDFs |
| `figures/supplementary/` | Figures S1-S14 as PNGs and vector PDFs |
| `R/` | Figure reproduction, metric extraction, and reference-point validation |
| `data/` | Compact performance data, reference-point cache, and example trajectory |
| `results/` | Raw summaries, scaled metrics, composite scores, and rankings |
| `simulation/` | Original scenario scripts, seeds, shared conditioning data, and legacy aggregation source |
| `src/reference_points/` | SPR evaluation and root-solving source |
| `environment/` | R session information for the tested figure pipeline |

## Reproduce the figures

The figure pipeline requires R and the packages `dplyr`, `tidyr`, `ggplot2`,
`scales`, and `patchwork`. Package versions are recorded in
[`environment/session-info.txt`](environment/session-info.txt).

From the repository root:

```sh
Rscript --vanilla R/reproduce_figures.R
Rscript --vanilla R/validate_static_risks.R
```

An optional destination keeps regenerated files in a separate directory:

```sh
Rscript --vanilla R/reproduce_figures.R output
```

The script regenerates 19 numerical figures from the included data. Figures 1-2
are archived illustrations and are copied to the destination when a separate
output directory is used. All 21 figures have PNG files. The figure pipeline
uses the included analysis data.

## Performance metrics and scoring

Catch and spawning stock biomass (SSB) are summarized within each realization
over the first and final ten management years. Catch/SSB variability and risk
probabilities use the full management period. Raw-metric figures show medians
across realizations; Figures S6-S11 also show interquartile ranges.

For Figures 5-7 and S12-S14, **raw realization-level metrics are averaged across
the 100 retained realizations before min-max scaling**. Scaling is performed
separately for each operating-model scenario, spatial scale, and metric across
the five allocation strategies. Higher catch and SSB receive higher scores;
lower variability and risk receive higher scores. The composite is the
arithmetic mean of the eight metric scores. Relative-to-Equal catch/SSB fields
are retained in the input files for provenance and are not used for scoring.

Exactly tied catch, SSB, or variability metrics receive a score of 1. For an
exactly tied risk probability `p`, all strategies receive `1 - p`. Scores are
descriptive comparisons within scenarios; raw differences determine the
magnitude and practical importance of strategy effects.

Risk metrics use terminal-year static operating-model reference points, held
fixed over 2024-2073. Overfishing is `F/F40% > 1`; depletion is
`SSB/SSB40% < 0.5`. Details of the data and reference-point cache are in
[`data/README.md`](data/README.md).

## Simulation provenance

[`simulation/`](simulation/) contains the study scenario scripts, seeds, shared
conditioning inputs, and aggregation source. The directory layout and execution
setup are described in [`simulation/README.md`](simulation/README.md).

The MSE framework is [SPASAM-MSE](https://lichengxue.github.io/SPASAM.MSE/),
which uses WHAM. Historical scripts retain the earlier `whamMSE` package name.

## Citation and source licenses

Please cite the associated manuscript when using these materials. Full author
and reference information is provided in the manuscript. Third-party source
licenses are described in [`NOTICE.md`](NOTICE.md); this repository does not
assign a new blanket license to the manuscript, datasets, or artwork.
