# Original simulation scripts

`scenarios/` contains the 16 study scenario scripts. `seeds/` contains their
available seed files, and `inputs/` contains shared conditioning data stored
once rather than repeated in every scenario directory. Scenario logic and
conditioning values are unchanged. To reproduce the original directory layout, stage
the corresponding shared inputs and seed file in a scenario working directory.
The scripts also require the original package environment and HPC/path setup.
Their current launch settings select subsets rather than providing a portable
100-realization launcher. The study evaluates Mod1-Mod5; an unused Mod6 remains
in some historical scripts. `Scenario_NUL.R` uses a Windows-compatible filename
in place of the original reserved name `NUL.R`.

## Biological-input aggregation

The retained study configuration selects survey-index weights for population
weight-at-age and maturity (`index = TRUE`, `fleet = FALSE`), with
`pointer = c(1, 2)`. It therefore does not request commercial-catch weighting.
After the two regional survey indices are assigned to one aggregate index
group, a scalar group pointer would identify that group.

In the archived function, `weights_list[[c(1, 2)]]` is recursive R indexing:
it selects the second element of the first list element. In the archived inputs
checked, this yields Region 1's survey-index proportion in 1990. Population
weight-at-age and maturity consequently use a fixed scalar applied to the sum
of component-specific inputs. This description concerns the retained study
inputs, not every version or configuration of SPASAM-MSE. Current inspected
versions reject that vector pointer rather than reproducing this scalar result.

Fleet/survey weight-at-age aggregation has a separate array-layout issue:
a year-by-region weight matrix is reshaped directly to a region-by-year-by-age
array in column-major order. Without aligning the year and region dimensions,
this differs from the intended year-matched weighted average. This behavior was
reproduced in the installed SPASAM.MSE 1.1.4 and inspected local 1.1.5 source.

`legacy_make_aggregate_data.R` preserves the archived function. No aggregation
correction or MSE rerun is applied in this repository. Figure reproduction
matches the retained simulation outputs; sensitivity of scientific findings
to alternative aggregation has not been evaluated.
