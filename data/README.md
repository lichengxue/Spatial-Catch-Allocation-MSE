# Analysis data

`performance/` contains one `performance2_*.RDS` file for each of 16 scenarios.
Each file has 100 retained realizations, five strategies, three spatial scales,
and three period rows (`short`, `long`, `all`), giving 72,000 rows in total.

Core fields are `Realization`, `Model`, `Period`, `Region`, `Catch`, `SSB`,
`AAV_Catch`, `AAV_SSB`, `p_overfishing`, and `p_overfished`. Catch/SSB values
are temporal means within a realization. AAV is the mean absolute annual change
divided by the mean over the management period. The two risk fields are annual
exceedance proportions over 2024-2073. Legacy `rel_Catch` and `rel_SSB` fields
are not used in the published scoring pipeline.

`reference_points/` contains the matched static-reference cache and compact
realized F/SSB trajectories. Static biological/fishery inputs are averaged over
2069-2073; true mean recruitment uses the full simulated series. Regional and
system realized F columns are paired with the corresponding reference columns.
SPR evaluation and bracketed root-solving source is provided in
`src/reference_points/`.
The included validation script checks all 48,000 risk probabilities against
these cached trajectories and static reference points.

`trajectories/NBL_Example_Trajectory.csv` supplies Figure S1.

## Scenario labels

| Data-file code | Manuscript label |
| --- | --- |
| NB | NBL |
| NBH | NBH |
| NU | NUL |
| NUH | NUH |
| NBI | NBLI |
| NBO | NBLO |
| NBC | NBLD |
| NUC | NULD |
| NBHS | NBH-BH |
| MB | MBL |
| MBH | MBH |
| MU | MUL |
| MUH | MUH |
| MBI | MBLI |
| MBO | MBLO |
| MBHS | MBH-BH |
