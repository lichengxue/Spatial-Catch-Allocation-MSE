# Simulation scripts

| Path | Contents |
| --- | --- |
| `scenarios/` | Sixteen operating-model scenario scripts |
| `seeds/` | Scenario-specific random-number seeds |
| `inputs/` | Shared conditioning data |

Shared conditioning files are stored once. To reproduce the original working
directory layout, place the corresponding inputs and seed file in the scenario
working directory. Execution uses the original R package environment and
HPC/path settings specified in the scripts. The manuscript evaluates Mod1-Mod5.

`Scenario_NUL.R` is the Windows-compatible filename for the NUL scenario.
The figure-reproduction commands are provided in the repository README.
