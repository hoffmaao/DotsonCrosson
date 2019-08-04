# DotsonCrosson
Model experiments and observation of seasonal grounding-line changes in response to seasonal thermocline changes.

# Model Workflow

## Diagnostic runs
1) aggregate data for inversions of ice fluidity parameter and basal resistance. These include Mathieu's latest bed-machine product, which matches gravity inversions for coastal bathymetry. Surface elevation data from REMA and Cryosat-2. Surface velocities from Ian Joughin. Depth variable meltrates derived from power lay temperature relation and glider data. 

This problem is interestingly dependent on groundinglines. We should discuss processing grounding lines for C-D with Ian.

2) Invert for basal resistance using both an effective pressure dependent friction law and a weertman sliding law assuming advected WAIS divide temperatures at the boundary (allow temperature to evolve to steady state).

3) Solve for distributed fluidity parameter assuming new basal resistance field then reinvert for final basal resistance assuming "more consistent" fluidity parameter.

(all of these steps are done with quarterly velocity product when melt rates are their slowest and groundingline is at it's furthest extent -likely the first quarter of the first year of available velocities).

## Prognostic runs
1) Use time series of meltrates from 2018-2019 observations to force melt in the cavity and model consequential changes in grounding line position.

















