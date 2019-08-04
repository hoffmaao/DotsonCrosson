# DotsonCrosson
This repository includes model experiments and observations of seasonal grounding-line response to thermocline changes.

# Background
Tropical sea surface temperatures have been linked to near surface wind stresses just off the shelf break in front of the Amundsen Sea. These tropically teleconnected Westerly wind patterns affect the thermocline depth beneath the ice shelf and the amount of circumpolar deep water that can melt the shelf cavity walls.


# Hypothesis
The seasonal ice velocity variations we observe across the Crosson and Dotson shelves is caused by seasonal movement of the grounding line due to seasonal melt.

# Model Workflow


## Diagnostic runs
1) First, we aggregate data for inversions of the ice fluidity parameter (rate factor) and basal resistance. These include Mathieu's latest bed-machine product, which matches gravity inversions for coastal bathymetry, surface elevation data from REMA and Cryosat-2, surface velocities from Ian Joughin, and depth variable meltrates derived from a power lay temperature relation and glider data. If we are able to make our own inversion of the bed topography in front of Crosson and Dotson all the better; however, these beds should be consistent with the cavity bathymetry we've mapped out with gliders (this is potentially a whole separate project).

Our hypothesis predicts that the smith, pope and kolher glacier grounding lines are modulated by seasonal melt rates and that there should be a lag between the grounding line response and ocean forcing that we cannot easily detect in quarterly averaged velocity products. To better constrain these two signals we will initialize teh

2) Invert for basal resistance using both an effective pressure dependent friction law and a weertman sliding law assuming advected WAIS divide temperatures at the boundary (allow temperature to evolve to steady state).

3) Solve for distributed fluidity parameter assuming new basal resistance field then reinvert for final basal resistance assuming the "more consistent" fluidity parameter.

(all of these steps are done with quarterly velocity products when melt rates are their slowest and the groundingline is at it's furthest extent -likely the first quarter of the first year of available velocities-January of 2015).

## Prognostic runs

Use time series of sub-shelf temperatures from 2018-2019 observations to force melt in the cavity and model any modulating affect these meltrates have on grounding line position.

We can then compare high frequency grounding-line response predicted by simple models (i.e. linearized marine outlet glacier model with prograde bed; Robel et al. 2018) with the power spectra predicted by our more sophisticated model of outlet glacier response.

















