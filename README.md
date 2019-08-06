# DotsonCrosson
This repository includes model experiments and observations of seasonal grounding-line response to thermocline changes.

# Background
The history of advances in satellite and moored ocean observations has made it clear that ice-sheet ocean interactions take place across a broad range of spatial and temporal scales and regulate global overturning circulation and global sea-level rise. Remote-sensing observations from the Amundsen Sea Embayment suggest that marine ice-sheet instability may well be underway for West Antarctica; however, the pacing of WAIS collapse and ensuing sea-level contribution in the coming decades remains highly uncertain, in part, due to variance in predicted glacier response to ocean heat content variability. Limited sub-annual observations have restricted studies to interannual and decadal syntheses of climate reanalysis, float observations, and satellite derived grounding-line observations. Here we present the first continuous observations of sub-annual melt rates in the ASE derived from a Seaglider-float occupation in the cavity of the Dotson and Crosson ice shelf from 2018-2019. The data gathered autonomously during the year-long campaign suggest (i) the observed amplitude of seasonal ocean heat content in the year 2018 is similar to the amplitude of multi-year temperature variability over the previous 10 years; (ii) mesoscale intrusions of open ocean water into the ice-shelf cavity drives large spatial and temporal temperature variability across the ice-shelf calving front; and (iii) near the ice-ocean interface, water properties (i.e. salinity, temperature, eddy covariant structure) and stratification vary spatially and evolve sub-annually. Using a higher order model for ice deformation and sliding, we simulate the influence of seasonally and spatially variable submarine melt on grounding line migration. Our model experiments and observations link seasonal variations in 10m zonal winds along the continental-shelf slope to seasonal changes in the cavity thermocline and sub-annual variations in Crosson and Dotson ice-shelf velocity.


# Hypothesis
The seasonal ice velocity variations we observe across the Crosson and Dotson shelves are caused by movement of the grounding line due to seasonal melt.

# Model Workflow


## Diagnostic runs
1) First, we aggregate data for inversions of the ice fluidity parameter (rate factor) and basal resistance. These include Mathieu's latest bed-machine product, which matches gravity inversions for coastal bathymetry, surface elevation data from REMA and Cryosat-2, surface velocities from Ian Joughin, and depth variable ocean meltrates derived from a power lay temperature relation and glider data. If we are able to make our own inversion of the bed topography in front of Crosson and Dotson all the better; however, these beds should be consistent with the cavity bathymetry we've mapped out with gliders (this is potentially a whole separate project).

Our hypothesis predicts that the smith, pope and kolher glacier grounding lines are modulated by seasonal melt rates and that there should be a lag between the grounding line response and ocean forcing that we cannot easily detect in quarterly averaged velocity products. To better constrain these two signals we will:

2) Invert for basal resistance using both an effective pressure dependent friction law and a weertman sliding law assuming advected WAIS divide temperatures at the boundary (allow temperature to evolve to steady state).

3) Solve for distributed fluidity parameter assuming new basal resistance field then reinvert for final basal resistance assuming the "more consistent" fluidity parameter.

(these steps are done with quarterly velocity products when melt rates are their slowest and the groundingline is at it's furthest extent -likely the first quarter of the first year of available velocities-January of 2015).

## Prognostic runs

Use time series of sub-shelf temperatures from 2018-2019 observations to force melt in the cavity and model any modulating affect these meltrates have on grounding line position.

We can then compare high frequency grounding-line response predicted by simple models (i.e. linearized marine outlet glacier model with prograde bed; Robel et al. 2018) with the power spectra predicted by our more sophisticated model of outlet glacier response.

















