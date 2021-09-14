r"""Physical constants
This module constains physical constants used throughout the library such
as the acceleration due to gravity, the universal gas constant, the density
of ice and water, etc.
"""


#: number of seconds in a year, for unit conversions
year = 365.2425 * 24 * 60 * 60

#: acceleration due to gravity (m / s^2)
gravity = 9.81

#: heat capacity of air (J kg^-1 K^-1)
heat_capacity_air = 1005.0

#: heat capacity of water (J kg^-1 K^-1)
heat_capacity_water = 4190.0

#: longwave emissivity of ice
emissivity_ice = 1.0

#: longwave emissivity of ocean
emissivity_ocean = .97

#: thermal conductivity of ice
thermal_conductivity_ice = 2.04

#: thermal conductivity of snow
thermal_conductivity_snow = 0.31

#: latent heat of fusion
latent_heat_fusion = 3.340e5

#: latent heat of sublimation
latent_heat_sublimation = 2.834e6

#: latent heat of vaporisation
latent_heat_vapor = 2.501e6

#: stefan-boltzmann constant
stefan_boltzmann_constant = 5.67e-8

#: mixed layer depth density threshold (kg m^-3)
mixed_layer_depth_threshold_density = .03

#: mixed layer depth temperature threshold (K)
mixed_layer_depth_threshold_temperature = 0.2

#: gradient richardson number
gradient_richardson_criteria = .25

#: static bulk richardson number
bulk_richardson_criteria = .65

#: albedo of ocean
albedo_ocean = 0.06

#: albedo of ice
albedo_ice = .5

#: albedo of snow
albedo_snow = .80

#: radiation coefficient 
beta1 = 0.6

#: radiation coefficient
beta2 = 20.0

#: rkz
rkz = 1e-6

#: sea ice density (kg m^-3)
density_ice = 930.0

#: air density (kg m^-3)
density_air = 1.275

#: turbulent transfer coefficient over sea ice
transfer_coefficient_ice = 0.0013

#: turbulent transfer coefficient over ocean
transfer_coefficient_ocean = 0.001

#: Stanton number for mixed layer to sea ice heat transfer
stanton_number_ice = 0.006

#:reference salinity
salinity_reference = 5.0

#: sea ice minimum thickness
sea_ice_minimum_thickness = 0.1



