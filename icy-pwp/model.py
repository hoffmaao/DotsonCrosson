import numpy as np
import seawater as sw
from constants import (
    gravity as g,
    beta1 as β1,
    beta2 as β2,
    gradient_richardson_criteria as rg,
    bulk_richardson_criteria as rb,
    mixed_layer_depth_threshold_density  as mld_thresh_density,
    mixed_layer_depth_threshold_temperature as mld_thresh_temperature,
    sea_ice_minimum_thickness as min_thickness,
    rkz as rkz,
    density_air as ρ_a,
    density_ice as ρ_i,
    heat_capacity_air as c_a,
    heat_capacity_water as c_w,
    albedo_snow as α_s,
    albedo_ice as α_i,
    albedo_ocean as α_o,
    transfer_coefficient_ocean as C_o,
    transfer_coefficient_ice as C_i,
    emissivity_ice as ϵ_i,
    emissivity_ocean as ϵ_o,
    thermal_conductivity_ice as k_i,
    thermal_conductivity_snow as k_s,
    latent_heat_sublimation as L_s,
    latent_heat_vapor as L_v,
    latent_heat_fusion as L_f,
    stanton_number_ice as c_h,
    stefan_boltzmann_constant as σ,
    salinity_reference as s_ref,
    year
    )
from scipy import optimize

def icy_pwp(
    sea_ice_thickness, sea_ice_concentration, temperature_relax,
    salinity_relax, temperature, salinity, density,
    uvel, vvel, taux, tauy, z, Q_sw, Q_lw, temperature_surface,
    temperature_ambient, pressure, wind_speed, specific_humidity,precipitation,
    ϕ, melt, dt, T_cdw=1.0, T_out=-1.84, S_cdw=34.51, S_out=34.03, half_depth=450, snow_thickness=0.1,
    A_max=.9, Rb=.75, divergence=1.0, lat = 74.0, domain = 1e4,shelf_length=4e4):
    r"""
    pwp model modified to include sea ice.


    """

    dz = z[1]-z[0]
    half_ind=int(dz*half_depth)

    f = sw.f(lat)
    ucon = (0.1*np.abs(f))
    temperature_freezing = sw.fp(salinity,z)
    density=sw.dens0(salinity,temperature)





    temperature_ocean=optimize.root(energy_balance_ocean,temperature_ambient,args=(temperature_ambient,density[0],temperature[0]+273.15,absorb(β1,β2,len(z),dz),Q_sw,Q_lw,wind_speed,pressure,specific_humidity),tol=10e-20)['x']-273.15
    temperature_ice=optimize.root(energy_balance_ice,temperature_ambient,args=(temperature_ambient,temperature_freezing[0]+273.15,sea_ice_thickness,snow_thickness,Q_sw,Q_lw,wind_speed,pressure,specific_humidity),tol=10e-20)['x']-273.15
    

    #if (temperature_ice>(temperature_freezing[0])) and (sea_ice_thickness>min_thickness):
    #    temperature_surface = temperature_freezing[0]+273.15
    #flux_balance_ocean=energy_balance_ocean(temperature_ocean,temperature_ambient,density[0],temperature[0]+273.15, absorb(β1,β2,len(z),dz),Q_sw,Q_lw,wind_speed,pressure,specific_humidity)
    #flux_balance_ice=energy_balance_ice(temperature_ice,temperature_ambient,temperature[0]+273.15,temperature_freezing[0]+273.15,sea_ice_thickness,snow_thickness,Q_sw,Q_lw,wind_speed,pressure,specific_humidity)


    ### profile of velocity
    profile = ((temperature - temperature_freezing) - np.mean(temperature-temperature_freezing))#/((temperature[0] - temperature_freezing[0]) - np.mean(temperature-temperature_freezing))
    velocity_relax=-ϕ*profile

    #h_interface=len(velocity_relax)-(np.abs(velocity_relax[half_ind:]).argmin()+half_ind)
    h_interface=(np.abs(z-half_depth)).argmin()

    ### relaxation profiles
    T_relax=np.tanh((z-z[h_interface])/200)*(T_cdw-T_out)/2+(T_cdw+T_out)/2
    S_relax=np.tanh((z-z[h_interface])/250)*(S_cdw-S_out)/2+(S_cdw+S_out)/2

    Q_profile = surface_flux(sea_ice_concentration,temperature[0],density[0],temperature_freezing[0],temperature_ocean,Q_sw,wind_speed,z,dz)


    inv_time_length = np.zeros(shape=len(z))
    inv_time_length = np.absolute(velocity_relax)/domain#+max(velocity_relax)/12/domain
    
    dTdt = inv_time_length*(T_relax-temperature)
    dSdt = inv_time_length*(S_relax-salinity)


    ### Absorb radiation at depth and flux melt at the surface
    temperature +=  Q_profile*dt/(dz*density*c_w) + dTdt*dt

    ### compute new density
    density = sw.dens0(salinity, temperature)



    ϕ_A=Φ_A(density[0],sea_ice_thickness,sea_ice_concentration,Rb,temperature_ocean,temperature_freezing[0],wind_speed,A_max)
    ϕ_h=Φ_h(density[0],snow_thickness,sea_ice_thickness,sea_ice_concentration,wind_speed,temperature_freezing[0],temperature_ice,temperature_ocean,temperature[0],dt,dz,Rb)
    ϕ_R=Φ_R(density[0],sea_ice_concentration,temperature_freezing[0],temperature_ocean,wind_speed,A_max)


    sea_ice_concentration += -divergence/year*dt + ϕ_A*dt
    sea_ice_thickness += (ϕ_h+ϕ_R)*dt
    sea_ice_thickness = np.maximum(sea_ice_thickness,min_thickness)
    sea_ice_concentration = np.maximum(np.minimum(sea_ice_concentration,A_max),0.0)

    Q_lat=latent(temperature_ocean+273.15,pressure,specific_humidity,wind_speed)

    salinity[0]+=(salt_flux(density[0],salinity[0],ϕ_h,ϕ_R,ϕ_A,sea_ice_thickness,sea_ice_concentration,melt,Q_lat,precipitation)/dz+dSdt[0])*dt
    salinity[1:] += dSdt[1:]*dt


    ### relieve static instability
    temperature, salinity, density, uvel, vvel = remove_si(temperature, salinity, density, uvel, vvel)
    ### advection
    temperature_freezing = sw.fp(salinity,z)
    ### compute temperature at the grounding zone.
    #temperature_grounding = temperature[ind_grounding]-temperature_freezing[ind_grounding]

    ### compute MLD
    mld_ind= np.flatnonzero(density-density[0]>mld_thresh_density)[0]
    #mld_ind = np.flatnonzero(np.abs(temperature-np.mean(temperature[:int(10/dz)]))>mld_thresh_temperature)[0]
    #mld_ind=min(mld_ind_density,mld_ind_temp)
    mld = z[mld_ind]


    ### rotate u,v, do wind input, rotate again, apply mixing
    ang = -f*dt/2
    uvel, vvel = rot(uvel, vvel, ang)

    if sea_ice_thickness>0.1:
        taux = taux/3.
        tauy = tauy/3.
    du = (taux/((mld+dz)*density[0]))*dt
    dv = (tauy/((mld+dz)*density[0]))*dt
    uvel[:mld_ind] = uvel[:mld_ind]+du
    vvel[:mld_ind] = vvel[:mld_ind]+dv
    ### apply drag
    if ucon > 1e-10:
        uvel = uvel*(1-dt*ucon)
        vvel = vvel*(1-dt*ucon)
    uvel, vvel = rot(uvel, vvel, ang)

	### Apply bulk richardson number instability form of mixing
    if rb > 1e-5:
        temperature, salinity, density, uvel, vvel = bulk_mix(temperature, salinity, density, uvel, vvel, g, rb, len(z), z, mld_ind)

	### Do the gradient Richardson number instability form of mixing
    if rg > 0.0:
        temperature, salinity, density, uvel, vvel = grad_mix(temperature, salinity, density, uvel, vvel, dz, g, rg, rb, len(z), z, mld_ind)

	### Apply diffusion
    diffusivity = dt*rkz/dz**2
    temperature = diffuse(len(z), temperature, diffusivity) 
    salinity = diffuse(len(z), salinity, diffusivity) 
    density = sw.dens0(salinity, temperature)
    uvel = diffuse(len(z), uvel, diffusivity)
    vvel = diffuse(len(z), vvel, diffusivity)



    ### integration of heat fluxes

    HE=ρ_i*L_f*sea_ice_thickness*sea_ice_concentration*divergence/year
    HR=np.sum((T_relax-temperature)*c_w*density*inv_time_length)*dz
    HS=np.sum(Q_profile)*dz

    ### integration of salt fluxes

    SE=sea_ice_thickness*sea_ice_concentration*(salinity[0]-s_ref)*divergence/year
    SR=np.sum((S_relax-salinity)*inv_time_length)*dz
    SP=(precipitation-Q_lat/L_v)*salinity[0]*(1-sea_ice_concentration)/1000
    SM=melt* melt/1000*(salinity[0])



    ### Output
    return temperature, salinity, density, uvel, vvel, mld, sea_ice_thickness, sea_ice_concentration, temperature_ocean, temperature_ice, Q_profile[0], velocity_relax, HE, HR, HS, SE, SR, SP, SM


def absorb(beta1, beta2, zlen, dz):
    
    # Compute solar radiation absorption profile. This
    # subroutine assumes two wavelengths, and a double
    # exponential depth dependence for absorption.
    # 
    # Subscript 1 is for red, non-penetrating light, and
    # 2 is for blue, penetrating light. rs1 is the fraction
    # assumed to be red.
    
    rs1 = 0.6
    rs2 = 1.0-rs1
    z1 = np.arange(0,zlen)*dz
    z2 = z1 + dz
    z1b1 = z1/beta1
    z2b1 = z2/beta1
    z1b2 = z1/beta2
    z2b2 = z2/beta2
    absrb = rs1*(np.exp(-z1b1)-np.exp(-z2b1))+rs2*(np.exp(-z1b2)-np.exp(-z2b2))
    
    return absrb


def pv(T):
    r"""
    calculate the

    """
    return 2.53e8*np.exp(-5420./T)

def qsat(T,P):
    r"""
    calculate the near surface specific humidity
    
    """
    return .622*pv(T)/(P - .378*pv(T))

def sensible(T,T_a,u,C=C_o):
    r"""
    calculate the sensible heat flux
    
    """
    return ρ_a*c_a*C*u*(T-T_a)

def latent(T,P,q,u,C=C_o,L=L_v):
    r"""
    claculate the latent heat flux
    
    """
    return ρ_a*L*C*u*(qsat(T,P)-q)


def ustar(u,C,density):
    r"""
    calculate the friction velocity
    
    NOTE: likely overestimates ustar in the area covered case.
    """

    return u*np.sqrt((ρ_a/density)*C)


def surface_flux(A,temperature,density,temperature_freezing,temperature_ocean,Q_sw,u,z,dz):


    F_mi = density*c_w*c_h*ustar(u,C_i,density)*(temperature-temperature_freezing)
    F_mo = density*c_w*ustar(u,C_o,density)*(temperature-temperature_ocean)


    Q_profile=np.zeros(len(absorb(β1,β2,len(z),dz)))
    Q_profile[1:] = (1-A)*(1-α_o)*Q_sw*absorb(β1,β2,len(z),dz)[1:]
    Q_profile[0]= -(1-A)*F_mo-A*F_mi
    

    return Q_profile


def energy_balance_ocean(temperature_ocean,temperature_ambient,density,temperature,I,Q_sw,Q_lw,u,P,q):
    r"""
    calculate ocean surface energy balance
    
    """

    surface_flux = density*c_w*ustar(u,C_o,density)*(temperature-temperature_ocean)
    Q_sens = sensible(temperature_ocean,temperature_ambient,u,C_o)
    Q_lat = latent(temperature_ocean,P,q,u,C_o,L_v)

    return Q_sens + Q_lat + ϵ_o*σ*temperature_ocean**4 - (1-α_o)*I[0]*Q_sw - ϵ_o*Q_lw - surface_flux


def energy_balance_ice(temperature_ice,temperature_ambient,freezing_temperature,h_i,h_s,Q_sw,Q_lw,u,P,q):
    r"""
    calculate ice surface energy balance
    
    """

    surface_flux = k_i*k_s*(freezing_temperature-temperature_ice)/(k_s*h_i+k_i*h_s)    
    Q_sens = sensible(temperature_ice,temperature_ambient,u,C_i)
    Q_lat = latent(temperature_ice,P,q,u,C_i,L_s)

    return Q_sens + Q_lat + ϵ_i*σ*temperature_ice**4 - (1-α_s)*Q_sw - ϵ_i*Q_lw  - surface_flux



def H_fr(density,temperature_freezing,temperature_ocean,u,C=C_o):
    r"""
    calculate ocean surface heat potential

    """
    return density*c_w*ustar(u,C,density)*(temperature_freezing-temperature_ocean)



def Φ_A(density,h,A,Rb,temperature_ocean,temperature_freezing,u,A_max):
    r"""
    calculate sea-ice concentration change balancing ocean surface
    heat potential with the latent heat released/absorbed
    by the ice growth or melt.

    """

    if (A<A_max) and (temperature_ocean<temperature_freezing):
        dAdt = H_fr(density,temperature_freezing, temperature_ocean,u)*(1-A)/(ρ_i*L_f*h)

    elif (A>0) and (temperature_ocean>temperature_freezing):
        dAdt = H_fr(density,temperature_freezing, temperature_ocean,u)*(1-Rb)*(1-A)/(ρ_i*L_f*h)
    else:
        dAdt=0
    return dAdt


def Φ_R(density,A,temperature_freezing,temperature_ocean,u,A_max,C=C_o):
    r"""
    calculate the sea-ice thickness change balancing the ocean
    surface heat potential after reaching max concentration.

    """
    if A==A_max:
        dhdt = H_fr(density,temperature_freezing,temperature_ocean,u)/(ρ_i*L_f)
    else:
        dhdt = 0
    return dhdt

def Φ_h(density,h_s,h_i,A,u,temperature_freezing,temperature_ice,temperature_ocean,temperature,dt,dz,Rb):
    r"""
    calculate the sea-ice thickness change due to differences
    in heat fluxes at the ice-ocean interface

    """

    if (A>0) and (temperature_ocean>temperature_freezing):
        F_sb = -H_fr(density,temperature_freezing,temperature_ocean,u)*Rb*(1-A)
    else:
        F_sb = 0

    F_mi = density*c_w*c_h*ustar(u,C_i,density)*(temperature-temperature_freezing)

    F_c = k_i*k_s*(temperature_freezing-temperature_ice)/(k_s*h_i+k_i*h_s)  

    return (F_c - F_mi - F_sb)/(ρ_i*L_f)

def salt_flux(density,salinity,ϕ_h,ϕ_R,ϕ_A,sea_ice_thickness,A,melt,Q_lat,precipitation):
    r"""
    calculate surface salt flux

    """
    evaporation = Q_lat/L_v

    return -ρ_i/density*(s_ref-salinity)*((ϕ_h+ϕ_R)*A+ϕ_A*sea_ice_thickness) -(precipitation-evaporation)*salinity*(1-A)/1000 - melt/1000*(salinity)

def remove_si(t, s, d, u, v):
    r"""
    Find and relieve static instability that may occur in the
    density array 'd'. This simulates free convection.
    ml_index is the index of the depth of the surface mixed layer after adjustment,
    
    """

    stat_unstable = True
      
    while stat_unstable:
        
        d_diff = np.diff(d)
        if np.any(d_diff<0):
            stat_unstable=True
            first_inst_idx = np.flatnonzero(d_diff<0)[0]
            d0 = d
            (t, s, d, u, v) = mix5(t, s, d, u, v, first_inst_idx+1)       
        else:
            stat_unstable = False
            
    return t, s, d, u, v

def mix5(t, s, d, u, v, j):
    r"""
    mix temperature, salinity, momentum down to level j.
    
    """
    j = j+1
    t[:j] = np.mean(t[:j])
    s[:j] = np.mean(s[:j])
    d[:j] = sw.dens0(s[:j], t[:j])
    u[:j] = np.mean(u[:j])
    v[:j] = np.mean(v[:j])
    
    return t, s, d, u, v

def rot(u, v, ang):
    r"""
    This subroutine rotates the vector (u,v) through an angle, ang
    
    """
    r = (u+1j*v)*np.exp(1j*ang)
    u = r.real
    v = r.imag
    
    return u, v   

def bulk_mix(t, s, d, u, v, g, rb, nz, z, mld_idx):
    r"""
    sub-routine to do bulk richardson mixing
    
    """
    rvc = rb #critical rich number??
    
    for j in range(mld_idx, nz):
        h   = z[j]
        #it looks like density and velocity are mixed from the surface down to the ML depth
        dd  = (d[j]-d[0])/d[0]
        dv  = (u[j]-u[0])**2+(v[j]-v[0])**2
        if dv == 0:
            rv = np.inf
        else:
            rv = g*h*dd/dv

        if rv > rvc:
            break
        else:
            t, s, d, u, v = mix5(t, s, d, u, v, j)
            
    return t, s, d, u, v

def grad_mix(t, s, d, u, v, dz, g, rg, rb, nz, z, mld_idx):
    r"""
    This function performs the gradeint Richardson Number relaxation
    by mixing adjacent cells just enough to bring them to a new Richardson Number.
    Compute the gradeint Richardson Number, taking care to avoid dividing by
    zero in the mixed layer.  The numerical values of the minimum allowable
    density and velocity differences are entirely arbitrary, and should not
    effect the calculations (except that on some occasions they evidently have!)

    """

    rc = rg #critical rich. number
    j1 = 0
    j2 = nz-1
    j_range = np.arange(j1,j2)
    i = 0 #loop count
    #debug_here()
    r_min=rc
    while r_min>rc:
        #TODO: find a better way to do implement this loop
        r = np.zeros(len(j_range))
        
        for j in j_range:
            
            dd = (d[j+1]-d[j])/d[j]
            dv = (u[j+1]-u[j])**2+(v[j+1]-v[j])**2
            if dv==0 or dv<1e-10:
                r[j] = np.inf
            else:
                r[j] = g*dz*dd/dv                
         
        r_min = np.nanmin(r)
        j_min_idx = np.nanargmin(r)
            
        #Mix the cells j_min_idx and j_min_idx+1 that had the smallest Richardson Number
        t, s, d, u, v = stir(t, s, d, u, v, rc, r_min, j_min_idx)

        t, s, d, u, v = bulk_mix(t, s, d, u, v, g, rb, nz, z, mld_idx)

        
        #recompute the rich number over the part of the profile that has changed
        j1 = j_min_idx-2
        if j1 < 1:
             j1 = 0
        
        j2 = j_min_idx+2
        if j2 > nz-1:
             j2 = nz-1
             
        i+=1
                     
    return t, s, d, u, v
                
def stir(t, s, d, u, v, rc, r, j):
    r"""
    This subroutine mixes cells j and j+1 just enough so that
    the Richardson number after the mixing is brought up to
    the value rnew. In order to have this mixing process
    converge, rnew must exceed the critical value of the
    richardson number where mixing is presumed to start. If
    r critical = rc = 0.25 (the nominal value), and r = 0.20, then
    rnew = 0.3 would be reasonable. If r were smaller, then a
    larger value of rnew - rc is used to hasten convergence.
    This subroutine was modified by JFP in Sep 93 to allow for an
    aribtrary rc and to achieve faster convergence.
    
    """
    
    rcon = 0.02+(rc-r)/2
    rnew = rc+rcon/5.
    f = 1-r/rnew
    
    dt = (t[j+1]-t[j])*f/2.
    t[j+1] = t[j+1]-dt
    t[j] = t[j]+dt
    
    ds = (s[j+1]-s[j])*f/2.
    s[j+1] = s[j+1]-ds
    s[j] = s[j]+ds
    d[[j,j+1]] = sw.dens0(s[[j,j+1]], t[[j,j+1]])
    
    du = (u[j+1]-u[j])*f/2
    u[j+1] = u[j+1]-du
    u[j] = u[j]+du
    
    dv = (v[j+1]-v[j])*f/2
    v[j+1] = v[j+1]-dv
    v[j] = v[j]+dv
    
    return t, s, d, u, v
    
def diffuse(nz,a,diffusivity):
    
    "finite difference implementation of diffusion equation"
     
    #matlab code:
    #a(2:nz-1) = a(2:nz-1) + dstab*(a(1:nz-2) - 2*a(2:nz-1) + a(3:nz));
    
    a[1:nz-1] = a[1:nz-1] + diffusivity*(a[0:nz-2] - 2*a[1:nz-1] + a[2:nz]) 
    return a


