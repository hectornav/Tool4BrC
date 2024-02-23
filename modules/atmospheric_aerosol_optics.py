"""
Aerosol Optical Properties Analysis Module

This module is dedicated to the computation of optical properties of atmospheric aerosols. 
It leverages detailed physical properties of various aerosol types to calculate key optical 
metrics like extinction, scattering, and asymmetry parameters. The module is designed to 
support atmospheric science research, particularly in the fields of climate modeling and 
environmental studies where understanding aerosol behavior is crucial.

The module includes functions to calculate optical properties for a broad range of aerosol types,
including specific computations for brown carbon aerosols. It utilizes data from NetCDF files 
and integrates with other modules, such as `physics`, for comprehensive property analysis. 
This approach allows for flexible, accurate calculations that can be tailored to specific research needs.

Key Functions:
- calculate_optical_properties: Computes optical properties for a variety of aerosol types based on 
  relative humidity, type identifiers, and wavelength considerations.
- calculate_optical_properties4brc: Specialized function for calculating the optical properties of 
  brown carbon aerosols.

The module emphasizes versatility and precision, making it a valuable tool for researchers 
studying the impact of aerosols on atmospheric phenomena and climate change.
"""

import numpy as np
import netCDF4 as nc
from modules.physics import calculate_physical_properties, calculate_brown_carbon_properties, calculate4oatot, calculatePhysicalParameters
import os 
import sys

#print(os.getcwd())

def calculate_optical_properties(hamb, tids, wvln, **kwargs):
    """Calculate the optical properties of atmospheric aerosols.

    Args:
        hamb (ndarray): Array of float values representing the relative
            humidities for which the optical properties are to be calculated.
        tids (list): List of strings representing the names of the aerosol
            types for which the optical properties are to be calculated.
        wvln: wavelength in microns. as default 0.37 which corresponds to 370 nm

    Returns:
        tuple: A tuple of arrays containing the calculated optical properties
        of the aerosols.

    """
    ri_gfs_poa = kwargs.get('ri_gfs_poa', None)
    ri_gfs_soa = kwargs.get('ri_gfs_soa', None)
    ri_res_poa = kwargs.get('ri_res_poa', None)
    ri_res_soa = kwargs.get('ri_res_soa', None)
    ri_shp_poa = kwargs.get('ri_shp_poa', None)
    ri_shp_soa = kwargs.get('ri_shp_soa', None)
    ri_trf_poa = kwargs.get('ri_trf_poa', None)
    ri_trf_soa = kwargs.get('ri_trf_soa', None)
    ri_oth_poa = kwargs.get('ri_oth_poa', None)
    ri_oth_soa = kwargs.get('ri_oth_soa', None)
    
    wvln = wvln/1000                          #pass wavelength to mm
    phys = calculate_physical_properties(ri_gfs_poa, ri_gfs_soa, ri_res_poa, ri_res_soa,
                            ri_shp_poa, ri_shp_soa, ri_trf_poa, ri_trf_soa,
                            ri_oth_poa, ri_oth_soa)
         
    with nc.Dataset('core_data/nkxrnel.nc', 'r') as inp:
        ngrd = inp.variables['NGRD'][:]
        kgrd = inp.variables['KGRD'][:]
        xvnb = inp.variables['XVNB'][:]
        qext = inp.variables['QEXT'][:,:,:]
        qsca = inp.variables['QSCA'][:,:,:]
        asym = inp.variables['ASYM'][:,:,:]

    rvnb = xvnb * wvln / (2.0 * np.pi)
    hlev = np.array([0.0, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99])

    nwat = phys['WATE'][0]
    kwat = phys['WATE'][1]

    hamb = hamb.flatten()
    H = hamb.shape[0]
    I = len(tids)

    eout, sout, gout, dout, rout, vout = [], [], [], [], [], []
    for h in np.arange(H):
        eens, sens, gens, dens, rens, vens = [], [], [], [], [], []
        for i in np.arange(I):
            rgeo = phys[tids[i]][0]
            sdev = phys[tids[i]][1]
            rmin = phys[tids[i]][2]
            rmax = phys[tids[i]][3]
            real = phys[tids[i]][4]
            imag = phys[tids[i]][5]
            dens.append(phys[tids[i]][6])
            alfa = phys[tids[i]][7]
            beta = phys[tids[i]][8]

            # dry geometric integration
            r0 = np.abs(rvnb - rmin).argmin()
            r1 = np.abs(rvnb - rmax).argmin() + 1
            rgrd = rvnb[r0:r1]
            rrri = rvnb[r0:r1]**3
            logn = (1.0/(np.sqrt(2.0*np.pi)*np.log(sdev)))*(1.0/rvnb[r0:r1])*np.exp(-0.5*(np.log(rvnb[r0:r1]/rgeo)/np.log(sdev))**2)
            drrr = np.trapz(rrri*logn, x=rgrd)

            # hygroscopic growth factors
            h0 = np.searchsorted(hlev, hamb[h]) - 1
            h1 = np.searchsorted(hlev, hamb[h])
            MA = (alfa[h1] - alfa[h0])/(hlev[h1] - hlev[h0])
            alfh = alfa[h0] + MA*(hamb[h] - hlev[h0])
            MB = (beta[h1] - beta[h0])/(hlev[h1] - hlev[h0])
            beth = beta[h0] + MB*(hamb[h] - hlev[h0])

            # wet geometric integration
            rgeo = alfh*rgeo**beth
            sdev = sdev**beth
            rmin = alfh*rmin**beth
            rmax = alfh*rmax**beth

            r0 = np.abs(rvnb - rmin).argmin()
            r1 = np.abs(rvnb - rmax).argmin() + 1
            rgrd = rvnb[r0:r1]
            rrii = rvnb[r0:r1]**2
            prri = np.pi*rvnb[r0:r1]**2
            rrri = rvnb[r0:r1]**3
            logn = (1.0/(np.sqrt(2.0*np.pi)*np.log(sdev)))*(1.0/rvnb[r0:r1])*np.exp(-0.5*(np.log(rvnb[r0:r1]/rgeo)/np.log(sdev))**2)
            wwrr = np.trapz(rrii*logn, x=rgrd)
            area = np.trapz(prri*logn, x=rgrd)
            wrrr = np.trapz(rrri*logn, x=rgrd)

            # wet geometric ratios
            rens.append(wrrr/wwrr)
            vens.append(wrrr/drrr)
            vfra = drrr/wrrr

            # wet optical integration
            real = real*vfra + nwat*(1.0 - vfra)
            imag = imag*vfra + kwat*(1.0 - vfra)

            n0 = np.searchsorted(ngrd, real) - 1
            n1 = np.searchsorted(ngrd, real)
            k0 = np.searchsorted(kgrd, imag) - 1
            k1 = np.searchsorted(kgrd, imag)

            qe00, qe01, qe10, qe11 = qext[n0,k0,r0:r1], qext[n0,k1,r0:r1], qext[n1,k0,r0:r1], qext[n1,k1,r0:r1]
            qs00, qs01, qs10, qs11 = qsca[n0,k0,r0:r1], qsca[n0,k1,r0:r1], qsca[n1,k0,r0:r1], qsca[n1,k1,r0:r1]
            qg00, qg01, qg10, qg11 = asym[n0,k0,r0:r1], asym[n0,k1,r0:r1], asym[n1,k0,r0:r1], asym[n1,k1,r0:r1]

            enqe = np.array([np.trapz(qe00*prri*logn, x=rgrd)/area, np.trapz(qe01*prri*logn, x=rgrd)/area,
                                np.trapz(qe10*prri*logn, x=rgrd)/area, np.trapz(qe11*prri*logn, x=rgrd)/area])
            enqs = np.array([np.trapz(qs00*prri*logn, x=rgrd)/area, np.trapz(qs01*prri*logn, x=rgrd)/area,
                                np.trapz(qs10*prri*logn, x=rgrd)/area, np.trapz(qs11*prri*logn, x=rgrd)/area])
            enqg = np.array([np.trapz(qg00*qs00*prri*logn, x=rgrd)/(enqs[0]*area), np.trapz(qg01*qs01*prri*logn, x=rgrd)/(enqs[1]*area),
                                np.trapz(qg10*qs10*prri*logn, x=rgrd)/(enqs[2]*area), np.trapz(qg11*qs11*prri*logn, x=rgrd)/(enqs[3]*area)])

            # kernel grid interpolation
            ddnk = (ngrd[n1] - ngrd[n0])*(kgrd[k1] - kgrd[k0])
            wt00 = (ngrd[n1] - real)*(kgrd[k1] - imag)/ddnk
            wt01 = (ngrd[n1] - real)*(imag - kgrd[k0])/ddnk
            wt10 = (real - ngrd[n0])*(kgrd[k1] - imag)/ddnk
            wt11 = (real - ngrd[n0])*(imag - kgrd[k0])/ddnk
            wgth = np.array([wt00, wt01, wt10, wt11])

            eens.append(np.average(enqe, weights=wgth))
            sens.append(np.average(enqs, weights=wgth))
            gens.append(np.average(enqg, weights=wgth))

        eout.append(eens)
        sout.append(sens)
        gout.append(gens)
        dout.append(dens)
        rout.append(rens)
        vout.append(vens)

    #make a dictionary with tids and the output arrays
    optical_parameters = {}
    for i in range(len(tids)):
        optical_parameters[tids[i].lower()] = {'qext': eout[0][i], 'qsca': sout[0][i], 
                                       'asym': gout[0][i], 'dens': dout[0][i], 
                                       'reff': rout[0][i], 'vphi': vout[0][i]}
    return optical_parameters
    #return np.array(eout), np.array(sout), np.array(gout), np.array(dout), np.array(rout), np.array(vout)

def calculate_optical_properties4brc(hamb, tids, wvln, **kwargs):
    """Calculate the optical properties of brown carbon aerosols, primary and secondary.

    Args:
        hamb (ndarray): Array of float values representing the relative
            humidities for which the optical properties are to be calculated.
        tids (list): List of strings representing the names of the aerosol
            types for which the optical properties are to be calculated.
        wvln: wavelength in microns. as default 0.37 which corresponds to 370 nm

    Returns:
        tuple: A tuple of arrays containing the calculated optical properties
        of the aerosols.calculate_optical_properties4brc

    """
    ri_brc_strng = kwargs.get('ri_brc_strng')
    ri_brc_blchd = kwargs.get('ri_brc_blchd')

    wvln = wvln/1000                          #pass wavelength to mm

    phys = calculate_brown_carbon_properties(ri_brc_strng, ri_brc_blchd)
         
    with nc.Dataset('core_data/nkxrnel.nc', 'r') as inp:
        ngrd = inp.variables['NGRD'][:]
        kgrd = inp.variables['KGRD'][:]
        xvnb = inp.variables['XVNB'][:]
        qext = inp.variables['QEXT'][:,:,:]
        qsca = inp.variables['QSCA'][:,:,:]
        asym = inp.variables['ASYM'][:,:,:]

    rvnb = xvnb * wvln / (2.0 * np.pi)
    hlev = np.array([0.0, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99])

    nwat = phys['WATE'][0]
    kwat = phys['WATE'][1]

    hamb = hamb.flatten()
    H = hamb.shape[0]
    I = len(tids)

    eout, sout, gout, dout, rout, vout = [], [], [], [], [], []
    for h in np.arange(H):
        eens, sens, gens, dens, rens, vens = [], [], [], [], [], []
        for i in np.arange(I):
            rgeo = phys[tids[i]][0]
            sdev = phys[tids[i]][1]
            rmin = phys[tids[i]][2]
            rmax = phys[tids[i]][3]
            real = phys[tids[i]][4]
            imag = phys[tids[i]][5]
            dens.append(phys[tids[i]][6])
            alfa = phys[tids[i]][7]
            beta = phys[tids[i]][8]

            # dry geometric integration
            r0 = np.abs(rvnb - rmin).argmin()
            r1 = np.abs(rvnb - rmax).argmin() + 1
            rgrd = rvnb[r0:r1]
            rrri = rvnb[r0:r1]**3
            logn = (1.0/(np.sqrt(2.0*np.pi)*np.log(sdev)))*(1.0/rvnb[r0:r1])*np.exp(-0.5*(np.log(rvnb[r0:r1]/rgeo)/np.log(sdev))**2)
            drrr = np.trapz(rrri*logn, x=rgrd)

            # hygroscopic growth factors
            h0 = np.searchsorted(hlev, hamb[h]) - 1
            h1 = np.searchsorted(hlev, hamb[h])
            MA = (alfa[h1] - alfa[h0])/(hlev[h1] - hlev[h0])
            alfh = alfa[h0] + MA*(hamb[h] - hlev[h0])
            MB = (beta[h1] - beta[h0])/(hlev[h1] - hlev[h0])
            beth = beta[h0] + MB*(hamb[h] - hlev[h0])

            # wet geometric integration
            rgeo = alfh*rgeo**beth
            sdev = sdev**beth
            rmin = alfh*rmin**beth
            rmax = alfh*rmax**beth

            r0 = np.abs(rvnb - rmin).argmin()
            r1 = np.abs(rvnb - rmax).argmin() + 1
            rgrd = rvnb[r0:r1]
            rrii = rvnb[r0:r1]**2
            prri = np.pi*rvnb[r0:r1]**2
            rrri = rvnb[r0:r1]**3
            logn = (1.0/(np.sqrt(2.0*np.pi)*np.log(sdev)))*(1.0/rvnb[r0:r1])*np.exp(-0.5*(np.log(rvnb[r0:r1]/rgeo)/np.log(sdev))**2)
            wwrr = np.trapz(rrii*logn, x=rgrd)
            area = np.trapz(prri*logn, x=rgrd)
            wrrr = np.trapz(rrri*logn, x=rgrd)

            # wet geometric ratios
            rens.append(wrrr/wwrr)
            vens.append(wrrr/drrr)
            vfra = drrr/wrrr

            # wet optical integration
            real = real*vfra + nwat*(1.0 - vfra)
            imag = imag*vfra + kwat*(1.0 - vfra)

            n0 = np.searchsorted(ngrd, real) - 1
            n1 = np.searchsorted(ngrd, real)
            k0 = np.searchsorted(kgrd, imag) - 1
            k1 = np.searchsorted(kgrd, imag)

            qe00, qe01, qe10, qe11 = qext[n0,k0,r0:r1], qext[n0,k1,r0:r1], qext[n1,k0,r0:r1], qext[n1,k1,r0:r1]
            qs00, qs01, qs10, qs11 = qsca[n0,k0,r0:r1], qsca[n0,k1,r0:r1], qsca[n1,k0,r0:r1], qsca[n1,k1,r0:r1]
            qg00, qg01, qg10, qg11 = asym[n0,k0,r0:r1], asym[n0,k1,r0:r1], asym[n1,k0,r0:r1], asym[n1,k1,r0:r1]

            enqe = np.array([np.trapz(qe00*prri*logn, x=rgrd)/area, np.trapz(qe01*prri*logn, x=rgrd)/area,
                                np.trapz(qe10*prri*logn, x=rgrd)/area, np.trapz(qe11*prri*logn, x=rgrd)/area])
            enqs = np.array([np.trapz(qs00*prri*logn, x=rgrd)/area, np.trapz(qs01*prri*logn, x=rgrd)/area,
                                np.trapz(qs10*prri*logn, x=rgrd)/area, np.trapz(qs11*prri*logn, x=rgrd)/area])
            enqg = np.array([np.trapz(qg00*qs00*prri*logn, x=rgrd)/(enqs[0]*area), np.trapz(qg01*qs01*prri*logn, x=rgrd)/(enqs[1]*area),
                                np.trapz(qg10*qs10*prri*logn, x=rgrd)/(enqs[2]*area), np.trapz(qg11*qs11*prri*logn, x=rgrd)/(enqs[3]*area)])

            # kernel grid interpolation
            ddnk = (ngrd[n1] - ngrd[n0])*(kgrd[k1] - kgrd[k0])
            wt00 = (ngrd[n1] - real)*(kgrd[k1] - imag)/ddnk
            wt01 = (ngrd[n1] - real)*(imag - kgrd[k0])/ddnk
            wt10 = (real - ngrd[n0])*(kgrd[k1] - imag)/ddnk
            wt11 = (real - ngrd[n0])*(imag - kgrd[k0])/ddnk
            wgth = np.array([wt00, wt01, wt10, wt11])

            eens.append(np.average(enqe, weights=wgth))
            sens.append(np.average(enqs, weights=wgth))
            gens.append(np.average(enqg, weights=wgth))

        eout.append(eens)
        sout.append(sens)
        gout.append(gens)
        dout.append(dens)
        rout.append(rens)
        vout.append(vens)

    #make a dictionary with tids and the output arrays
    optical_parameters = {}
    for i in range(len(tids)):
        optical_parameters[tids[i].lower()] = {'qext': eout[0][i], 'qsca': sout[0][i], 
                                       'asym': gout[0][i], 'dens': dout[0][i], 
                                       'reff': rout[0][i], 'vphi': vout[0][i]}
    return optical_parameters
    #return np.array(eout), np.array(sout), np.array(gout), np.array(dout), np.array(rout), np.array(vout)

def calculate_optical_properties4oa(hamb, tids, wvln, **kwargs):
    """Calculate the optical properties of organic aerosols, primary and secondary.

    Args:
        hamb (ndarray): Array of float values representing the relative
            humidities for which the optical properties are to be calculated.
        tids (list): List of strings representing the names of the aerosol
            types for which the optical properties are to be calculated.
        wvln: wavelength in microns. as default 0.37 which corresponds to 370 nm

    Returns:
        tuple: A tuple of arrays containing the calculated optical properties
        of the aerosols.calculate_optical_properties4brc

    """
    ri_oa = kwargs.get('ri_oa')

    wvln = wvln/1000                          #pass wavelength to mm

    phys = calculate4oatot(ri_oa)

         
    with nc.Dataset('core_data/nkxrnel.nc', 'r') as inp:
        ngrd = inp.variables['NGRD'][:]
        kgrd = inp.variables['KGRD'][:]
        xvnb = inp.variables['XVNB'][:]
        qext = inp.variables['QEXT'][:,:,:]
        qsca = inp.variables['QSCA'][:,:,:]
        asym = inp.variables['ASYM'][:,:,:]

    rvnb = xvnb * wvln / (2.0 * np.pi)
    hlev = np.array([0.0, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99])

    nwat = phys['WATE'][0]
    kwat = phys['WATE'][1]

    hamb = hamb.flatten()
    H = hamb.shape[0]
    I = len(tids)

    eout, sout, gout, dout, rout, vout = [], [], [], [], [], []
    for h in np.arange(H):
        eens, sens, gens, dens, rens, vens = [], [], [], [], [], []
        for i in np.arange(I):
            rgeo = phys[tids[i]][0]
            sdev = phys[tids[i]][1]
            rmin = phys[tids[i]][2]
            rmax = phys[tids[i]][3]
            real = phys[tids[i]][4]
            imag = phys[tids[i]][5]
            dens.append(phys[tids[i]][6])
            alfa = phys[tids[i]][7]
            beta = phys[tids[i]][8]

            # dry geometric integration
            r0 = np.abs(rvnb - rmin).argmin()
            r1 = np.abs(rvnb - rmax).argmin() + 1
            rgrd = rvnb[r0:r1]
            rrri = rvnb[r0:r1]**3
            logn = (1.0/(np.sqrt(2.0*np.pi)*np.log(sdev)))*(1.0/rvnb[r0:r1])*np.exp(-0.5*(np.log(rvnb[r0:r1]/rgeo)/np.log(sdev))**2)
            drrr = np.trapz(rrri*logn, x=rgrd)

            # hygroscopic growth factors
            h0 = np.searchsorted(hlev, hamb[h]) - 1
            h1 = np.searchsorted(hlev, hamb[h])
            MA = (alfa[h1] - alfa[h0])/(hlev[h1] - hlev[h0])
            alfh = alfa[h0] + MA*(hamb[h] - hlev[h0])
            MB = (beta[h1] - beta[h0])/(hlev[h1] - hlev[h0])
            beth = beta[h0] + MB*(hamb[h] - hlev[h0])

            # wet geometric integration
            rgeo = alfh*rgeo**beth
            sdev = sdev**beth
            rmin = alfh*rmin**beth
            rmax = alfh*rmax**beth

            r0 = np.abs(rvnb - rmin).argmin()
            r1 = np.abs(rvnb - rmax).argmin() + 1
            rgrd = rvnb[r0:r1]
            rrii = rvnb[r0:r1]**2
            prri = np.pi*rvnb[r0:r1]**2
            rrri = rvnb[r0:r1]**3
            logn = (1.0/(np.sqrt(2.0*np.pi)*np.log(sdev)))*(1.0/rvnb[r0:r1])*np.exp(-0.5*(np.log(rvnb[r0:r1]/rgeo)/np.log(sdev))**2)
            wwrr = np.trapz(rrii*logn, x=rgrd)
            area = np.trapz(prri*logn, x=rgrd)
            wrrr = np.trapz(rrri*logn, x=rgrd)

            # wet geometric ratios
            rens.append(wrrr/wwrr)
            vens.append(wrrr/drrr)
            vfra = drrr/wrrr

            # wet optical integration
            real = real*vfra + nwat*(1.0 - vfra)
            imag = imag*vfra + kwat*(1.0 - vfra)

            n0 = np.searchsorted(ngrd, real) - 1
            n1 = np.searchsorted(ngrd, real)
            k0 = np.searchsorted(kgrd, imag) - 1
            k1 = np.searchsorted(kgrd, imag)

            qe00, qe01, qe10, qe11 = qext[n0,k0,r0:r1], qext[n0,k1,r0:r1], qext[n1,k0,r0:r1], qext[n1,k1,r0:r1]
            qs00, qs01, qs10, qs11 = qsca[n0,k0,r0:r1], qsca[n0,k1,r0:r1], qsca[n1,k0,r0:r1], qsca[n1,k1,r0:r1]
            qg00, qg01, qg10, qg11 = asym[n0,k0,r0:r1], asym[n0,k1,r0:r1], asym[n1,k0,r0:r1], asym[n1,k1,r0:r1]

            enqe = np.array([np.trapz(qe00*prri*logn, x=rgrd)/area, np.trapz(qe01*prri*logn, x=rgrd)/area,
                                np.trapz(qe10*prri*logn, x=rgrd)/area, np.trapz(qe11*prri*logn, x=rgrd)/area])
            enqs = np.array([np.trapz(qs00*prri*logn, x=rgrd)/area, np.trapz(qs01*prri*logn, x=rgrd)/area,
                                np.trapz(qs10*prri*logn, x=rgrd)/area, np.trapz(qs11*prri*logn, x=rgrd)/area])
            enqg = np.array([np.trapz(qg00*qs00*prri*logn, x=rgrd)/(enqs[0]*area), np.trapz(qg01*qs01*prri*logn, x=rgrd)/(enqs[1]*area),
                                np.trapz(qg10*qs10*prri*logn, x=rgrd)/(enqs[2]*area), np.trapz(qg11*qs11*prri*logn, x=rgrd)/(enqs[3]*area)])

            # kernel grid interpolation
            ddnk = (ngrd[n1] - ngrd[n0])*(kgrd[k1] - kgrd[k0])
            wt00 = (ngrd[n1] - real)*(kgrd[k1] - imag)/ddnk
            wt01 = (ngrd[n1] - real)*(imag - kgrd[k0])/ddnk
            wt10 = (real - ngrd[n0])*(kgrd[k1] - imag)/ddnk
            wt11 = (real - ngrd[n0])*(imag - kgrd[k0])/ddnk
            wgth = np.array([wt00, wt01, wt10, wt11])

            eens.append(np.average(enqe, weights=wgth))
            sens.append(np.average(enqs, weights=wgth))
            gens.append(np.average(enqg, weights=wgth))

        eout.append(eens)
        sout.append(sens)
        gout.append(gens)
        dout.append(dens)
        rout.append(rens)
        vout.append(vens)

    #make a dictionary with tids and the output arrays
    optical_parameters = {}
    for i in range(len(tids)):
        optical_parameters[tids[i].lower()] = {'qext': eout[0][i], 'qsca': sout[0][i], 
                                       'asym': gout[0][i], 'dens': dout[0][i], 
                                       'reff': rout[0][i], 'vphi': vout[0][i]}
    
    return optical_parameters


def calculateOpticalPropertiesNoSecondary(hamb, tids, wvln, **kwargs):
    """Calculate the optical properties of atmospheric aerosols.

    Args:
        hamb (ndarray): Array of float values representing the relative
            humidities for which the optical properties are to be calculated.
        tids (list): List of strings representing the names of the aerosol
            types for which the optical properties are to be calculated.
        wvln: wavelength in microns. as default 0.37 which corresponds to 370 nm

    Returns:
        tuple: A tuple of arrays containing the calculated optical properties
        of the aerosols.

    """
    ri_gfas = kwargs.get('ri_gfas', None)
    ri_resi = kwargs.get('ri_resi', None)
    ri_ship = kwargs.get('ri_ship', None)
    ri_traf = kwargs.get('ri_traf', None)
    ri_othr = kwargs.get('ri_othr', None)
    
    wvln = wvln/1000                          #pass wavelength to mm
    phys = calculatePhysicalParameters(ri_gfas, ri_resi, ri_ship, ri_traf,
                                       ri_othr)
         
    with nc.Dataset('core_data/nkxrnel.nc', 'r') as inp:
        ngrd = inp.variables['NGRD'][:]
        kgrd = inp.variables['KGRD'][:]
        xvnb = inp.variables['XVNB'][:]
        qext = inp.variables['QEXT'][:,:,:]
        qsca = inp.variables['QSCA'][:,:,:]
        asym = inp.variables['ASYM'][:,:,:]

    rvnb = xvnb * wvln / (2.0 * np.pi)
    hlev = np.array([0.0, 0.50, 0.70, 0.80, 0.90, 0.95, 0.99])

    nwat = phys['WATE'][0]
    kwat = phys['WATE'][1]

    hamb = hamb.flatten()
    H = hamb.shape[0]
    I = len(tids)

    eout, sout, gout, dout, rout, vout = [], [], [], [], [], []
    for h in np.arange(H):
        eens, sens, gens, dens, rens, vens = [], [], [], [], [], []
        for i in np.arange(I):
            rgeo = phys[tids[i]][0]
            sdev = phys[tids[i]][1]
            rmin = phys[tids[i]][2]
            rmax = phys[tids[i]][3]
            real = phys[tids[i]][4]
            imag = phys[tids[i]][5]
            dens.append(phys[tids[i]][6])
            alfa = phys[tids[i]][7]
            beta = phys[tids[i]][8]

            # dry geometric integration
            r0 = np.abs(rvnb - rmin).argmin()
            r1 = np.abs(rvnb - rmax).argmin() + 1
            rgrd = rvnb[r0:r1]
            rrri = rvnb[r0:r1]**3
            logn = (1.0/(np.sqrt(2.0*np.pi)*np.log(sdev)))*(1.0/rvnb[r0:r1])*np.exp(-0.5*(np.log(rvnb[r0:r1]/rgeo)/np.log(sdev))**2)
            drrr = np.trapz(rrri*logn, x=rgrd)

            # hygroscopic growth factors
            h0 = np.searchsorted(hlev, hamb[h]) - 1
            h1 = np.searchsorted(hlev, hamb[h])
            MA = (alfa[h1] - alfa[h0])/(hlev[h1] - hlev[h0])
            alfh = alfa[h0] + MA*(hamb[h] - hlev[h0])
            MB = (beta[h1] - beta[h0])/(hlev[h1] - hlev[h0])
            beth = beta[h0] + MB*(hamb[h] - hlev[h0])

            # wet geometric integration
            rgeo = alfh*rgeo**beth
            sdev = sdev**beth
            rmin = alfh*rmin**beth
            rmax = alfh*rmax**beth

            r0 = np.abs(rvnb - rmin).argmin()
            r1 = np.abs(rvnb - rmax).argmin() + 1
            rgrd = rvnb[r0:r1]
            rrii = rvnb[r0:r1]**2
            prri = np.pi*rvnb[r0:r1]**2
            rrri = rvnb[r0:r1]**3
            logn = (1.0/(np.sqrt(2.0*np.pi)*np.log(sdev)))*(1.0/rvnb[r0:r1])*np.exp(-0.5*(np.log(rvnb[r0:r1]/rgeo)/np.log(sdev))**2)
            wwrr = np.trapz(rrii*logn, x=rgrd)
            area = np.trapz(prri*logn, x=rgrd)
            wrrr = np.trapz(rrri*logn, x=rgrd)

            # wet geometric ratios
            rens.append(wrrr/wwrr)
            vens.append(wrrr/drrr)
            vfra = drrr/wrrr

            # wet optical integration
            real = real*vfra + nwat*(1.0 - vfra)
            imag = imag*vfra + kwat*(1.0 - vfra)

            n0 = np.searchsorted(ngrd, real) - 1
            n1 = np.searchsorted(ngrd, real)
            k0 = np.searchsorted(kgrd, imag) - 1
            k1 = np.searchsorted(kgrd, imag)

            qe00, qe01, qe10, qe11 = qext[n0,k0,r0:r1], qext[n0,k1,r0:r1], qext[n1,k0,r0:r1], qext[n1,k1,r0:r1]
            qs00, qs01, qs10, qs11 = qsca[n0,k0,r0:r1], qsca[n0,k1,r0:r1], qsca[n1,k0,r0:r1], qsca[n1,k1,r0:r1]
            qg00, qg01, qg10, qg11 = asym[n0,k0,r0:r1], asym[n0,k1,r0:r1], asym[n1,k0,r0:r1], asym[n1,k1,r0:r1]

            enqe = np.array([np.trapz(qe00*prri*logn, x=rgrd)/area, np.trapz(qe01*prri*logn, x=rgrd)/area,
                                np.trapz(qe10*prri*logn, x=rgrd)/area, np.trapz(qe11*prri*logn, x=rgrd)/area])
            enqs = np.array([np.trapz(qs00*prri*logn, x=rgrd)/area, np.trapz(qs01*prri*logn, x=rgrd)/area,
                                np.trapz(qs10*prri*logn, x=rgrd)/area, np.trapz(qs11*prri*logn, x=rgrd)/area])
            enqg = np.array([np.trapz(qg00*qs00*prri*logn, x=rgrd)/(enqs[0]*area), np.trapz(qg01*qs01*prri*logn, x=rgrd)/(enqs[1]*area),
                                np.trapz(qg10*qs10*prri*logn, x=rgrd)/(enqs[2]*area), np.trapz(qg11*qs11*prri*logn, x=rgrd)/(enqs[3]*area)])

            # kernel grid interpolation
            ddnk = (ngrd[n1] - ngrd[n0])*(kgrd[k1] - kgrd[k0])
            wt00 = (ngrd[n1] - real)*(kgrd[k1] - imag)/ddnk
            wt01 = (ngrd[n1] - real)*(imag - kgrd[k0])/ddnk
            wt10 = (real - ngrd[n0])*(kgrd[k1] - imag)/ddnk
            wt11 = (real - ngrd[n0])*(imag - kgrd[k0])/ddnk
            wgth = np.array([wt00, wt01, wt10, wt11])

            eens.append(np.average(enqe, weights=wgth))
            sens.append(np.average(enqs, weights=wgth))
            gens.append(np.average(enqg, weights=wgth))

        eout.append(eens)
        sout.append(sens)
        gout.append(gens)
        dout.append(dens)
        rout.append(rens)
        vout.append(vens)

    #make a dictionary with tids and the output arrays
    optical_parameters = {}
    for i in range(len(tids)):
        optical_parameters[tids[i].lower()] = {'qext': eout[0][i], 'qsca': sout[0][i], 
                                       'asym': gout[0][i], 'dens': dout[0][i], 
                                       'reff': rout[0][i], 'vphi': vout[0][i]}
    return optical_parameters
