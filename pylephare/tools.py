
import numpy as np
from astropy import constants

# --------------------------- #
# - Conversion Tools        - #
# --------------------------- #
def flux_to_mag(flux, dflux=None, wavelength=None, zp=None, inhz=False):
    """
    Convert fluxes into AB or zp magnitudes.

    Parameters
    ----------
    flux : [float or array]
        Flux(es).
    
    dflux : [float or array or None]
        Flux uncertainties if any.
        Default is None.

    wavelength : [float or array]
        // Ignored if inhz=True //
        Central wavelength(s) [in AA] of the photometric filter(s).
        Default is None.

    zp : [float or array]
        // Ignored if inhz=True //
        // Ignored if wavelength is provided //
        Zero point of the flux(es).
        Default is None.

    inhz : [bool]
        Set to true if the flux(es) (and uncertainties) is/are given in erg/s/cm2/Hz
        False means in erg/s/cm2/AA.
        Default is False.
        
    Returns
    -------
    - float or array (if magerr is None)
    - float or array, float or array (if magerr provided)
    """
    if inhz:
        zp = -48.598 # instaad of -48.60 such that hz to aa is correct
        wavelength = 1
    else:
        if zp is None and wavelength is None:
            raise ValueError("'zp' or 'wavelength' must be provided.")
        if zp is None:
            zp = -2.406 
        else:
            wavelength=1
            
    mag_ab = -2.5*np.log10(flux*wavelength**2) + zp
    if dflux is None:
        return mag_ab, None
    
    dmag_ab = -2.5/np.log(10) * dflux / flux
    return mag_ab, dmag_ab

def mag_to_flux(mag, magerr=None, wavelength=None, zp=None, inhz=False):
    """
    Convert AB or zp magnitudes into fluxes

    Parameters
    ----------
    mag : [float or array]
        AB magnitude(s)

    magerr : [float or array or None]
        Magnitude uncertainties if any.
        Default is None.

    wavelength : [float or array]
        // Ignored if inhz=True //
        Central wavelength(s) [in AA] of the photometric filter(s).
        Default is None.

    zp : [float or array]
        // Ignored if inhz=True //
        // Ignored if wavelength is provided //
        Zero point of the flux(es).
        Default is None.

    inhz : [bool]
        Set to true if the flux(es) (and uncertainties) is/are returned in erg/s/cm2/Hz
        False means in erg/s/cm2/AA.
        Default is False.

    Returns
    -------
    - float or array (if magerr is None)
    - float or array, float or array (if magerr provided)
    """
    if inhz:
        zp = -48.598 # instaad of -48.60 such that hz to aa is correct
        wavelength = 1
    else:
        if zp is None and wavelength is None:
            raise ValueError("zp or wavelength must be provided")
        if zp is None:
            zp = -2.406 
        else:
            wavelength=1

    flux = 10**(-(mag-zp)/2.5) / wavelength**2
    if magerr is None:
        return flux, None
    
    dflux = np.abs(flux*(-magerr/2.5*np.log(10))) # df/f = dcount/count
    return flux, dflux


def flux_aa_to_hz(flux, wavelength):
    """
    Convert fluxes from erg/s/cm2/AA to erg/s/cm2/Hz.
    
    Parameters
    ----------
    flux_aa : [float or array]
        Flux(es) in erg/s/cm2/AA.
    
    wavelength : [float or array]
        Central wavelength(s) [in AA] of the photometric filter(s).
    
    
    Returns
    -------
    float or array
    """
    return flux * (wavelength**2 / constants.c.to("AA/s").value)
    
def flux_hz_to_aa(flux, wavelength):
    """
    Convert fluxes from erg/s/cm2/Hz to erg/s/cm2/AA.
    
    Parameters
    ----------
    flux_aa : [float or array]
        Flux(es) in erg/s/cm2/Hz.
    
    wavelength : [float or array]
        Central wavelength(s) [in AA] of the photometric filter(s).
    
    
    Returns
    -------
    float or array
    """
    return flux / (wavelength**2 / constants.c.to("AA/s").value)
