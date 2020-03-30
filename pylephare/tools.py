
import numpy as np
from astropy import constants

# --------------------------- #
# - Conversion Tools        - #
# --------------------------- #
def flux_to_mag(flux, dflux, wavelength=None, zp=None, inhz=False):
    """ Converts fluxes (erg/s/cm2/A) into AB or zp magnitudes


    Parameters
    ----------
    flux, fluxerr: [float or array]
        flux and its error 

    wavelength: [float or array] -optional-
        central wavelength [in AA] of the photometric filter.
        // Ignored if inhz=True //

    zp: [float or array] -optional-
        zero point of for flux;
        // Ignored if inhz=True //
        // Ignored if wavelength is provided //

    inhz:
        set to true if the flux (and flux) are given in erg/s/cm2/Hz
        (False means in erg/s/cm2/AA)
        
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
            
    mag_ab = -2.5*np.log10(flux*wavelength**2) + zp
    if dflux is None:
        return mag_ab, None
    
    dmag_ab = -2.5/np.log(10) * dflux / flux
    return mag_ab, dmag_ab

def mag_to_flux(mag, magerr=None, wavelength=None, zp=None, inhz=False):
    """ converts magnitude into flux

    Parameters
    ----------
    mag: [float or array]
        AB magnitude(s)

    magerr: [float or array] -optional-
        magnitude error if any

    wavelength: [float or array] -optional-
        central wavelength [in AA] of the photometric filter.
        // Ignored if inhz=True //

    zp: [float or array] -optional-
        zero point of for flux;
        // Ignored if inhz=True //
        // Ignored if wavelength is provided //

    inhz:
        set to true if the flux (and flux) are given in erg/s/cm2/Hz
        (False means in erg/s/cm2/AA)

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


def flux_aa_to_hz(flux_aa, wavelength):
    """ """
    return flux_aa * (wavelength**2 / constants.c.to("AA/s").value)
    
def flux_hz_to_aa(flux_hz, wavelength):
    """ """
    return flux_hz / (wavelength**2 / constants.c.to("AA/s").value)
