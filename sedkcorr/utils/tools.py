
import numpy as np
import pandas
from astropy import units

# Filled thanks to http://svo2.cab.inta-csic.es/svo/theory/fps3/
FILTER_BANDS = {"galex.FUV":{"lbda":1545.8,
                             "color":"xkcd:purple",
                             "ZP":3619.9, #np.log10(1.40e-15) + 0.4 * 18.82,
                             "bandpass_file":"filter_bandpass/GALEX_GALEX.FUV.dat",
                             "prospector_name":"galex_FUV",
                             "error_scale":0.052,
                             "lephare_name":"galex/FUV.pb"},
                "galex.NUV":{"lbda":2344.9,
                             "color":"xkcd:violet",
                             "ZP":3801.4, #np.log10(2.06e-16) + 0.4 * 20.08,
                             "bandpass_file":"filter_bandpass/GALEX_GALEX.NUV.dat",
                             "prospector_name":"galex_NUV",
                             "error_scale":0.026,
                             "lephare_name":"galex/NUV.pb"},
                "sdss.u":{"lbda":3561.8,
                          "color":"xkcd:blue",
                          "ZP":3767.2, #-8.056,
                          "bandpass_file":"filter_bandpass/SLOAN_SDSS.u.dat",
                          "prospector_name":"sdss_u0",
                          "error_scale":0.05,
                          "lephare_name":"sdss/up.pb"},
                "sdss.g":{"lbda":4718.9,
                          "color":"xkcd:green",
                          "ZP":3631, #-8.326,
                          "bandpass_file":"filter_bandpass/SLOAN_SDSS.g.dat",
                          "prospector_name":"sdss_g0",
                          "error_scale":0.02,
                          "lephare_name":"sdss/gp.pb"},
                "sdss.r":{"lbda":6185.2,
                          "color":"xkcd:red",
                          "ZP":3631, #-8.555,
                          "bandpass_file":"filter_bandpass/SLOAN_SDSS.r.dat",
                          "prospector_name":"sdss_r0",
                          "error_scale":0.02,
                          "lephare_name":"sdss/rp.pb"},
                "sdss.i":{"lbda":7499.7,
                          "color":"xkcd:cyan",
                          "ZP":3631, #-8.732,
                          "bandpass_file":"filter_bandpass/SLOAN_SDSS.i.dat",
                          "prospector_name":"sdss_i0",
                          "error_scale":0.02,
                          "lephare_name":"sdss/ip.pb"},
                "sdss.z":{"lbda":8961.5,
                          "color":"xkcd:magenta",
                          "ZP":3564.7, #-8.882,
                          "bandpass_file":"filter_bandpass/SLOAN_SDSS.z.dat",
                          "prospector_name":"sdss_z0",
                          "error_scale":0.03,
                          "lephare_name":"sdss/zp.pb"},
                "ps1.g":{"lbda":4866.5,
                         "color":"xkcd:green",
                         "ZP":3631,
                         "bandpass_file":"filter_bandpass/PAN-STARRS_PS1.g.dat",
                         "prospector_name":"",
                         "error_scale":0.0,
                         "lephare_name":"ps1/g_ps.pb"},
                "ps1.r":{"lbda":6214.6,
                         "color":"xkcd:red",
                         "ZP":3631,
                         "bandpass_file":"filter_bandpass/PAN-STARRS_PS1.r.dat",
                         "prospector_name":"",
                         "error_scale":0.0,
                         "lephare_name":"ps1/r_ps.pb"},
                "ps1.i":{"lbda":7544.6,
                         "color":"xkcd:cyan",
                         "ZP":3631,
                         "bandpass_file":"filter_bandpass/PAN-STARRS_PS1.i.dat",
                         "prospector_name":"",
                         "error_scale":0.0,
                         "lephare_name":"ps1/i_ps.pb"},
                "ps1.z":{"lbda":8679.5,
                         "color":"xkcd:magenta",
                         "ZP":3631,
                         "bandpass_file":"filter_bandpass/PAN-STARRS_PS1.z.dat",
                         "prospector_name":"",
                         "error_scale":0.0,
                         "lephare_name":"ps1/z_ps.pb"},
                "ps1.y":{"lbda":9633.3,
                         "color":"xkcd:pink",
                         "ZP":3631,
                         "bandpass_file":"filter_bandpass/PAN-STARRS_PS1.y.dat",
                         "prospector_name":"",
                         "error_scale":0.0,
                         "lephare_name":"ps1/y_ps.pb"}
                }

def lbda_z_to_z0(lbda, z):
    """
    Shift wavelength to redshift zero.
    lbda_out = lbda_in / ( 1 + z )

    Parameters
    ----------
    lbda : [float or np.array]
        Wavelength.

    z : [float or np.array]
        Redshift.


    Returns
    -------
    float or np.array
    """
    return lbda / ( 1 + z )

def lbda_z0_to_z(lbda, z):
    """
    Shift wavelength from redshift zero to the observed one.
    lbda_out = lbda_out * ( 1 + z )
    
    Parameters
    ----------
    lbda : [float or np.array]
        Wavelength.
    
    z : [float or np.array]
        Redshift.
    
    
    Returns
    -------
    float or np.array
    """
    return lbda * ( 1 + z )
    
def flux_z_to_z0(flux, z):
    """
    Shift flux to redshift zero.
    flux_out = flux_in / ( ( 1 + z ) ** 3 )
    
    Parameters
    ----------
    flux : [float or np.array]
        Flux.
    
    z : [float or np.array]
        Redshift.
    
    
    Returns
    -------
    float or np.array
    """
    return flux / ( ( 1 + z )**3 )

def flux_z0_to_z(flux, z):
    """
    Shift flux from redshift zero to the observed one.
    flux_out = flux_in * ( ( 1 + z ) ** 3 )
    
    Parameters
    ----------
    flux : [float or np.array]
        Flux.
    
    z : [float or np.array]
        Redshift.
    
    
    Returns
    -------
    float or np.array
    """
    return flux * ( ( 1 + z )**3 )

def mag_ab_zeropoint(flux, lbda):
    """
    Return the AB magnitude zero point.
    
    Parameters
    ----------
    flux : [float]
        Zeropoint flux in Jansky.
    
    lbda : [float]
        Wavelength.
    
    
    Returns
    -------
    float
    """
    flux *= 10**-23
    zp = convert_flux_unit(flux=flux, lbda=lbda, unit_in="Hz", unit_out="AA")
    return np.log10(zp)

def mag_to_flux(mag, mag_err=None, band=None, flux_unit="Hz", opt_mAB0=False):
    """
    Convert magnitude to flux.
    Return the flux and its error.
    
    Parameters
    ----------
    mag : [float or np.array]
        Magnitude.
    
    mag_err : [float or np.array or None]
        Magnitude error.
        Default is None.
    
    band : [string or np.array]
        The output flux is converted in the input unit with the wavelength based on this given 'band' filter.
        You can also directly give the wavelength.
    
    flux_unit : [string]
        Define the output flux unit :
        - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
        - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
        - "mgy" : mgy (mgy = maggies)
    
    Options
    -------
    opt_mAB0 : [bool]
        If True, take in account the filter calibration (AB zero magnitude, depend on the instrument (SDSS, Galex, etc.).
        Only works if the given 'band' is a filter name "project.band" (ex: "sdss.u", ...).
        Default is False.
    
    
    Returns
    -------
    float or np.array, float or np.array
    """
    if len(np.atleast_1d(mag)) > 1 and (np.all(np.atleast_1d(mag_err)==0) or mag_err is None):
        mag_err = None
    if opt_mAB0:
        flux_out = 10**(KCorrection.mag_ab_zeropoint(FILTER_BANDS[band]["ZP"], FILTER_BANDS[band]["lbda"]) - 0.4*mag)
        unit_in = "AA"
    else:
        flux_out = 10**((mag + 48.585)/(-2.5))
        unit_in = "Hz"
    flux_err_out = None if mag_err is None else (0.4 * np.log(10) * flux_out) * mag_err**2
    
    if flux_unit in ("AA", "Hz", "mgy"):
        flux_out = convert_flux_unit(flux_out, FILTER_BANDS[band]["lbda"] if type(band)==str else band, unit_in=unit_in, unit_out=flux_unit)
        if flux_err_out is not None:
            flux_err_out = convert_flux_unit(flux_err_out, FILTER_BANDS[band]["lbda"] if type(band)==str else band,
                                             unit_in=unit_in, unit_out=flux_unit)
    else:
        raise ValueError("{} is not a valid flux unit.".format(flux_unit))
    
    if mag_err is None:
        flux_err_out = np.zeros(len(mag)) if len(np.atleast_1d(mag)) > 1 else 0.
    
    return flux_out, flux_err_out

def flux_to_mag(flux, flux_err=None, band=None, flux_unit="Hz", opt_mAB0=False):
    """
    Convert flux to magnitude.
    Return the magnitude and its error.
    
    Parameters
    ----------
    flux : [float or np.array]
        Flux.
    
    flux_err : [float or np.array or None]
        Flux error.
        Default is None.
    
    band : [string or np.array]
        The input flux is converted in a conversion convenient unit from the input unit with the wavelength based on this given 'band' filter.
        You can also directly give the wavelength.
    
    flux_unit : [string]
        Define the input flux unit :
        - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
        - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
        - "mgy" : mgy (mgy = maggies)
    
    Options
    -------
    opt_mAB0 : [bool]
        If True, take in account the filter calibration (AB zero magnitude), depending on the instrument (SDSS, Galex, etc.).
        Only works if the given 'band' is a filter name "project.band" (ex: "sdss.u", ...).
        Default is False.
    
    
    Returns
    -------
    float or np.array, float or np.array
    """
    if (flux_err == 0. or flux_err is None) and len(np.atleast_1d(flux)) > 1:
        flux_err = None
    if opt_mAB0:
        unit_out = "AA"
    else:
        unit_out = "Hz"
    if flux_unit in ("AA", "Hz", "mgy"):
        flux = convert_flux_unit(flux, FILTER_BANDS[band]["lbda"] if type(band)==str else band, unit_in=flux_unit, unit_out=unit_out)
        if flux_err is not None:
            flux_err = convert_flux_unit(flux_err, FILTER_BANDS[band]["lbda"] if type(band)==str else band, unit_in=flux_unit, unit_out=unit_out)
    else:
        raise ValueError("{} is not a valid flux unit.".format(flux_unit))

    if opt_mAB0:
        mag_out = -2.5 * (np.log10(flux) - mag_ab_zeropoint(FILTER_BANDS[band]["ZP"], FILTER_BANDS[band]["lbda"]))
    else:
        mag_out = -2.5 * np.log10(flux) - 48.585
    
    if flux_err is None:
        mag_err_out = np.zeros(len(flux)) if len(np.atleast_1d(flux)) > 1 else 0.
    else:
        mag_err_out = (2.5 / np.log(10)) * (flux_err / flux)
        
    return mag_out, mag_err_out

def convert_flux_unit(flux, lbda, unit_in, unit_out):
    """
    Convert the flux unit.
    
    Parameters
    ----------
    flux : [float or np.array]
        Input flux.
    
    lbda : [float or np.array]
        Wavelength of the flux. Have to be the same size than 'flux'.
    
    unit_in : [string]
        Unit of the input flux :
        - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
        - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
        - "mgy" : mgy (mgy = maggies)
    
    unit_in : [string]
        Unit of the output flux :
        - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
        - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
        - "mgy" : mgy (mgy = maggies)
    
    
    Returns
    -------
    float or np.array
    """
    unit_base = units.erg / units.s / units.cm**2
    flux = np.asarray(flux) if len(np.atleast_1d(flux)) > 1 else float(flux)
    lbda = np.asarray(lbda) if len(np.atleast_1d(lbda)) > 1 else float(lbda)
    equiv = units.spectral_density(lbda * units.AA)
    if unit_in not in ("Hz", "AA", "mgy"):
        raise ValueError("{} is not a valid flux unit.".format(unit_in))
    if unit_out not in ("Hz", "AA", "mgy"):
        raise ValueError("{} is not a valid flux unit.".format(unit_out))
    if unit_in == unit_out:
        return flux
    elif "mgy" not in (unit_in, unit_out):
        if unit_in == "Hz" and unit_out == "AA":
            unit_in, unit_out = unit_base / units.Hz, unit_base / units.AA
        elif unit_in == "AA" and unit_out == "Hz":
            unit_in, unit_out = unit_base / units.AA, unit_base / units.Hz
        return ((flux*unit_in).to(unit_out, equivalencies=equiv)).value
    elif unit_out == "mgy":
        if unit_in == "AA":
            flux = ((flux*unit_base/units.AA).to(unit_base/units.Hz, equivalencies=equiv)).value
        return flux / (3631e-23)
    elif unit_in == "mgy":
        flux = flux * (3631e-23) * unit_base / units.Hz
        if unit_out == "AA":
            flux = flux.to(unit_base/units.AA, equivalencies=equiv)
        return flux.value
    else:
        raise ValueError("-- ERROR 404 -- : There is a big problem in the conversion process.")
