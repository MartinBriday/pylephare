
import numpy as np
import pandas
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
import pkg_resources

from astropy import units
from sncosmo import bandpasses

from pyifu import spectroscopy
from propobject import BaseObject


LIST_BANDS = ["FUV", "NUV", "u", "g", "r", "i", "z"]
FILTER_BANDS = {"FUV":{"lbda":1539,
                       "color":"xkcd:purple",
                       "mAB0":np.log10(1.40e-15) + 0.4 * 18.82,
                       "context_id":0,
                       "prospector_name":"galex_FUV"},
                "NUV":{"lbda":2316,
                       "color":"xkcd:violet",
                       "mAB0":np.log10(2.06e-16) + 0.4 * 20.08,
                       "context_id":1,
                       "prospector_name":"galex_NUV"},
                "u":{"lbda":3562,
                     "color":"xkcd:blue",
                     "mAB0":-8.056,
                     "context_id":2,
                     "prospector_name":"sdss_u0"},
                "g":{"lbda":4719,
                     "color":"xkcd:green",
                     "mAB0":-8.326,
                     "context_id":3,
                     "prospector_name":"sdss_g0"},
                "r":{"lbda":6185,
                     "color":"xkcd:red",
                     "mAB0":-8.555,
                     "context_id":4,
                     "prospector_name":"sdss_r0"},
                "i":{"lbda":7500,
                     "color":"xkcd:cyan",
                     "mAB0":-8.732,
                     "context_id":5,
                     "prospector_name":"sdss_i0"},
                "z":{"lbda":8961,
                     "color":"xkcd:magenta",
                     "mAB0":-8.882,
                    "context_id":6,
                    "prospector_name":"sdss_z0"},
                }



# ----------------- #
#                   #
#    Photometric    #
#     Functions     #
#                   #
# ------------------#
def lbda_z0(lbda, z):
    return lbda/(1+z)

def flux_z0(flux, z):
    return flux/(1+z)**3

def mag_to_flux(mag, band=None, flux_unit="Hz"):
    if band is None:
        return 10**((mag + 48.585)/(-2.5))
    elif band in FILTER_BANDS:
        flux_lbda = 10**(FILTER_BANDS[band]["mAB0"] - 0.4*mag)
        if flux_unit in ("AA", "A", "Angstrom"):
            return flux_lbda
        elif flux_unit in ("Hz", "Herz"):
            return convert_flux_unit(flux_lbda, FILTER_BANDS[band]["lbda"], unit_in="AA", unit_out="Hz")
        elif flux_unit in ("mgy", "maggy"):
            return convert_flux_unit(flux_lbda, FILTER_BANDS[band]["lbda"], unit_in="AA", unit_out="mgy")
        else:
            raise ValueError("{} is not a valid flux unit.".format(flux_unit))
    else:
        raise ValueError("{} is not an existing filter band.".format(band))

def flux_to_mag(flux, band=None, flux_unit="Hz"):
    if band is None:
        return -2.5 * np.log10(flux) - 48.585
    elif band in FILTER_BANDS:
        if flux_unit in ("AA", "A", "Angstrom"):
            flux = flux
        elif flux_unit in ("Hz", "Herz"):
            flux = convert_flux_unit(flux, FILTER_BANDS[band]["lbda"], unit_in="Hz", unit_out="AA")
        elif flux_unit in ("mgy", "maggy"):
            flux =  convert_flux_unit(flux, FILTER_BANDS[band]["lbda"], unit_in="mgy", unit_out="AA")
        else:
            raise ValueError("{} is not a valid flux unit.".format(flux_unit))
        return -2.5 * (np.log10(flux) - FILTER_BANDS[band]["mAB0"])
    else:
        raise ValueError("{} is not an existing filter band.".format(band))

def sed_mag_to_flux_err(mag, mag_err, mag0_err=0.005):
    return (0.4 * np.log(10) * mag_to_flux(mag))*np.sqrt(mag_err**2 + mag0_err**2)

def band_flux_to_mag(flux, band):
    flux = flux_nu_to_flux_lbda(flux, FILTER_BANDS[band]["lbda"])
    return -2.5 * (np.log10(flux) - FILTER_BANDS[band]["mAB0"])

def band_flux_to_mag_err(flux, flux_err, band):
    flux = flux_nu_to_flux_lbda(flux, FILTER_BANDS[band]["lbda"])
    flux_err = flux_nu_to_flux_lbda(flux_err, FILTER_BANDS[band]["lbda"])
    return (2.5 / (np.log(10) * flux)) * flux_err

def band_mag_to_flux(mag, band):
    flux = 10**(FILTER_BANDS[band]["mAB0"] - 0.4*mag)
    return flux_lbda_to_flux_nu(flux, FILTER_BANDS[band]["lbda"])

def band_mag_to_flux_err(mag, mag_err, band):
    flux = band_mag_to_flux(mag, band)
    flux_err = 0.4 * np.log(10) * flux_nu_to_flux_lbda(flux, FILTER_BANDS[band]["lbda"]) * mag_err
    return flux_lbda_to_flux_nu(flux_err, FILTER_BANDS[band]["lbda"])

def flux_nu_to_flux_lbda(flux_nu, band_lbda):
    flux_nu = flux_nu * units.erg / units.s / units.cm**2 / units.Hz
    flux_lbda = flux_nu.to(units.erg / units.s / units.cm**2 / units.AA, equivalencies=units.spectral_density(band_lbda * units.AA))
    return flux_lbda.value

def flux_lbda_to_flux_nu(flux_lbda, band_lbda):
    flux_lbda = flux_lbda * units.erg / units.s / units.cm**2 / units.AA
    flux_nu = flux_lbda.to(units.erg / units.s / units.cm**2 / units.Hz, equivalencies=units.spectral_density(band_lbda * units.AA))
    return flux_nu.value

def flux_nu_to_mgy(flux_nu):
    return flux_nu / (3631e-23)




# ----------------- #
#                   #
#     SED class     #
#                   #
# ------------------#

class SED( BaseObject ):
    """
    
    """
    
    PROPERTIES         = ["data_sed", "data_meas", "z"]
    SIDE_PROPERTIES    = ["list_bands", "filter_bandpass"]
    DERIVED_PROPERTIES = ["data_sed_shifted", "data_kcorr"]
    
    def __init__(self, **kwargs):
        """
            
        """
        if kwargs != {}:
            self.set_data(**kwargs)
    
    def set_data(self, **kwargs):
        """
            
        """
        self.set_data_sed(**kwargs)
        self.set_data_meas(**kwargs)
            
    def set_data_sed(self, data_sed=None, **extras):
        """
        
        """
        if type(data_sed) is pandas.DataFrame:
            self._properties["data_sed"] = data_sed
        elif type(data_sed) is dict:
            self._properties["data_sed"] = pandas.DataFrame(data_sed)
        else:
            raise TypeError("data_sed must be a DataFrame or a dict")
        
        if "mag" in self.data_sed.keys():
            self.data_sed["flux"] = sed_mag_to_flux(self.data_sed["mag"])
        elif "flux" in self.data_sed.keys():
            self.data_sed["mag"] = sed_flux_to_mag(self.data_sed["flux"])

    def set_data_meas(self, data_meas=None, z=None, **extras):
        """
        
        """
        self._properties["data_meas"] = data_meas
        
        ######### Option to set : mag or flux ##############
        for band in self.list_bands:
            if "mag" in self.data_meas[band].keys() and "flux" not in self.data_meas[band].keys():
                self.data_meas[band]["flux"] = band_mag_to_flux(self.data_meas[band]["mag"], band)
                self.data_meas[band]["flux.err"] = band_mag_to_flux_err(self.data_meas[band]["mag"], self.data_meas[band]["mag.err"], band)
            elif "flux" in self.data_meas[band].keys() and "mag" not in self.data_meas[band].keys():
                self.data_meas[band]["mag"] = band_flux_to_mag(self.data_meas[band]["flux"], band)
                self.data_meas[band]["mag.err"] = band_flux_to_mag_err(self.data_meas[band]["flux"], self.data_meas[band]["flux.err"], band)
        
        self._properties["z"] = z
    
    def get_filter_bandpass_path(self):
        """
        Set the path to filter bandpasses data, included in the package.
        
        
        Returns
        -------
        dict
        """
        dict_path = {band:pkg_resources.resource_filename(__name__, "filter_bandpass/SLOAN_SDSS."+band+".dat")
                     for band in ["u", "g", "r", "i", "z"]}
        dict_path["FUV"] = pkg_resources.resource_filename(__name__, "filter_bandpass/GALEX_GALEX.FUV.dat")
        dict_path["NUV"] = pkg_resources.resource_filename(__name__, "filter_bandpass/GALEX_GALEX.NUV.dat")
        return dict_path
    
    def set_filter_bandpass(self, opt_bands=None, from_sncosmo=False):
        """
        Load filter bandpasses data.
        
        Options
        -------
        opt_bands : [list[string] or None]
        Optionnal list of bands to load filter bandpass data.
        
        from_sncosmo : [bool]
        If True, filter bandpass data will be loaded from sncosmo package.
        Else, loaded from tables included in the packages, downloaded from "http://svo2.cab.inta-csic.es/svo/theory/fps3/"
        Default is False.
        
        
        Returns
        -------
        Void
        """
        opt_bands = [b for b in LIST_BANDS] if opt_bands is None else opt_bands if type(opt_bands)==list else [opt_bands]
        list_sdss_bands = ["u", "g", "r", "i", "z"]
        if from_sncosmo:
            self._side_properties["filter_bandpass"] = {band:pandas.DataFrame({"lbda":bandpasses.get_bandpass("sdss"+band).wave,
                                                                               "trans":bandpasses.get_bandpass("sdss"+band).trans})
                                                        for band in opt_bands if band in list_sdss_bands}
        else:
            self._side_properties["filter_bandpass"] = {band:pandas.read_table(self.filter_bandpass_path[band], sep=' ', names=["lbda", "trans"])
                                                        for band in opt_bands if band in list_sdss_bands}
        if "FUV" in opt_bands:
            self.filter_bandpass["FUV"] = pandas.read_table(self.filter_bandpass_path["FUV"], sep=' ', names=["lbda", "trans"])
        if "NUV" in opt_bands:
            self.filter_bandpass["NUV"] = pandas.read_table(self.filter_bandpass_path["NUV"], sep=' ', names=["lbda", "trans"])

    def shift_sed(self):
        """
        Shift the SED spectrum to redshift zero.
        
        
        Returns
        -------
        Void
        """
        self._derived_properties["data_sed_shifted"] = pandas.DataFrame({"lbda":lbda_z0(self.data_sed["lbda"], self.z),
                                                                     "flux":flux_z0(self.data_sed["flux"], self.z)})
        self.data_sed_shifted["mag"] = pandas.Series(sed_flux_to_mag(self.data_sed_shifted["flux"]))
    
    def k_correction(self):
        """
        Recover the integrated flux from every filter bands from the shifted SED.
        Then convert them into magnitudes.
        
        
        Returns
        -------
        Void
        """
        self._derived_properties["data_kcorr"] = {band:{"flux":spectroscopy.synthesize_photometry(self.data_sed_shifted["lbda"], self.data_sed_shifted["flux"],
                                                                                                  self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"])}
                                                  for band in LIST_BANDS}
        for band in LIST_BANDS:
            self.data_kcorr[band]["mag"] = band_flux_to_mag(self.data_kcorr[band]["flux"], band)
    
    def show(self, y_plot="flux", sed_shifted=True, plot_bandpasses=False, plot_filter_points=True,
             xlim=(None, None), ylim=(None, None), xscale="linear", yscale="linear", savefile=None):
        """
        Plot method.
        
        Parameters
        ----------
        y_plot : [string]
            Choice to plot "flux" or "mag".
        
        sed_shifted : [bool]
            If True, the SED is shifted to redshift zero.
            Else, the LePhare fitted SED is plotted.
        
        Options
        -------
        plot_bandpasses : [bool]
            If True, plot the filter bandpass transmitions.
        
        plot_filter_points : [bool]
            If True, plot the filter points with errors, either in flux or in magnitude.
        
        xlim : [tuple[float or None]]
            Set the limits on the x axis.
            If (None, None), the figure has free x axis limits.
        
        ylim : [tuple[float or None]]
            Set the limits on the y axis.
            If (None, None), the figure has free y axis limits.
        
        savefile : [string or None]
            If None, the figure won't be saved.
            To save it, input a path directory + filename.
        
        Returns
        -------
        dict(fig, ax)
        """
        x_sed = self.data_sed_shifted["lbda"] if sed_shifted else self.data_sed["lbda"]
        y_sed = self.data_sed_shifted[y_plot] if sed_shifted else self.data_sed[y_plot]
        
        fig, ax = plt.subplots()
        opt_sed = {"ls":"-", "marker":"", "color":"0.4"}
        ax.plot(x_sed, y_sed, label="_nolegend_", **opt_sed)
        
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax_ylim = ax.get_ylim()
        if y_plot=="flux" and ylim[0] is None:
            ax_ylim = (0., ax_ylim[1])
        
        ax.axhline(ax_ylim[0], color="black", zorder=5)
        
        if plot_bandpasses:
            for band in LIST_BANDS:
                ax.plot(self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"]*(ax_ylim[1]-ax_ylim[0]) + ax_ylim[0],
                        ls='-', marker='', color=FILTER_BANDS[band]["color"], label="_nolegend_")
        
        if plot_filter_points:
            for band in (LIST_BANDS if sed_shifted else self.list_bands):
                x_point = FILTER_BANDS[band]["lbda"]
                y_point = self.data_kcorr[band][y_plot] if sed_shifted else self.data_meas[band][y_plot]
                y_err_point = self.data_kcorr[band][y_plot+".err"] if sed_shifted else self.data_meas[band][y_plot+".err"]
                ax.errorbar(x_point, y_point, yerr=y_err_point, ls='', marker='o', color=FILTER_BANDS[band]["color"], label=band)
        
            ax.set_xlabel(r"$\lambda$ [\AA]", fontsize="large")
            ax.set_ylabel(("${m}_{AB}$" if y_plot=="mag" else r"${f}_{\nu}$ $[erg.{s}^{-1}.{cm}^{-2}.{Hz}^{-1}]$"), fontsize="large")
            ax.legend(loc="upper right", ncol=1)

        ax.set_xscale(xscale)
        if y_plot=="flux":
            ax.set_yscale(yscale)
                
        if savefile is not None:
            fig.savefig(savefile)
                        
        return {"ax":ax, "fig":fig}


    # ================ #
    #  Properties      #
    # ================ #
    @property
    def data_sed(self):
        """ DataFrame of the the SED fit data """
        return self._properties["data_sed"]
    
    @property
    def filter_bandpass(self):
        """ Dictionnary of DataFrame of the bandpass data for every filter band """
        if self._side_properties["filter_bandpass"] is None:
            self.set_filter_bandpass()
        return self._side_properties["filter_bandpass"]
    
    @property
    def data_meas(self):
        """ Table like of the measurements """
        if self._properties["data_meas"] is None:
            self._properties["data_meas"] = {}
        return self._properties["data_meas"]
    
    @property
    def z(self):
        """ SNeIa host redshift """
        return self._properties["z"]
    
    @property
    def filter_bandpass_path(self):
        """ Dictionnary of paths for every filter badpass data """
        return self.get_filter_bandpass_path()
    
    @property
    def list_bands(self):
        """ List of the filter bands used in SED fitting """
        return self._side_properties["list_bands"]
    
    @property
    def data_sed_shifted(self):
        """ Dataframe of the shifted to redshift zero SED """
        if self._derived_properties["data_sed_shifted"] is None:
            self.shift_sed()
        return self._derived_properties["data_sed_shifted"]
    
    @property
    def data_kcorr(self):
        """ Table of k corrected data """
        if self._derived_properties["data_kcorr"] is None:
            self.k_correction()
        return self._derived_properties["data_kcorr"]
