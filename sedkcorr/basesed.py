
import numpy as np
import pandas
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
import pkg_resources

from astropy import units
from sncosmo import bandpasses

from pyifu import spectroscopy
from propobject import BaseObject




MAG_0_INFO = {"FUV":np.log10(1.40e-15) + 0.4 * 18.82,
    "NUV":np.log10(2.06e-16) + 0.4 * 20.08,
    "u":-8.056, "g":-8.326, "r":-8.555, "i":-8.732, "z":-8.882}

FILTER_BANDS = ["FUV", "NUV", "u", "g", "r", "i", "z"]

LBDA_INFO = {"FUV":1539, "NUV":2316, "u":3562, "g":4719, "r":6185, "i":7500, "z":8961}

COLOR_INFO = {"FUV":"xkcd:purple", "NUV":"xkcd:violet", "u":"xkcd:blue", "g":"xkcd:green", "r":"xkcd:red", "i":"xkcd:cyan", "z":"xkcd:magenta"}



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

def sed_mag_to_flux(mag):
    return np.power(10, (mag + 48.585)/(-2.5))

def sed_mag_to_flux_err(mag, mag_err, mag0_err=0.005):
    return (0.4 * np.log(10) * mag_to_flux(mag))*np.sqrt(mag_err**2 + mag0_err**2)

def sed_flux_to_mag(flux):
    return -2.5 * np.log10(flux) - 48.585

def band_flux_to_mag(flux, band):
    flux = flux_nu_to_flux_lbda(flux, LBDA_INFO[band])
    return -2.5 * (np.log10(flux) - MAG_0_INFO[band])

def band_flux_to_mag_err(flux, flux_err, band):
    flux = flux_nu_to_flux_lbda(flux, LBDA_INFO[band])
    flux_err = flux_nu_to_flux_lbda(flux_err, LBDA_INFO[band])
    return (2.5 / (np.log(10) * flux)) * flux_err

def band_mag_to_flux(mag, band):
    flux = 10**(MAG_0_INFO[band] - 0.4*mag)
    return flux_lbda_to_flux_nu(flux, LBDA_INFO[band])

def band_mag_to_flux_err(mag, mag_err, band):
    flux = band_mag_to_flux(mag, band)
    flux_err = 0.4 * np.log(10) * flux_nu_to_flux_lbda(flux, LBDA_INFO[band]) * mag_err
    return flux_lbda_to_flux_nu(flux_err, LBDA_INFO[band])

def flux_nu_to_flux_lbda(flux_nu, band_lbda):
    flux_nu = flux_nu * units.erg / units.s / units.cm**2 / units.Hz
    flux_lbda = flux_nu.to(units.erg / units.s / units.cm**2 / units.AA, equivalencies=units.spectral_density(band_lbda * units.AA))
    return flux_lbda.value

def flux_lbda_to_flux_nu(flux_lbda, band_lbda):
    flux_lbda = flux_lbda * units.erg / units.s / units.cm**2 / units.AA
    flux_nu = flux_lbda.to(units.erg / units.s / units.cm**2 / units.Hz, equivalencies=units.spectral_density(band_lbda * units.AA))
    return flux_nu.value




# ----------------- #
#                   #
#     SED class     #
#                   #
# ------------------#

class SED( BaseObject ):
    """
    
    """
    
    PROPERTIES         = ["data_sed", "list_bands", "data_meas", "z", "filter_bandpass"]
    SIDE_PROPERTIES    = ["filter_bandpass_path"]
    DERIVED_PROPERTIES = ["data_sed_shifted", "data_kcorr"]
    
    def set_data(self, data_sed, data_meas, z, list_bands=None, **kwargs):
        """
            
        """
        if type(data_sed) is pandas.DataFrame:
            self._properties["data_sed"] = data_sed
        elif type(data_sed) is dict:
            self._properties["data_sed"] = pandas.DataFrame(data_sed)
        else:
            raise TypeError("data_sed must be a DataFrame or a dict")
        self.data_sed["flux"] = sed_mag_to_flux(self.data_sed["mag"])

        self._properties["list_bands"] = [band for band in data_meas] if list_bands is None else list_bands
        
        self._properties["data_meas"] = data_meas
        ######### Option to set : mag or flux ##############
        for band in self.list_bands:
            self.data_meas[band]["flux"] = band_mag_to_flux(self.data_meas[band]["mag"], band)
            self.data_meas[band]["flux.err"] = band_mag_to_flux_err(self.data_meas[band]["mag"], self.data_meas[band]["mag.err"], band)
        
        self._properties["z"] = z
    
    def set_filter_bandpass_path(self):
        """
        Set the path to filter bandpasses data, included in the package.
        
        
        Returns
        -------
        Void
        """
        self._side_properties["filter_bandpass_path"] = {band:pkg_resources.resource_filename(__name__, "filter_bandpass/SLOAN_SDSS."+band+".dat")
                                                         for band in ["u", "g", "r", "i", "z"]}
        self.filter_bandpass_path["FUV"] = pkg_resources.resource_filename(__name__, "filter_bandpass/GALEX_GALEX.FUV.dat")
        self.filter_bandpass_path["NUV"] = pkg_resources.resource_filename(__name__, "filter_bandpass/GALEX_GALEX.NUV.dat")
    
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
        opt_bands = FILTER_BANDS if opt_bands is None else opt_bands if type(opt_bands)==list else [opt_bands]
        list_sdss_bands = ["u", "g", "r", "i", "z"]
        if from_sncosmo:
            self._properties["filter_bandpass"] = {band:pandas.DataFrame({"lbda":bandpasses.get_bandpass("sdss"+band).wave,
                                                                      "trans":bandpasses.get_bandpass("sdss"+band).trans})
                                                   for band in opt_bands if band in list_sdss_bands}
        else:
            self._properties["filter_bandpass"] = {band:pandas.read_table(self.filter_bandpass_path[band], sep=' ', names=["lbda", "trans"])
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
                                                  for band in FILTER_BANDS}
        for band in FILTER_BANDS:
            self.data_kcorr[band]["mag"] = band_flux_to_mag(self.data_kcorr[band]["flux"], band)
    
    def show(self, y_plot="flux", sed_shifted=True, plot_bandpasses=False, plot_filter_points=True, xlim=(None, None), ylim=(None, None), savefile=None):
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
            for band in FILTER_BANDS:
                ax.plot(self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"]*(ax_ylim[1]-ax_ylim[0]) + ax_ylim[0],
                        ls='-', marker='', color=COLOR_INFO[band], label="_nolegend_")
        
        if plot_filter_points:
            for band in (FILTER_BANDS if sed_shifted else self.list_bands):
                x_point = LBDA_INFO[band]
                y_point = self.data_kcorr[band][y_plot] if sed_shifted else self.data_meas[band][y_plot]
                y_err_point = self.data_kcorr[band][y_plot+".err"] if sed_shifted else self.data_meas[band][y_plot+".err"]
                ax.errorbar(x_point, y_point, yerr=y_err_point, ls='', marker='o', color=COLOR_INFO[band], label=band)
        
            ax.set_xlabel(r"$\lambda$ [\AA]", fontsize="large")
            ax.set_ylabel(("${m}_{AB}$" if y_plot=="mag" else r"${f}_{\nu}$ $[erg.{s}^{-1}.{cm}^{-2}.{Hz}^{-1}]$"), fontsize="large")
            ax.legend(loc="upper right", ncol=1)
                
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
        if self._properties["filter_bandpass"] is None:
            self.set_filter_bandpass()
        return self._properties["filter_bandpass"]
    
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
        if self._side_properties["filter_bandpass_path"] is None:
            self.set_filter_bandpass_path()
        return self._side_properties["filter_bandpass_path"]
    
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
