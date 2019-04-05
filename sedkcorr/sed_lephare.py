
import pandas as pd
import numpy as np
from astropy import units
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
import os
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
#  Functions        #
#                   #
# ------------------#
def lbda_z0(lbda, z):
    return lbda/(1+z)

def flux_z0(flux, z):
    return flux*(1+z)**3

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






class SED( BaseObject ):
    """
    LePhare SED fits
    """
    
    PROPERTIES         = ["sed_data", "filter_bandpass", "meas_data", "z"]
    SIDE_PROPERTIES    = ["filter_bandpass_path", "list_bands"]
    DERIVED_PROPERTIES = ["sed_shifted", "kcorr_data"]
    
    def set_sed_data(self, sed_index=None, sed_data=None, sed_dir=None):
        """
        Load the SED spectrum.
        
        Parameters
        ----------
        sed_index : [int or None]
            Object SED file index.
            A SED directory is needed to be usable (sed_dir)
            You can not use both sed_index and sed_data.
        
        sed_data : [pandas.DataFrame or None]
            Object SED data input in a DataFrame.
            You can not use both sed_data and sed_index.
        
        sed_dir : [string]
            Path to the SED files folder.
            Needed if sed_index used.
        
        
        Returns
        -------
        Void
        """
        if sed_index is not None and sed_data is None:
            if sed_dir is None:
                raise ValueError("You did not input the sprectrum directory.")
            sed_index = str(sed_index)
            while len(sed_index)<9:
                sed_index = "0"+sed_index
            sed_filename = "Id"+sed_index+".spec"
            self._properties["sed_data"] = pd.read_table(os.path.expanduser(sed_dir+sed_filename), skiprows=20, names=["lbda", "mag"], sep="  ", engine="python", nrows=1057)
        elif sed_index is None and sed_data is not None:
            self._properties["sed_data"] = sed_data
        else:
            raise ValueError("You must input one (and only one) data set option : the spectrum index or a DataFrame with columns = ['lbda', 'mag'].")
        self.sed_data["flux"] = sed_mag_to_flux(self.sed_data["mag"])
    
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
            self._properties["filter_bandpass"] = {band:pd.DataFrame({"lbda":bandpasses.get_bandpass("sdss"+band).wave, 
                                                                      "trans":bandpasses.get_bandpass("sdss"+band).trans}) 
                                                   for band in opt_bands if band in list_sdss_bands}
        else:
            self._properties["filter_bandpass"] = {band:pd.read_table(self.filter_bandpass_path[band], sep=' ', names=["lbda", "trans"])
                                                   for band in opt_bands if band in list_sdss_bands}
        if "FUV" in opt_bands:
            self.filter_bandpass["FUV"] = pd.read_table(self.filter_bandpass_path["FUV"], sep=' ', names=["lbda", "trans"])
        if "NUV" in opt_bands:
            self.filter_bandpass["NUV"] = pd.read_table(self.filter_bandpass_path["NUV"], sep=' ', names=["lbda", "trans"])
    
    def context_filters(self, context):
        """
        Return a list of the concerned filter bands relative to the given context.
        
        Parameters
        ----------
        context : [int]
            LePhare type context, it defines the used filter bands for the SED fitting.
        
        
        Returns
        -------
        list(string)
        """
        idx = []
        for ii in range(len(FILTER_BANDS)-1,-1,-1):
            if (context - 2**ii) >= 0:
                context = context - 2**ii
                idx.append(ii)
        return [FILTER_BANDS[ii] for ii in reversed(idx)]
    
    def set_meas_data(self, meas_data=None, z=None, col_syntax=["mag_band", "mag_band_err"], list_bands=None):
        """
        Set the host redshift and the measured magnitudes for every filter bands used in SED fitting.
        
        Parameters
        ----------
        meas_data : [table like]
            Table like (eg : DataFrame line) of the measurements.
        
        z : [float or None]
            Redshift of the SNeIa host.
            If None, the redshift is supposed to be in the meas_data table under the name "Z-SPEC".
            Default is None.
        
        col_syntax : [list[string]]
            Syntax of measurements and errors column names in the meas_data table.
            Replace the filter band in the column names with the word "band" (eg: ["mag_band", "mag_band_err"]).
        
        Options
        -------
        list_bands : [list[string] or None]
            List of the filter bands used in SED fitting.
            If None, the LePhare context set list_bands. The context is supposed to be in the meas_data table under the name "CONTEXT".
        
        
        Returns
        -------
        Void
        """
        self._side_properties["list_bands"] = self.context_filters(meas_data["CONTEXT"]) if list_bands is None else list_bands
        z = z if z is not None else meas_data["Z-SPEC"]
        col_syntax = ["band", "band.err"] if col_syntax is None else col_syntax
        
        self._properties["meas_data"] = {band:{"mag":meas_data[col_syntax[0].replace("band",band)], 
                                               "mag.err":meas_data[col_syntax[1].replace("band",band)], 
                                               "flux":band_mag_to_flux(meas_data[col_syntax[0].replace("band",band)], band), 
                                               "flux.err":band_mag_to_flux_err(meas_data[col_syntax[0].replace("band",band)], 
                                                                               meas_data[col_syntax[1].replace("band",band)], band)}
                                         for band in self.list_bands}
        
        self._properties["z"] = z
    
    def shift_sed(self):
        """
        Shift the SED spectrum to redshift zero.
        
        
        Returns
        -------
        Void
        """
        self._derived_properties["sed_shifted"] = pd.DataFrame({"lbda":lbda_z0(self.sed_data["lbda"], self.z), 
                                                                "flux":flux_z0(self.sed_data["flux"], self.z)})
        self.sed_shifted["mag"] = pd.Series(sed_flux_to_mag(self.sed_shifted["flux"]))
    
    def k_correction(self):
        """
        Recover the integrated flux from every filter bands from the shifted SED.
        Then convert them into magnitudes.
        
        
        Returns
        -------
        Void
        """
        self._derived_properties["kcorr_data"] = {band:{"flux":spectroscopy.synthesize_photometry(self.sed_shifted["lbda"], self.sed_shifted["flux"], 
                                                                                                  self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"])}
                                                  for band in FILTER_BANDS}
        for band in self.list_bands:
            self.kcorr_data[band]["mag.err"] = self.meas_data[band]["mag.err"]
            self.kcorr_data[band]["flux.err"] = self.meas_data[band]["flux.err"]
        
        for band in FILTER_BANDS:
            self.kcorr_data[band]["mag"] = band_flux_to_mag(self.kcorr_data[band]["flux"], band)
            if band not in self.list_bands:
                self.kcorr_data[band]["flux.err"] = 0.
                self.kcorr_data[band]["mag.err"] = 0.
        
    
    def show(self, y_plot="flux", sed_shifted=True, plot_bandpasses=False, plot_filter_points=True, xlim=(None, None), ylim=(None, None), save_fig=False, save_dir=None):
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
        
        save_fig : [bool]
            If True, the plot will be saved in the given directory.
        
        save_dir : [string or None]
            Directory where the figure will be saved.
            If None, the default value is the Desktop (Unix users...).
        
        Returns
        -------
        
        """
        x_sed = self.sed_shifted["lbda"] if sed_shifted else self.sed_data["lbda"]
        y_sed = self.sed_shifted[y_plot] if sed_shifted else self.sed_data[y_plot]
        
        fig, ax = plt.subplots()
        opt_sed = {"ls":"-", "marker":"", "color":"0.4"}
        ax.plot(x_sed, y_sed, label="_nolegend_", **opt_sed)
        
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax_ylim = ax.get_ylim()
        if y_plot=="flux":
            ax_ylim = (0., ax_ylim[1])
            
        ax.axhline(ax_ylim[0], color="black", zorder=5)
            
        if plot_bandpasses:
            for band in FILTER_BANDS:
                ax.plot(self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"]*(ax_ylim[1]-ax_ylim[0]) + ax_ylim[0], 
                        ls='-', marker='', color=COLOR_INFO[band])
        
        if plot_filter_points:
            for band in (FILTER_BANDS if sed_shifted else self.list_bands):
                x_point = LBDA_INFO[band]
                y_point = self.kcorr_data[band][y_plot] if sed_shifted else self.meas_data[band][y_plot]
                y_err_point = self.kcorr_data[band][y_plot+".err"] if sed_shifted else self.meas_data[band][y_plot+".err"]
                ax.errorbar(x_point, y_point, yerr=y_err_point, ls='', marker='o', color=COLOR_INFO[band], label=band)
        
        ax.set_xlabel(r"$\lambda$ [\AA]", fontsize="large")
        ax.set_ylabel(("${m}_{AB}$" if y_plot=="mag" else r"${f}_{\nu}$ $[erg.{s}^{-1}.{cm}^{-2}.{Hz}^{-1}]$"), fontsize="large")
        ax.legend(loc="upper right", ncol=1)
        
        plt.show()
        
        if save_fig:
            save_dir = "~/Desktop/" if save_dir is None else save_dir
            directory_out = os.path.expanduser(save_dir)
            filename_out = "/sed_"+y_plot+("_shifted" if sed_shifted else "")+("_bandpasses" if plot_bandpasses else "")+("_filterpoints" if plot_filter_points else "")+".pdf"
            fig.savefig(directory_out + filename_out)
        
        
    
    
    
    
    
    
    # ================ #
    #  Properties      #
    # ================ #
    @property
    def sed_data(self):
        """ DataFrame of the the SED fit data """
        return self._properties["sed_data"]
    
    @property
    def filter_bandpass(self):
        """ Dictionnary of DataFrame of the bandpass data for every filter band """
        if self._properties["filter_bandpass"] is None:
            self.set_filter_bandpass()
        return self._properties["filter_bandpass"]
    
    @property
    def meas_data(self):
        """ Table like of the measurements """
        if self._properties["meas_data"] is None:
            self._properties["meas_data"] = {}
        return self._properties["meas_data"]
    
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
    def sed_shifted(self):
        """ Dataframe of the shifted to redshift zero SED """
        if self._derived_properties["sed_shifted"] is None:
            self.shift_sed()
        return self._derived_properties["sed_shifted"]
    
    @property
    def kcorr_data(self):
        """ Table of k corrected data """
        if self._derived_properties["kcorr_data"] is None:
            self.k_correction()
        return self._derived_properties["kcorr_data"]
    
    


