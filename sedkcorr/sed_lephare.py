
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
import os
import pkg_resources

from astropy import units
from sncosmo import bandpasses

from pyifu import spectroscopy
from propobject import BaseObject



# ----------------- #
#                   #
#  Functions        #
#                   #
# ------------------#
def lbda_z0(lbda, z):
    return lbda/(1+z)

def flux_z0(flux, z):
    return flux*(1+z)**2

def mag_to_flux(mag):
    return np.power(10, (mag + 48.585)/(-2.5))

def mag_to_flux_err(mag, mag_err, mag0_err=0.005):
    return (mag_to_flux(mag)*np.log(10)/(-2.5))*np.sqrt(mag_err**2 + mag0_err**2)






class SED( BaseObject ):
    """
    LePhare SED fits
    """
    
    PROPERTIES         = ["sed_data", "filter_bandpass", "meas_data", "z"]
    SIDE_PROPERTIES    = ["filter_bandpass_path"]
    DERIVED_PROPERTIES = ["sed_shifted", "kcorr_data"]
    
    FILTER_BANDS = ["FUV", "NUV", "u", "g", "r", "i", "z"]

    LBDA_INFO = {"FUV":1539, "NUV":2316, "u":3562, "g":4719, "r":6185, "i":7500, "z":8961}

    COLOR_INFO = {"FUV":"xkcd:purple", "NUV":"xkcd:violet", "u":"xkcd:blue", "g":"xkcd:green", "r":"xkcd:red", "i":"xkcd:cyan", "z":"xkcd:magenta"}
    
    def set_sed_data(self, sed_filename=None, sed_data=None, nrows=1057):
        """
        
        """
        if sed_filename is not None and sed_data is None:
            self._properties["sed_data"] = pd.read_table(sed_filename, skiprows=20, names=["lbda", "mag"], sep="  ", engine="python", nrows=nrows)
        elif sed_filename is None and sed_data is not None:
            self._properties["sed_data"] = sed_pandas
        else:
            raise ValueError("You must input one (and only one) data set option : the filename directory or a DataFrame with columns = ['lbda', 'mag'].")
        self.sed_data["flux"] = mag_to_flux(self.sed_data["mag"])
    
    def set_filter_bandpass_path(self):
        """
        
        """
        self._side_properties["filter_bandpass_path"] = {band:pkg_resources.resource_filename(__name__, "filter_bandpass/SLOAN_SDSS."+band+".dat") 
                                                         for band in ["u", "g", "r", "i", "z"]}
        self.filter_bandpass_path["FUV"] = pkg_resources.resource_filename(__name__, "filter_bandpass/GALEX_GALEX.FUV.dat")
        self.filter_bandpass_path["NUV"] = pkg_resources.resource_filename(__name__, "filter_bandpass/GALEX_GALEX.NUV.dat")
    
    def set_filter_bandpass(self, opt_bands=None, from_sncosmo=False):
        """
        
        """
        opt_bands = self.FILTER_BANDS if opt_bands is None else opt_bands if type(opt_bands)==list else [opt_bands]
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
    
    def set_meas_data(self, meas_data=None, z=None, col_syntax=["mag_band", "mag_band_err"]):
        """
        
        """
        z = z if z is not None else meas_data["Z-SPEC"]
        col_syntax = ["band", "band.err"] if col_syntax is None else col_syntax
        
        self._properties["meas_data"] = {band:{"mag":meas_data[col_syntax[0].replace("band",band)], 
                                               "mag.err":meas_data[col_syntax[1].replace("band",band)], 
                                               "flux":mag_to_flux(meas_data[col_syntax[0].replace("band",band)]), 
                                               "flux.err":mag_to_flux_err(meas_data[col_syntax[0].replace("band",band)], 
                                                                          meas_data[col_syntax[1].replace("band",band)])}
                                         for band in self.FILTER_BANDS}
        
        self._properties["z"] = z
    
    def shift_sed(self):
        """
        
        """
        self._derived_properties["sed_shifted"] = pd.DataFrame({"lbda":lbda_z0(self.sed_data["lbda"], self.z), 
                                                                "flux":flux_z0(self.sed_data["flux"], self.z)})
    
    def k_correction(self):
        """
        
        """
        self._derived_properties["kcorr_data"] = {band:{"mag.err":self.meas_data[band]["mag.err"], 
                                                        "flux":spectroscopy.synthesize_photometry(self.sed_shifted["lbda"], self.sed_shifted["flux"], 
                                                                                                  self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"]), 
                                                        "flux.err":self.meas_data[band]["flux.err"]}
                                                  for band in self.FILTER_BANDS}
        for band in self.FILTER_BANDS:
            self.kcorr_data[band]["mag"] = -9999999. #TO BE DONE
    
    def show(self, y_plot="flux", sed_shifted=True, plot_bandpasses=False, plot_filter_points=True, xlim=(None, None), ylim=(None, None), save_fig=False, save_dir=None):
        """
        
        """
        x_sed = self.sed_shifted["lbda"] if sed_shifted else self.sed_data["lbda"]
        y_sed = self.sed_shifted[y_plot] if sed_shifted else self.sed_data[y_plot]
        
        fig, ax = plt.subplots()
        opt_sed = {"ls":"-", "marker":"", "color":"0.4"}
        ax.plot(x_sed, y_sed, **opt_sed)
        
        ax.set_xlim(xlim)
        ax.set_xlim(ylim)
        ax_ylim = ax.get_ylim()
        if y_plot=="flux":
            ax_ylim = (0., ax_ylim[1])
            
        ax.axhline(ax_ylim[0], color="black", zorder=5)
            
        if plot_bandpasses:
            for band in self.FILTER_BANDS:
                ax.plot(self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"]*(ax_ylim[1]-ax_ylim[0]) + ax_ylim[0], 
                        ls='-', marker='', color=self.COLOR_INFO[band])
        
        if plot_filter_points:
            for band in self.FILTER_BANDS:
                x_point = self.LBDA_INFO[band]
                y_point = self.kcorr_data[band][y_plot] if sed_shifted else self.meas_data[band][y_plot]
                y_err_point = self.kcorr_data[band][y_plot+".err"] if sed_shifted else self.meas_data[band][y_plot+".err"]
                ax.errorbar(x_point, y_point, yerr=y_err_point, ls='', marker='o', color=self.COLOR_INFO[band])
        
        ax.set_xlabel(r"$\lambda$ [\AA]", fontsize="large")
        ax.set_ylabel(("${m}_{AB}$" if y_plot=="mag" else r"${f}_{\nu}$ $[erg.{s}^{-1}.{cm}^{-2}.{Hz}^{-1}]$"), fontsize="large")
        
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
        """  """
        return self._properties["sed_data"]
    
    @property
    def filter_bandpass(self):
        """  """
        if self._properties["filter_bandpass"] is None:
            self.set_filter_bandpass()
        return self._properties["filter_bandpass"]
    
    @property
    def meas_data(self):
        """  """
        if self._properties["meas_data"] is None:
            self._properties["meas_data"] = {}
        return self._properties["meas_data"]
    
    @property
    def z(self):
        """  """
        return self._properties["z"]
    
    @property
    def filter_bandpass_path(self):
        """  """
        if self._side_properties["filter_bandpass_path"] is None:
            self.set_filter_bandpass_path()
        return self._side_properties["filter_bandpass_path"]
    
    @property
    def sed_shifted(self):
        """  """
        if self._derived_properties["sed_shifted"] is None:
            self.shift_sed()
        return self._derived_properties["sed_shifted"]
    
    @property
    def kcorr_data(self):
        """  """
        if self._derived_properties["kcorr_data"] is None:
            self.k_correction()
        return self._derived_properties["kcorr_data"]
    
    


