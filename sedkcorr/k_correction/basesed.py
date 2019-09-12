
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
                       "prospector_name":"galex_FUV",
                       "lephare_name":"galex/FUV.pb"},
                "NUV":{"lbda":2316,
                       "color":"xkcd:violet",
                       "mAB0":np.log10(2.06e-16) + 0.4 * 20.08,
                       "context_id":1,
                       "prospector_name":"galex_NUV",
                       "lephare_name":"galex/NUV.pb"},
                "u":{"lbda":3562,
                     "color":"xkcd:blue",
                     "mAB0":-8.056,
                     "context_id":2,
                     "prospector_name":"sdss_u0",
                     "lephare_name":"sdss/up.pb"},
                "g":{"lbda":4719,
                     "color":"xkcd:green",
                     "mAB0":-8.326,
                     "context_id":3,
                     "prospector_name":"sdss_g0",
                     "lephare_name":"sdss/gp.pb"},
                "r":{"lbda":6185,
                     "color":"xkcd:red",
                     "mAB0":-8.555,
                     "context_id":4,
                     "prospector_name":"sdss_r0",
                     "lephare_name":"sdss/rp.pb"},
                "i":{"lbda":7500,
                     "color":"xkcd:cyan",
                     "mAB0":-8.732,
                     "context_id":5,
                     "prospector_name":"sdss_i0",
                     "lephare_name":"sdss/ip.pb"},
                "z":{"lbda":8961,
                     "color":"xkcd:magenta",
                     "mAB0":-8.882,
                    "context_id":6,
                    "prospector_name":"sdss_z0",
                    "lephare_name":"sdss/zp.pb"},
                }



class KCorrection( BaseObject ):
    """
    This class redshifts an SED spectrum and recover the k-corrected photometry.
    """
    
    PROPERTIES         = ["data_sed", "data_meas"]
    SIDE_PROPERTIES    = ["list_bands", "filter_bandpass"]
    DERIVED_PROPERTIES = ["data_sed_shifted", "data_kcorr"]
    
    def __init__(self, **kwargs):
        """
        The constructor can automatically call 'set_data'.
        
        Options
        -------
        ### set_data_sed ###
        data_sed : [pandas.DataFrame or dict]
            Data table of the SED.
        
        ### set_data_meas ###
        data_meas : [pandas.DataFrame or table]
            Data table of the measurements.
        
        z : [list(float) or table]
            Redshift table.
        
        
        Returns
        -------
        Void
        """
        if kwargs != {}:
            self.set_data(**kwargs)
    
    def set_data(self, **kwargs):
        """
        Execute 'set_data_sed' and 'set_data_meas'.
        Basically set up data attributes.
        
        Parameters
        ----------
        ### set_data_sed ###
        data_sed : [pandas.DataFrame or dict]
            Data table of the SED.
        
        ### set_data_meas ###
        data_meas : [pandas.DataFrame or table]
            Data table of the measurements.
        
        z : [list(float) or table]
            Redshift table.
        
        
        Returns
        -------
        Void
        """
        self.set_data_sed(**kwargs)
        self.set_data_meas(**kwargs)
            
    def set_data_sed(self, data_sed=None, flux_unit="Hz", **extras):
        """
        Turn the input SED data into DataFrame, unless already in pandas.DataFrame type.
        Automatically add either flux or magnitudes columns if one is missing.
        
        Parameters
        ----------
        data_sed : [pandas.DataFrame or dict or table like]
            Data table of the SED.
        
        flux_unit : [string]
            If 'data_sed' given in flux, define here the flux unit :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        
        Returns
        -------
        Void
        """
        self._properties["data_sed"] = self._get_dataframe_(data_sed)
        input_type = "mag" in self.data_sed.keys() else "flux"
        conv_type = "flux" if input_type=="mag" else "mag"
        f_conv = self.mag_to_flux if input_type=="mag" else self.flux_to_mag
        
        if input_type+".err_low" not in self.data_sed.keys() and input_type+".err" in self.data_sed.keys():
            self.data_sed[input_type+".err_low"] = self.data_sed[input_type+".err"]
            self.data_sed[input_type+".err_up"] = self.data_sed[input_type+".err"]
            self.data_sed.drop(columns=[input_type+".err"], inplace=True)
        elif input_type+".err_low" not in self.data_sed.keys() and input_type+".err" not in self.data_sed.keys():
            raise ValueError("'data_sed' needs either '{0}.err' or '{0}.err_low' and '{0}.err_up' columns.")
        
        if input_type == "flux":
            self.data_sed["flux"] = convert_flux_unit(self.data_sed["flux"], lbda=self.data_sed["lbda"], unit_in=flux_unit, unit_out="Hz")
            self.data_sed["flux.err_low"] = convert_flux_unit(self.data_sed["flux.err_low"], lbda=self.data_sed["lbda"], unit_in=flux_unit, unit_out="Hz")
            self.data_sed["flux.err_up"] = convert_flux_unit(self.data_sed["flux.err_up"], lbda=self.data_sed["lbda"], unit_in=flux_unit, unit_out="Hz")
        
        self.data_sed[conv_type] = f_conv(self.data_sed[input_type], np.zeros(len(data_sed)), band=self.data_sed["lbda"], flux_unit="Hz", opt_mAB0=False)
        self.data_sed[conv_type+".err_low"] = self.data_sed[conv_type] - f_conv(self.data_sed[input_type]-self.data_sed[input_type+".err_low"],
                                                                                np.zeros(len(data_sed)), band=self.data_sed["lbda"], flux_unit="Hz", opt_mAB0=False)
        self.data_sed[conv_type+".err_up"] = - self.data_sed[conv_type] + f_conv(self.data_sed[input_type]+self.data_sed[input_type+".err_up"],
                                                                                 np.zeros(len(data_sed)), band=self.data_sed["lbda"], flux_unit="Hz", opt_mAB0=False)

    def set_data_meas(self, data_meas=None, z=None, flux_unit="Hz", **extras):
        """
        Set up measurements attributes.
        
        Parameters
        ----------
        data_meas : [pandas.DataFrame or dict or table like]
            Data table of the measurements.
        
        z : [list(float) or table]
            Redshift data. Only used if no 'z' column name in 'data_meas'.
        
        flux_unit : [string]
            If 'data_sed' given in flux, define here the flux unit :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        
        Returns
        -------
        Void
        """
        self._properties["data_meas"] = self._get_dataframe_(data_meas)
        if "z" not in self.data_meas.keys() and z is not None:
            self.data_meas["z"] = z
        elif "z" not in self.data_meas.keys() and z is None:
            raise AttributeError("You have to input redshift data, either in 'data_meas' (with 'z' column name) or in 'z').")
        
        ######### Option to set : mag or flux ############## !!!!!!!!!!!!!!! To be changed !!!!!!!!!!!!!!!!!
        for band in self.list_bands:
            if "mag" in self.data_meas[band].keys() and "flux" not in self.data_meas[band].keys():
                self.data_meas[band]["flux"], self.data_meas[band]["flux.err"] = mag_to_flux(self.data_meas[band]["mag"], self.data_meas[band]["mag.err"],
                                                                                             band=band, flux_unit="Hz", opt_mAB0=True)
            elif "flux" in self.data_meas[band].keys() and "mag" not in self.data_meas[band].keys():
                self.data_meas[band]["mag"], self.data_meas[band]["mag.err"] = flux_to_mag(self.data_meas[band]["flux"], self.data_meas[band]["flux.err"],
                                                                                           band=band, flux_unit=flux_unit, opt_mAB0=True)
    
    def _get_dataframe_(self, data):
        """
        Return a pandas.DataFrame of the data.
        
        Parameters
        ----------
        data : [pandas.DataFrame or dict or table like]
            Data to transform into pandas.DataFrame.
        
        
        Returns
        -------
        pandas.DataFrame
        """
        if type(data) is pandas.DataFrame:
            data = data
        else:
            try:
                data = pandas.DataFrame(data)
            except:
                raise TypeError("'data' must be transformable into a pandas.DataFrame (eg: pandas.DataFrame, dict, etc.)")
        return data
    
    def _get_full_data_(self, data):
        """
        Return the data containing magnitudes and flux within the unit treated in the class code, and the errors in a good shape considering asymetric errors if there are any.
        
        Parameters
        ----------
        data : [pandas.DataFrame]
            Input data.
        
        
        Returns
        -------
        pandas.DataFrame
        """
        if ("mag" in self.data_sed.keys() and "flux" not in self.data_sed.keys()):
            input_type = "mag"
        elif ("mag" not in self.data_sed.keys() and "flux" in self.data_sed.keys()):
            input_type = "flux"
        elif ("mag" in self.data_sed.keys() and "flux" in self.data_sed.keys()):
            input_type = None
        else:
            raise ValueError("'data' neither contain flux data nor magnitude data.")
        
        if input_type in ["mag", "flux"]:
            conv_type = "flux" if input_type=="mag" else "mag"
            f_conv = self.mag_to_flux if input_type=="mag" else self.flux_to_mag
        return
    
    def _set_list_bands_(self):
        """
        Set as an attribute the list of the bands contained in the 'data_meas' photometry.
        
        
        Returns
        -------
        Void
        """
        list_bands
        self._side_properties["list_bands"]
        return
    
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
        self._derived_properties["data_sed_shifted"] = pandas.DataFrame({"lbda":lbda_z0(self.data_sed["lbda"], self.data_meas["z"]),
                                                                         "flux":flux_z0(self.data_sed["flux"], self.data_meas["z"]),
                                                                         "flux.err":flux_z0(self.data_sed["flux.err"], self.data_meas["z"])})
        self.data_sed_shifted["mag"], self.data_sed_shifted["mag.err"] = flux_to_mag(self.data_sed_shifted["flux"], self.data_sed_shifted["flux.err"], band=None, flux_unit="Hz")
    
    def k_correction(self, kcorr_flux_error=None):
        """
        Recover the integrated flux from every filter bands from the shifted SED.
        Then convert them into magnitudes.
        
        Options
        -------
        kcorr_flux_error : [dict or table or None]
            Error on k-corrected flux.
        
        
        Returns
        -------
        Void
        """
        data_kcorr = self.get_phot(data_sed=self.data_sed_shifted, bands=None)
        data_kcorr = {k:{"flux":v} for k, v in data_kcorr.items()}
        self._derived_properties["data_kcorr"] = data_kcorr
        for band in LIST_BANDS:
            self.data_kcorr[band]["flux.err"] = 0. if kcorr_flux_error is None else kcorr_flux_error[band]
            self.data_kcorr[band]["mag"], self.data_kcorr[band]["mag.err"] = flux_to_mag(self.data_kcorr[band]["flux"], self.data_kcorr[band]["flux.err"], band=band, flux_unit="Hz")
    
    def show(self, ax=None, y_unit="Hz", sed_shifted=True, plot_bandpasses=False, plot_photo_points=True,
             xlim=(None, None), ylim=(None, None), xscale="linear", yscale="linear", savefile=None):
        """
        Plot method.
        
        Parameters
        ----------
        y_unit : [string]
            Choice to plot "mag" or flux with :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        sed_shifted : [bool]
            If True, the SED is shifted to redshift zero.
            Else, the LePhare fitted SED is plotted.
        
        Options
        -------
        ax : [matplotlib.axes]
            Already existing axes you want to add stuff in.
            Else, None.
        
        plot_bandpasses : [bool]
            If True, plot the filter bandpass transmitions.
        
        plot_photo_points : [bool]
            If True, plot the filter points with errors, either in flux or in magnitude.
        
        xlim : [tuple[float or None]]
            Set the limits on the x axis.
            If (None, None), the figure has free x axis limits.
        
        ylim : [tuple[float or None]]
            Set the limits on the y axis.
            If (None, None), the figure has free y axis limits.
        
        xscale : [string]
            Scale of the x axis : "linear", "log", ...
            
        yscale : [string]
            Scale of the y axis : "linear", "log", ...
        
        savefile : [string or None]
            If None, the figure won't be saved.
            To save it, input a path directory + filename.
        
        Returns
        -------
        dict(fig, ax)
        """
        if ax is not None:
            fig = ax.figure
        else:
            fig, ax = plt.subplots()
        
        data_sed = self.data_sed_shifted if sed_shifted else self.data_sed
        data_phot = self.data_kcorr if sed_shifted else self.data_meas
        y_plot = "flux" if y_plot in ["Hz", "AA", "mgy"] else "mag"
        
        #SED
        x_sed = data_sed["lbda"]
        y_sed = data_sed[y_plot]
        y_sed_err_low = data_sed[y_plot+".err_low"] if y_plot+".err_low" in data_sed.keys() else data_sed[y_plot+".err"]
        y_sed_err_up = data_sed[y_plot+".err_up"] if y_plot+".err_up" in data_sed.keys() else data_sed[y_plot+".err"]
        if y_plot in ["Hz", "AA", "mgy"]:
            y_sed, y_sed_err_low, y_sed_err_up = convert_flux_unit((y_sed, y_sed_err_low, y_sed_err_up),
                                                                   lbda=x_sed, unit_in="Hz", unit_out=flux_unit)
        
        opt_sed = {"ls":"-", "marker":"", "color":"0.4"}
        ax.plot(x_sed, y_sed, label="_nolegend_", **opt_sed)
        ax.fill_between(x_sed, y_sed - y_sed_err_low, y_sed + y_sed_err_up, alpha=0.5, color="0.7")
        
        #Photopoints
        if plot_photo_points:
            for band in (LIST_BANDS if sed_shifted else self.list_bands):
                x_point = FILTER_BANDS[band]["lbda"]
                y_point = data_phot[band][y_plot]
                y_point_err_low = data_phot[y_plot+".err_low"] if y_plot+".err_low" in data_phot.keys() else data_phot[y_plot+".err"]
                y_point_err_up = data_phot[y_plot+".err_up"] if y_plot+".err_up" in data_phot.keys() else data_phot[y_plot+".err"]
                if y_plot == "flux":
                    y_point, y_point_err_low, y_point_err_up = convert_flux_unit((y_point, y_point_err_low, y_point_err_up),
                                                                                 lbda=FILTER_BANDS[band]["lbda"], unit_in="Hz", unit_out=flux_unit)
                ax.errorbar(x_point, y_point, yerr=[[y_point_err_low], [y_point_err_up]],
                            ls="", marker="o", color=FILTER_BANDS[band]["color"], label=band)
        
        ax.set_xlabel(r"$\lambda$ [\AA]", fontsize="large")
        if y_plot == "flux":
            ylabel = "mgy" if flux_unit == "mgy" else \
                     r"${{f}}_{{\nu}}$ $[erg.{{s}}^{{-1}}.{{cm}}^{{-2}}.{Hz}^{{-1}}]$" if flux_unit == "Hz" else \
                     r"${{f}}_{{\lambda}}$ $[erg.{{s}}^{{-1}}.{{cm}}^{{-2}}.{\AA}^{{-1}}]$"
        ax.set_ylabel(ylabel, fontsize="large")
        ax.legend(loc="upper right", ncol=1)
                
        ax.set_xlim(xlim)
        if ylim == (None, None):
            xmin, xmax = ax.get_xlim()
            mask = (xmin < np.asarray(x_sed)) * (np.asarray(x_sed) < xmax)
            ymin, ymax = np.min((np.asarray(y_sed)-np.asarray(y_sed_err))[mask]), np.max((np.asarray(y_sed)+np.asarray(y_sed_err))[mask])
            ylim = (ymin if ymin>0 else 0, ymax)
        ax.set_ylim(ylim)
        ax_ylim = ax.get_ylim()
        if y_plot=="flux" and ylim[0] is None:
            ax_ylim = (0., ax_ylim[1])

        ax.axhline(ax_ylim[0], color="black", zorder=5)
        
        #Photometry
        if plot_bandpasses:
            for band in LIST_BANDS:
                ax.plot(self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"]*(ax_ylim[1]-ax_ylim[0]) + ax_ylim[0],
                        ls="-", marker="", color=FILTER_BANDS[band]["color"], label="_nolegend_")

        ax.set_xscale(xscale)
        if y_plot=="flux":
            ax.set_yscale(yscale)
                
        if savefile is not None:
            fig.savefig(savefile)
                        
        return {"ax":ax, "fig":fig}
    
    
    
    
    # ================ #
    #  Static methods  #
    # ================ #
    @staticmethod
    def lbda_z0(lbda, z):
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
    
    @staticmethod
    def flux_z0(flux, z):
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
    
    @staticmethod
    def mag_to_flux(mag, mag_err=0., band=None, flux_unit="Hz", opt_mAB0=True):
        """
        Convert magnitude to flux.
        Return the flux and its error.
        
        Parameters
        ----------
        mag : [float or np.array]
            Magnitude.
        
        mag_err : [float or np.array]
            Magnitude error.
        
        band : [string or None]
            The output flux is converted in the input unit with the wavelength based on this given 'band' filter.
        
        flux_unit : [string]
            Define the output flux unit :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        Options
        -------
        opt_mAB0 : [bool]
            If True, take in account the filter calibration (AB zero magnitude, depend on the instrument (SDSS, Galex, etc.).
        
        
        Returns
        -------
        float or np.array, float or np.array
        """
        if opt_mAB0:
            flux_out = 10**(FILTER_BANDS[band]["mAB0"] - 0.4*mag)
            flux_err_out = (0.4 * np.log(10) * flux_out) * mag_err**2
            unit_in = "AA"
        else:
            flux_out = 10**((mag + 48.585)/(-2.5))
            flux_err_out = (0.4 * np.log(10) * flux_out)*np.sqrt(mag_err**2 + 0.005**2)
            unit_in = "Hz"
        
        if flux_unit in ("AA", "Hz", "mgy"):
            flux_out, flux_err_out = convert_flux_unit((flux_out, flux_err_out), FILTER_BANDS[band]["lbda"] if type(band)==str else band, unit_in=unit_in, unit_out=flux_unit)
        else:
            raise ValueError("{} is not a valid flux unit.".format(flux_unit))
        
        return flux_out, flux_err_out
    
    @staticmethod
    def flux_to_mag(flux, flux_err=0., band=None, flux_unit="Hz", opt_mAB0=True):
        """
        Convert flux to magnitude.
        Return the magnitude and its error.
        
        Parameters
        ----------
        flux : [float or np.array]
            Flux.
        
        flux_err : [float or np.array]
            Flux error.
        
        band : [string or None]
            The input flux is converted in a conversion convenient unit from the input unit with the wavelength based on this given 'band' filter.
        
        flux_unit : [string]
            Define the input flux unit :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        Options
        -------
        opt_mAB0 : [bool]
            If True, take in account the filter calibration (AB zero magnitude), depending on the instrument (SDSS, Galex, etc.).
        
        
        Returns
        -------
        float or np.array, float or np.array
        """
        if opt_mAB0:
            unit_out = "AA"
        else:
            unit_out = "Hz"
        if flux_unit in ("AA", "Hz", "mgy"):
            flux, flux_err = convert_flux_unit((flux, flux_err), FILTER_BANDS[band]["lbda"] if type(band)==str else band, unit_in=flux_unit, unit_out=unit_out)
        else:
            raise ValueError("{} is not a valid flux unit.".format(flux_unit))

        if opt_mAB0:
            mag_out = -2.5 * (np.log10(flux) - FILTER_BANDS[band]["mAB0"])
            mag_err_out = (2.5 / np.log(10)) * (flux_err / flux)
        else:
            mag_out = -2.5 * np.log10(flux) - 48.585
            mag_err_out = np.sqrt(((2.5/np.log(10))*(flux_err/flux))**2 + 0.005**2)
            
        return mag_out, mag_err_out
    
    @staticmethod
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
        flux = np.asarray(flux) if type(flux) != float else flux
        lbda = np.asarray(lbda) if type(lbda) != float else lbda
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
    
    @staticmethod
    def get_phot(self, data_sed=None, bands=None):
        """
        
        
        Parameters
        ----------
        data_sed : [pandas.DataFrame]
            Data table of the SED.
        
        Options
        -------
        bands : [list(string) or None]
            List of bands you want to get photometric data.
        
        
        Returns
        -------
        dict
        """
        list_bands = LIST_BANDS if bands is None else bands if type(bands)==list else [bands]
        data_kcorr = {band:spectroscopy.synthesize_photometry(data_sed["lbda"], data_sed["flux"], self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"])
                      for band in list_bands}
        return data_kcorr

    
    
    


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
