
import numpy as np
import pandas
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
import pkg_resources

from astropy import units
from sncosmo import bandpasses

from pyifu import spectroscopy
from propobject import BaseObject


LIST_BANDS = ["FUV", "NUV", "u", "g", "r", "i", "z", "y"]

0.052,0.026,0.05,0.02,0.02,0.02,0.03

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



class KCorrection( BaseObject ):
    """
    This class redshifts an SED spectrum and recover the k-corrected photometry.
    """
    
    PROPERTIES         = ["data_sed", "z", "data_phot"]
    SIDE_PROPERTIES    = ["list_bands", "filter_bandpass"]
    DERIVED_PROPERTIES = ["data_sed_kcorr", "data_phot_kcorr"]
    
    def __init__(self, **kwargs):
        """
        The constructor can automatically call 'set_data'.
        
        Options
        -------
        data : [dict or pandas.DataFrame]
            SED data, containing a Monte Carlo set of spectra.
            You can also give only one SED spectrum, but no error will be calculated on k-corrections.
        
        unit : [string]
            Define 'data' input unit :
            - "mag" : magnitudes
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        z : [float]
            Redshift.
        
        filters : [list(string) or None]
            List of band filters you want to get k-corrected photometry.
            If None, you must execute yourself 'set_list_bands'.
        
        
        Returns
        -------
        Void
        """
        if kwargs != {}:
            self.set_data(**kwargs)
    
    def set_data(self, data=None, unit="mag", z=0, filters=None, **kwargs):
        """
        Set up SED data and redshift attributes.
        
        Parameters
        ----------
        data : [dict or pandas.DataFrame]
            SED data, containing a Monte Carlo set of spectra.
            You can also give only one SED spectrum, but no error will be calculated on k-corrections.
            
        unit : [string]
            Define 'data' input unit :
            - "mag" : magnitudes
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        z : [float]
            Redshift.
        
        Options
        -------
        filters : [list(string) or None]
            List of band filters you want to get k-corrected photometry.
            If None, you must execute yourself 'set_list_bands'.
        
        
        Returns
        -------
        Void
        """
        self._properties["z"] = z
        self.set_list_bands(filters=filters)
        
        input_type = "mag" if unit == "mag" else "flux"
        self._properties["data_sed"] = self._get_dataframe_(data)
        self.data_sed.rename(columns={k:input_type+"_id_{}".format(ii) for ii, k in enumerate(self.data_sed.keys()) if k!="lbda"}, inplace=True)
        if unit == "mag":
            for k in self.data_sed.keys():
                if k == "lbda":
                    continue
                self.data_sed[k], _ = self.mag_to_flux(mag=self.data_sed[k], mag_err=0., band=self.data_sed["lbda"], flux_unit="Hz", opt_mAB0=False)
            self.data_sed.rename(columns={k:k.replace("mag", "flux") for k in self.data_sed.keys() if k != "lbda"}, inplace=True)
        elif unit in ["AA", "mgy"]:
            for k in self.data_sed.keys():
                if k == "lbda":
                    continue
                self.data_sed[k] = self.convert_flux_unit(flux=self.data_sed[k], lbda=self.data_sed["lbda"], unit_in=unit, unit_out="Hz")

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
    
    def set_list_bands(self, filters=["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"]):
        """
        Set as an attribute the given list of the band.
        
        Parameters
        ----------
        filters : [list(string)]
            List of the filters to get k-correction.
        
        
        Returns
        -------
        Void
        """
        self._side_properties["list_bands"] = filters
    
    def set_filter_bandpass(self, filters=None, from_sncosmo=False):
        """
        Load filter bandpasses data.
        
        Options
        -------
        filters : [list[string] or None]
            Optionnal list of bands to load filter bandpass data.
        
        from_sncosmo : [bool]
            If True, filter bandpass data will be loaded from sncosmo package.
            Else, loaded from tables included in the packages, downloaded from "http://svo2.cab.inta-csic.es/svo/theory/fps3/"
            Default is False.
        
        
        Returns
        -------
        Void
        """
        filters = FILTER_BANDS.keys() if filters is None else filters if type(filters)==list else [filters]
        list_sdss_bands = ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"]
        if from_sncosmo:
            self._side_properties["filter_bandpass"] = {_filt:pandas.DataFrame({"lbda":bandpasses.get_bandpass("".join(_filt.split("."))).wave,
                                                                                "trans":bandpasses.get_bandpass("".join(_filt.split("."))).trans})
                                                        for _filt in filters if _filt in list_sdss_bands}
            for _filt in filters:
                if _filt not in list_sdss_bands:
                    self.filter_bandpass[_filt] = pandas.read_csv(self.bandpass_files[_filt], sep=" ", names=["lbda", "trans"])
        else:
            self._side_properties["filter_bandpass"] = {_filt:pandas.read_csv(self.bandpass_files[_filt], sep=" ", names=["lbda", "trans"])
                                                        for _filt in filters}

    def shift_sed(self):
        """
        Shift the SED spectrum to redshift zero.
        
        
        Returns
        -------
        Void
        """
        data = self.data_sed.copy()
        for k, v in data.items():
            if k == "lbda":
                data[k] = self.lbda_z0(v, self.z)
            else:
                data[k] = self.flux_z0(v, self.z)
        self._derived_properties["data_sed_kcorr"] = data
            
    def _get_photopoint_(self, data_sed=None, bands=None):
        """
        Return the integrated photometry given a spectrum in the given bands.
        
        Parameters
        ----------
        data_sed : [pandas.DataFrame]
            Data table of the SED.
        
        Options
        -------
        bands : [string or list(string) or None]
            List of bands you want to get photometric data.
            If one band is given, the function return the value, instead of a dictionary with one key.
        
        
        Returns
        -------
        dict
        """
        list_bands = LIST_BANDS if bands is None else bands if type(bands)==list else [bands]
        data_kcorr = {band:spectroscopy.synthesize_photometry(data_sed["lbda"], data_sed["flux"], self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"])
                      for band in list_bands}
        return data_kcorr if type(bands)!=str else data_kcorr[bands]
    
    def _get_data_phot_(self, data_sed):
        """
        Return the integrated flux from every filter bands from a given SED.
        
        Parameters
        ----------
        data_sed : [pandas.DataFrame]
            SED data to get the photometry.
        
        
        Returns
        -------
        dict
        """
        data_phot = {}
        phot = {}
        for _band in self.list_bands:
            phot["flux"] = [self._get_photopoint_(data_sed=pandas.DataFrame({"lbda":data_sed["lbda"], "flux":v}), bands=_band)
                            for k, v in data_sed.items() if k != "lbda"]
            phot["mag"], _ = self.flux_to_mag(flux=phot["flux"], flux_err=0., band=data_sed["lbda"], flux_unit="Hz", opt_mAB0=False)
            data_phot[_band] = {}
            for k in ["flux", "mag"]:
                quants = np.quantile(phot[k], [0.16, 0.5, 0.84]) if len(phot[k])>1 else phot[k]*3
                data_phot[_band][k] = quants[1]
                data_phot[_band][k+".err_low"] = quants[1] - quants[0]
                data_phot[_band][k+".err_up"] = quants[2] - quants[1]
        return data_phot

    def k_correction(self):
        """
        Apply k-correction on the input SED.
        Firstly, measure the photometry from the input SED.
        Then shift the SED to redshift zero and measure the new (k-corrected) photometry.
        
        
        Returns
        -------
        Void
        """
        self._properties["data_phot"] = self._get_data_phot_(data_sed=self.data_sed)
        self.shift_sed()
        self._derived_properties["data_phot_kcorr"] = self._get_data_phot_(data_sed=self.data_sed_kcorr)

    def _get_fit_quantiles_(self, data_sed, quants=[0.16, 0.5, 0.84], y_unit="AA"):
        """
        Return the given quantiles on the Monte Carlo fitted spectra.
        
        Parameters
        ----------
        data_sed : [pandas.DataFrame]
            SED data to get the quantiles
        
        quants : [list(float)]
            Quantiles you want to get on the spectra. Each quantile must be between 0 and 1.
        
        y_unit : [string]
            Unit of the output spectrum quantiles:
            - "mag" : magnitudes
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        
        Returns
        -------
        Void
        """
        data_sed = data_sed.copy()
        lbda = data_sed["lbda"]
        data_sed.set_index("lbda", inplace=True)
        flux = np.quantile(data_sed, quants, axis=1)
        return {k:self.convert_flux_unit(_f, lbda=lbda, unit_in="Hz", unit_out=y_unit) if y_unit in ["Hz", "AA", "mgy"] else
                  self.flux_to_mag(_f, 0., band=lbda, flux_unit="Hz", opt_mAB0=False)[0]
                for k, _f in zip(quants, flux)}
    
    def show(self, ax=None, y_unit="Hz", sed_shifted=True, plot_bandpasses=False, plot_photo_points=True,
             xlim=(None, None), ylim=(None, None), xscale="linear", yscale="linear", show_sigmas=[1, 2], savefile=None):
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
            
        show_sigmas : [int or list(int) or None]
            Show 1 and/or 2 (or None) sigmas error of the SED.
        
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
        
        data_sed = self.data_sed_kcorr if sed_shifted else self.data_sed
        data_phot = self.data_phot_kcorr if sed_shifted else self.data_phot
        y_plot = "flux" if y_unit in ["Hz", "AA", "mgy"] else "mag"
        mask = np.array([True] * len(data_sed))
        if xlim != (None, None):
            mask = (xlim[0] < np.asarray(data_sed["lbda"])) if xlim[0] is not None else 1
            mask *= (np.asarray(data_sed["lbda"]) < xlim[1]) if xlim[1] is not None else 1
        
        #SED
        x_sed = np.asarray(data_sed["lbda"])[mask]
        y_sed = np.asarray(self._get_fit_quantiles_(data_sed, quants=[0.5], y_unit=y_unit)[0.5])[mask]
        opt_sed = {"ls":"-", "marker":"", "color":"0.4"}
        ax.plot(x_sed, y_sed, label="_nolegend_", **opt_sed)

        nsigmas = len(np.atleast_1d(show_sigmas))
        if 2 in show_sigmas:
            ff = self._get_fit_quantiles_(data_sed, quants=[0.05, 0.95], y_unit=y_unit)
            ax.fill_between(x_sed, np.asarray(ff[0.05])[mask], np.asarray(ff[0.95])[mask], alpha=0.3/nsigmas, color="C0", lw=0, zorder=1)
        if 1 in show_sigmas:
            ff = self._get_fit_quantiles_(data_sed, quants=[0.16, 0.84], y_unit=y_unit)
            ax.fill_between(x_sed, np.asarray(ff[0.16])[mask], np.asarray(ff[0.84])[mask], alpha=0.3/nsigmas, color="C0", lw=0, zorder=2)
        
        #Photopoints
        if plot_photo_points:
            for _band in self.list_bands:
                x_point = FILTER_BANDS[_band]["lbda"]
                y_point = data_phot[_band][y_plot]
                y_point_err_low = data_phot[_band][y_plot+".err_low"]
                y_point_err_up = data_phot[_band][y_plot+".err_up"]
                if y_plot == "flux":
                    y_point, y_point_err_low, y_point_err_up = self.convert_flux_unit((y_point, y_point_err_low, y_point_err_up),
                                                                                      lbda=FILTER_BANDS[_band]["lbda"], unit_in="Hz", unit_out=y_unit)
                ax.errorbar(x_point, y_point, yerr=[[y_point_err_low], [y_point_err_up]],
                            ls="", marker="o", color=FILTER_BANDS[_band]["color"], label=_band)
        
        ax.set_xlabel(r"$\lambda$ [\AA]", fontsize="large")
        ylabel = "mag"
        if y_plot == "flux":
            ylabel = "mgy" if y_unit == "mgy" else \
                     r"${{f}}_{{\nu}}$ $[erg.{{s}}^{{-1}}.{{cm}}^{{-2}}.{Hz}^{{-1}}]$" if y_unit == "Hz" else \
                     r"${{f}}_{{\lambda}}$ $[erg.{{s}}^{{-1}}.{{cm}}^{{-2}}.{\AA}^{{-1}}]$"
        ax.set_ylabel(ylabel, fontsize="large")
                
        ax.set_xlim(xlim)
        if ylim == (None, None):
            if ax.get_ylim()[0] <= 0:
                ylim = (np.min(y_sed), ylim[1])
        ax.set_ylim(ylim)
#        ax_ylim = ax.get_ylim()
#        if y_plot == "flux" and ylim[0] is None:
#            ax_ylim = (0., ax_ylim[1])

#        ax.axhline(ax_ylim[0], color="black", lw=0.5, zorder=5)

        #Photometry
        if plot_bandpasses:
            for band in self.list_bands:
                ax.plot(self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"]*(ax.get_ylim()[1]-ax.get_ylim()[0]) + ax.get_ylim()[0],
                        ls="-", marker="", color=FILTER_BANDS[band]["color"], label=("_nolegend_" if plot_photo_points else band))

        if plot_bandpasses or plot_photo_points:
            ax.legend(loc="upper right", ncol=1)

        ax.set_xscale(xscale)
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
        zp = KCorrection.convert_flux_unit(flux, lbda, "Hz", "AA")
        return np.log10(zp)
    
    @staticmethod
    def mag_to_flux(mag, mag_err=0., band=None, flux_unit="Hz", opt_mAB0=False):
        """
        Convert magnitude to flux.
        Return the flux and its error.
        
        Parameters
        ----------
        mag : [float or np.array]
            Magnitude.
        
        mag_err : [float or np.array]
            Magnitude error.
        
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
        
        
        Returns
        -------
        float or np.array, float or np.array
        """
        if type(mag_err) == float and mag_err == 0. and len(np.atleast_1d(mag)) > 1:
            mag_err = np.zeros(len(mag))
        if opt_mAB0:
            flux_out = 10**(KCorrection.mag_ab_zeropoint(FILTER_BANDS[band]["ZP"], FILTER_BANDS[band]["lbda"]) - 0.4*mag)
            flux_err_out = (0.4 * np.log(10) * flux_out) * mag_err**2
            unit_in = "AA"
        else:
            flux_out = 10**((mag + 48.585)/(-2.5))
            flux_err_out = (0.4 * np.log(10) * flux_out)*np.sqrt(mag_err**2 + 0.005**2)
            unit_in = "Hz"
        
        if flux_unit in ("AA", "Hz", "mgy"):
            flux_out, flux_err_out = KCorrection.convert_flux_unit((flux_out, flux_err_out), FILTER_BANDS[band]["lbda"] if type(band)==str else band,
                                                                   unit_in=unit_in, unit_out=flux_unit)
        else:
            raise ValueError("{} is not a valid flux unit.".format(flux_unit))
        
        return flux_out, flux_err_out
    
    @staticmethod
    def flux_to_mag(flux, flux_err=0., band=None, flux_unit="Hz", opt_mAB0=False):
        """
        Convert flux to magnitude.
        Return the magnitude and its error.
        
        Parameters
        ----------
        flux : [float or np.array]
            Flux.
        
        flux_err : [float or np.array]
            Flux error.
        
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
        
        
        Returns
        -------
        float or np.array, float or np.array
        """
        if type(flux_err) == float and flux_err == 0. and len(np.atleast_1d(flux)) > 1:
            flux_err = np.zeros(len(flux))
        if opt_mAB0:
            unit_out = "AA"
        else:
            unit_out = "Hz"
        if flux_unit in ("AA", "Hz", "mgy"):
            flux, flux_err = KCorrection.convert_flux_unit((flux, flux_err), FILTER_BANDS[band]["lbda"] if type(band)==str else band,
                                                           unit_in=flux_unit, unit_out=unit_out)
        else:
            raise ValueError("{} is not a valid flux unit.".format(flux_unit))

        if opt_mAB0:
            mag_out = -2.5 * (np.log10(flux) - KCorrection.mag_ab_zeropoint(FILTER_BANDS[band]["ZP"], FILTER_BANDS[band]["lbda"]))
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
    
    
    


    # ================ #
    #  Properties      #
    # ================ #
    @property
    def data_sed(self):
        """ DataFrame of the the SED fit data """
        return self._properties["data_sed"]
    
    @property
    def z(self):
        """ Host redshift """
        return self._properties["z"]
    
    @property
    def filter_bandpass(self):
        """ Dictionnary of DataFrame of the bandpass data for every filter band """
        if self._side_properties["filter_bandpass"] is None:
            self.set_filter_bandpass()
        return self._side_properties["filter_bandpass"]
    
    @property
    def data_phot(self):
        """ Table of photopoints of the input SED """
        return self._properties["data_phot"]
    
    @property
    def bandpass_files(self):
        """ Dictionnary of paths for every filter badpass data """
        return {_filt:pkg_resources.resource_filename(__name__, _dict["bandpass_file"]) for _filt, _dict in FILTER_BANDS.items()}
    
    @property
    def list_bands(self):
        """ List of the filter bands used in SED fitting """
        if self._side_properties["list_bands"] is None:
            raise AttributeError("'list_bands' is None, you must execute 'set_list_bands' or give a value to 'filters' parameter in 'set_dat' method.")
        return self._side_properties["list_bands"]
    
    @property
    def data_sed_kcorr(self):
        """ Dataframe of the shifted to redshift zero SED """
        return self._derived_properties["data_sed_kcorr"]
    
    @property
    def data_phot_kcorr(self):
        """ Table of k corrected photopoints """
        return self._derived_properties["data_phot_kcorr"]
