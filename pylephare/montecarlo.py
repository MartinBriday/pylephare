""" Handles the target SEDfitting Monte Carlo """

import pandas
import numpy as np
from . import io, base

class MCLePhare( base._FilterHolder_ ):
    """ """
    def __init__(self, serie=None, ndraw=500):
        """ """
        self.set_data(serie)
        if self.has_data() and (ndraw is not None and ndraw>1):
            self.draw_mc(ndraw)
            
    # ============== #
    #  Methods       #
    # ============== #
    def run(self, configfile = None, **kwargs):
        """ """
        if not self.has_lephare() or configfile is not None:
            self.load_lephare(configfile)
            
        self._lephare_out = self.lephare.run(**kwargs)

    # ------- #
    # GETTER  #
    # ------- #
    def get_spectrum(self, index):
        """ Returns the index-th LePhareSpectrum """
        return self.spectra.spectra[index]
    
    def get_synthetic_photometry(self, filter_, restframe=False, influx=True, inhz=False):
        """ get photometry synthesized trought the given filter/bandpass

        Parameters
        ----------
        filter_: [string, sncosmo.BandPass, 2D array]
            the filter through which the spectrum will be synthesized.
            accepted input format:
            - string: name of a know filter, the actual bandpass will be grab using io.get_filter_bandpass(filter_)
            - 2D array: bandpass = sncosmo.bandpass.BandPass(*filter_)
            - sncosmo.bandpass.BandPass
            
        restframe: [bool] -optional-
            The spectrum is first deredshifted before doing the synthetic photometry

        influx: [bool] -optional-
            Shall the result returned in flux (True) or in ABmag (False)

        inhz: [bool] -optional-
            Should the returned flux be in  erg/s/cm2/A or erg/s/cm2/Hz
            // ignored if influx = False //

        Returns
        -------
        effective wavelength, synthesize flux/mag (see influx)
        """
        return self.spectra.get_synthetic_photometry(filter_, restframe=restframe, influx=influx, inhz=inhz)
    
    # ------- #
    # SETTER  #
    # ------- #
    def set_data(self, serie):
        """ """
        if not type(serie) == pandas.core.series.Series:
            raise TypeError("insut serie must be pandas.Serie.")
        self._data = serie
        filters = io.keys_to_filters(serie.keys())
        self.set_filters(filters) 

    def set_intrinsic_scatter(self, magerr, onconfig=False, redraw=True):
        """ This sets the intrinisic scatter
        This scatter could be added during the lephare fit itself or when drawing the MonteCarlo.

        Parameters
        ----------
        magerr: [1D array]
            same number as magerr as you have filters
            
        onconfig: [bool] -optional-
            - False: Individual data error will be increased by the given magnitude error
                     and the ERR_SCALE config parameter will be set to None.
            - True: Individual data error are unchanged, magerr passed to ERR_SCALE.
            
        redraw: [bool] -optional-
            If redraw the mcdata.
            // Ignored if onconfig is True //
            
        Returns
        -------
        Void

        """
        if onconfig:
            if not self.has_lephare():
                raise AttributeError("lephare is not loaded.")
            self.lephare.config.set_intrinsic_error(magerr)
            self._sigma_int = None
        elif len( magerr) != self.nfilters:
            raise ValueError("You must provide exactly one magerr per filterband")
        
        #dflux = np.abs(flux*(-magerr/2.5*np.log(10)))
        self._sigma_int = {k:v for k,v in zip(self.filters,magerr)}
        
        if redraw:
            if self.ndraw is not None:
                self.draw_mc(self.ndraw)
            else:
                self.draw_mc()
            
    def draw_mc(self, ndraw=500):
        """ """
        data_mc = {}
        for filt_ in self.filters:
            if self.sigma_int is None:
                filt_err = self.data[filt_+".err"]
            else:
                filt_err = np.sqrt(self.data[filt_+".err"]**2  + self.sigma_int[filt_]**2)
            data_mc[filt_] = np.random.normal(loc=self.data[filt_], scale= filt_err, size=ndraw)
            data_mc[filt_+".err"] = filt_err
                
        mcdata = pandas.DataFrame(data_mc)
        for k in self.data.keys():
            if k not in mcdata.columns:
                mcdata[k] = self.data[k]
                
        self._mcdata = mcdata

    # ------- #
    # LOADER  #
    # ------- #
    def load_lephare(self, configfile):
        """ """
        from sedkcorr import lephare
        self._lephare = lephare.LePhare(self.mcdata, configfile)
        self._lephare.config.set_value("CAT_FMT","MMEE")
        
    def _load_results_(self, verbose=True):
        """ """
        if verbose:
            print("loading results...")
            
        from pylephare import configparser, spectrum, lephare
        self._config = configparser.ConfigParser(self._lephare_out["config"])
        self._catin = pandas.read_csv(self._lephare_out["catin"], sep=" ")
        self._catout = lephare.read_catout(self._lephare_out["catout"], self.filters)
        self._spectra = spectrum.LePhareSpectrumCollection.read_files(self._lephare_out["spec"], lbda_range=[3000,10000])
        
    # ------- #
    # PLOTTER #
    # ------- #
    def show(self, ax=None, nmc=10, add_spec=True, inhz=False):
        """ """
        if inhz:
            raise NotImplementedError("inhz plotting option is not implemented yet.")
            
        import matplotlib.pyplot as mpl
        if ax is None:
            fig = mpl.figure(figsize=[7,3.5])
            ax = fig.add_axes([0.1,0.2,0.8,0.7])
        else:
            fig = ax.figure
        
        prop = dict(ls="None", marker="o", ecolor="0.7", mfc="C0", mec="w", ms=10, zorder=5)
        propmc = dict(ls="None", marker="o", mfc="0.7", mec="None", ms=5, alpha=0.3, zorder=3)
        idmc = np.random.choice(np.arange(self.ndraw), replace=False, size=nmc)
        for f in self.filters:
            ax.errorbar(self.filter_bandpasses[f].wave_eff, self.data[f], 
                       yerr=self.data[f+".err"], **prop)
            
            ax.plot([self.filter_bandpasses[f].wave_eff]*nmc, self.mcdata[f][idmc],  **propmc)
            
        if add_spec:
            _ = [self.spectra.spectra[i_].show(ax=ax, showdata=False, showmagmodel=False,
                                        set_label=False, lw=1, alpha=0.3)
                for i_ in idmc]
            
        
        ax.set_xlabel(r"Wavelentgh [$\AA$]", fontsize="large")
        ax.set_ylabel(r"flux [$\mathrm{erg\,s^{-1}\,cm^{-2}\,%s}$]"%("Hz^{-1}" if inhz else "\AA^{-1}"),  fontsize="large")

    # ============== #
    #  Properties    #
    # ============== #
    # - Data
    @property
    def data(self):
        """ DataFrame containing the data """
        if not hasattr(self,"_data"):
            self._data = None

        return self._data

    def has_data(self):
        """ """
        return self.data is not None and len(self.data)>0

    @property
    def sigma_int(self):
        """ Intrinsic magnitude scatter to be applied to the data when drawing """
        if not hasattr(self, "_sigma_int"):
            self._sigma_int = None
        return self._sigma_int
        
    @property
    def mcdata(self):
        """ DataFrame containing the data """
        if not hasattr(self,"_mcdata"):
            self._mcdata = None

        return self._mcdata

    def has_mcdata(self):
        """ """
        return self.mcdata is not None and len(self.mcdata)>0

    @property
    def ndraw(self):
        """ """
        return len(self.mcdata) if self.has_mcdata() else None
        
    # - LePhare
    @property
    def lephare(self):
        """ """
        if not self.has_lephare():
            raise AttributeError("LePhare is not loaded, see self.load_lephare()")
        return self._lephare
    
    def has_lephare(self):
        """ """
        return hasattr(self, "_lephare") and self._lephare is not None

    def _did_run_(self):
        """ test if you ran self.run()"""
        return hasattr(self, "_lephare_out")

    # - lephare.run() results
    @property
    def catout(self):
        """ Output catalog"""
        if not hasattr(self, "_catout"):
            if self._did_run_():
                self._load_results_()
            else:
                raise AttributeError("you did not run self.run()")
                
        return self._catout
    
    @property
    def config(self):
        """ config used to run()"""
        if not hasattr(self, "_config"):
            if self._did_run_():
                self._load_results_()
            else:
                raise AttributeError("you did not run self.run()")
                
        return self._config

    @property
    def spectra(self):
        """ spectra estimated from run()"""
        if not hasattr(self, "_spectra"):
            if self._did_run_():
                self._load_results_()
            else:
                raise AttributeError("you did not run self.run()")
                
        return self._spectra
