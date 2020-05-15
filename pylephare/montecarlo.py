""" Handles the target SEDfitting Monte Carlo """

import pandas
import warnings
import numpy as np
from scipy import stats
from . import io, base

class MCLePhare( base._FilterHolder_ ):
    """ Class handling LePhare SED fitting coupled with a Monte Carlo drawing """
    
    def __init__(self, serie=None, ndraw=500):
        """
        Class builder.
        Set data as attribute.
        Draw Monte Carlo data.
        
        Parameters
        ----------
        serie : [pandas.Series]
            Input data in a compatible LePhare format. Here are the keys style to adopt:
                filtername0, filtername0.err, filtername1, filtername1.err, ... filternameN, filternameN.err, CONTEXT, Z-SPEC, STRING
            where 'filtername{}' is a known filter by LePhare with the format instrument.band, for instance sdss.u or ps1.z, given in flux.
            The 'CONTEXT' sets the filters to use among the provided ones in input data for the SED fitting in each row.
            It is the sum of 2 to the {i}, {i} being the filtername numbers to use in each row.
            For example, using the filters = ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] in data input:
                * ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3 + 2^4 = 31
                * ["sdss.g", "sdss.r", "sdss.i", "sdss.z"]           --> CONTEXT = 2^1 + 2^2 + 2^3 + 2^4       = 30
                * ["sdss.u", "sdss.g", "sdss.r", "sdss.i"]           --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3       = 15
                * ["sdss.u", "sdss.i", "sdss.z"]                     --> CONTEXT = 2^0 + 2^3 + 2^4             = 25
            If 'CONTEXT' column doesn't exist, automatically create one with the maximum CONTEXT value.
            'Z-SPEC' is the spectroscopic redshift, if available.
            'STRING' and all following columns are not used for the SED fitting. Allow the user to store additional informations.
        
        ndraw : [int]
            Number of Monte Carlo draws.
        
        
        Returns
        -------
        Void
        """
        if serie is not None:
            self.set_data(serie)
        if self.has_data() and (ndraw is not None and ndraw>1):
            self.draw_mc(ndraw)

    @classmethod
    def load_fromdir(cls, dirin, verbose=True):
        """
        Build an object from MCLePhare output results directory.
        Load the resulting spectra and the class attributes.
        
        Parameters
        ----------
        dirin : [string]
            MCLePhare output results directory.
            The folder must contain 'spec', 'config', 'catin' and 'catout'.
        
        Options
        -------
        verbose : [bool]
            Print informations.
            Default is True.
        
        
        Returns
        -------
        MCLePhare
        """
        import os
        try:
            datain = {"spec":[dirin+"/"+l for l in os.listdir(dirin) if l.endswith(".spec")],
                      "config":[dirin+"/"+l for l in os.listdir(dirin) if "config" in l][0],
                      "catin":[dirin+"/"+l for l in os.listdir(dirin) if "data" in l or "catin" in l][0],
                      "catout":[dirin+"/"+l for l in os.listdir(dirin) if "catout" in l][0]}
        except:
            raise IOError("Cannot find the spec, config, catin and catout from this directory")

        this = cls()
        this._lephare_out = datain
        this._load_results_(updatefilters=True, verbose=verbose)
        return this
        
    # ============== #
    #  Methods       #
    # ============== #
    def run(self, configfile=None, dirout=None, onwhat="gal", gallib="BC03", **kwargs):
        """
        Run the SED fitting on Monte Carlo generated data.
        If LePhare has not be loaded with the draws, automatically do it.
        
        Options
        -------
        configfile : [string or None]
            New configuration file directory.
            If None, use the one already set for this object.
            Default is None.
        
        dirout : [string or None]
            Modify the output folder directory.
            If None, use the one already set for this object.
            Default is None.
        
        onwhat : [string or list or None]
            // ignored if None //
            Defines on which kind of templates to run the SED fitting.
            List which could contain "gal", "star", "qso" (any or combination of).
            For example: ["gal", "qso"]
            Default is ["gal","star","qso"].
        
        gallib : [string]
            If 'gal' in 'onwhat', define the used galaxy library among the available ones ($LEPHAREDIR/sed/GAL/).
            Default is "BC03".
        
        **kwargs
        filters : [list or None]
            List of filters on which to execute the SED fitting, meaning that the context is modified.
            Must be known filters by LePhare with the format instrument.band, for instance sdss.u or ps1.z.
            Go with 'contextid' parameter.
            If None, the context is not changed from what has been set before this point (you can check with {self}.data).
            Default is None.
        
        contextid : [int or list(int) or None]
            Row index(es) to set the new context.
            Go with 'filters' parameter (setting the new context value).
            If None, set the new context value for the whole data table.
            Default is None.
        
        catinfile : [string or None]
            New input data catalog (automatically change it in the configuration file too).
            Go with 'configfile' parameter (must be None).
            If None, use the one already set for this object.
            Default is None.
        
        originalconfig : [bool]
            If True, use the original configuration file ; False acounts for any brought modifications at this point.
            Go with 'configfile' parameter (must be None).
            Default is False.
        
        update_init : [bool]
            If True, run the LePhare initialization, even if the templates already exist.
            Default is False.
        
        verbose : [bool]
            Print informations.
            Default is True.
        
        
        Returns
        -------
        Void
        """
        if not self.has_lephare() or configfile is not None:
            self.load_lephare(configfile=configfile, dirout=dirout)
            
        self._lephare_out = self.lephare.run(onwhat=onwhat, gallib=gallib, configfile=configfile, dirout=dirout, **kwargs)
    
    # ------- #
    # GETTER  #
    # ------- #
    def get_spectrum(self, index):
        """
        Returns the index-th LePhareSpectrum
        
        Parameters
        ----------
        index : [int]
            Index of the desired spectrum.
        
        
        Returns
        -------
        LePhareSpectrum
        """
        return self.spectra.spectra[index]
    
    def get_synthetic_photometry(self, filter_, restframe=False, influx=True, inhz=False):
        """
        Return photometry synthesized through the given filter/bandpass.
        The returned data are (effective wavelength, synthesize flux/mag) in an array with same size as for the number of draws.

        Parameters
        ----------
        filter_ : [string, sncosmo.BandPass, 2D array]
            The filter through which the spectrum will be synthesized.
            Accepted input format:
            - string: name of a known filter (instrument.band), the actual bandpass will be grabbed using io.get_filter_bandpass(filter_)
            - 2D array: bandpass = sncosmo.bandpass.BandPass(*filter_)
            - sncosmo.bandpass.BandPass
        
        Options
        -------
        restframe : [bool]
            If True, the spectra are first deredshifted before doing the synthetic photometry.
            Default is False.

        influx : [bool]
            Should the results be returned in flux (True) or in ABmag (False).
            Default is True (thus in flux).

        inhz : [bool]
            // ignored if influx = False //
            Should the flux be returned in erg/s/cm2/Hz (True) or erg/s/cm2/AA (False).
            Default is False (thus in erg/s/cm2/AA).
        
        
        Returns
        -------
        np.array, np.array
        """
        return self.spectra.get_synthetic_photometry(filter_, restframe=restframe, influx=influx, inhz=inhz)
    
    # ------- #
    # SETTER  #
    # ------- #
    def set_data(self, serie):
        """
        Set the input data and the available filters as attribute.
        
        Parameters
        ----------
        serie : [pandas.Series]
            Input data in a compatible LePhare format. Here are the keys style to adopt:
                filtername0, filtername0.err, filtername1, filtername1.err, ... filternameN, filternameN.err, CONTEXT, Z-SPEC, STRING
            where 'filtername{}' is a known filter by LePhare with the format instrument.band, for instance sdss.u or ps1.z, given in flux.
            The 'CONTEXT' sets the filters to use among the provided ones in input data for the SED fitting in each row.
            It is the sum of 2 to the {i}, {i} being the filtername numbers to use in each row.
            For example, using the filters = ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] in data input:
                * ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3 + 2^4 = 31
                * ["sdss.g", "sdss.r", "sdss.i", "sdss.z"]           --> CONTEXT = 2^1 + 2^2 + 2^3 + 2^4       = 30
                * ["sdss.u", "sdss.g", "sdss.r", "sdss.i"]           --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3       = 15
                * ["sdss.u", "sdss.i", "sdss.z"]                     --> CONTEXT = 2^0 + 2^3 + 2^4             = 25
            If 'CONTEXT' column doesn't exist, automatically create one with the maximum CONTEXT value.
            'Z-SPEC' is the spectroscopic redshift, if available.
            'STRING' and all following columns are not used for the SED fitting. Allow the user to store additional informations.
        
        
        Returns
        -------
        Void
        """
        if not type(serie) == pandas.core.series.Series:
            raise TypeError("insut serie must be pandas.Serie.")
        self._data = serie
        filters = io.keys_to_filters(serie.keys())
        self.set_filters(filters) 

    def set_intrinsic_scatter(self, magerr, onconfig=False, inhz=False, redraw=True):
        """
        Set the intrinsic scatter.
        This scatter could be added during the LePhare fit itself or when drawing the Monte Carlo.
        For the latter be careful of which flux unit is used (you can play on 'inhz')

        Parameters
        ----------
        magerr : [1D array]
            Magnitude intrinsic scatter.
            Applied independently on each filter, must be on the same size as for self.filters.
        
        Options
        -------
        onconfig : [bool]
            Can change the 'ERR_SCALE' parameter in the LePhare configuration file.
            - False: Individual data error will be increased by the given magnitude error
                     and the ERR_SCALE config parameter will be set to None.
            - True: Individual data error are unchanged, magerr passed to ERR_SCALE.
            Default is False.
        
        inhz : [bool]
            // ignored if onconfig = True //
            Are the input fluxes in erg/s/cm2/Hz (True) or erg/s/cm2/AA (False).
            Default is False (thus in erg/s/cm2/AA).
            
        redraw: [bool]
            // Ignored if 'onconfig' is True //
            If True, redraw the mcdata.
            Default is True.
        
        
        Returns
        -------
        Void
        """
        if len(magerr) != self.nfilters:
            raise ValueError("You must provide exactly one magerr per filterband\n"+
                             "(magerr = {} vs. filters = {}).".format(magerr, self.filters))
        if onconfig:
            if not self.has_lephare():
                raise AttributeError("lephare is not loaded.")
            self.lephare.config.set_intrinsic_error(magerr)
            self._sigma_int = None
        else:
            #raise ValueError("Not yet complete. Mag have to be convert in flux...")
            from . import tools
            _lbda = [self.filter_bandpasses[_f].wave_eff for _f in self.filters]
            fluxerr = tools.mag_to_flux(magerr, magerr=None, wavelength=_lbda, inhz=inhz)
            self._sigma_int = {k:v for k,v in zip(self.filters,magerr)}
            if redraw:
                if self.ndraw is not None:
                    self.draw_mc(self.ndraw)
                else:
                    self.draw_mc()
            
    def draw_mc(self, ndraw=500):
        """
        Draw Monte Carlo on the loaded fluxes.
        Apply intrinsic scatter (if any).
        Format the data into LePhare compatible dataframe and set it as attribute.
        
        Parameters
        ----------
        ndraw : [int]
            Number of Monte Carlo draws.
        
        
        Returns
        -------
        Void
        """
        data_mc = {}
        for filt_ in self.filters:
            if self.sigma_int is None:
                filt_err = self.data[filt_+".err"]
            else:
                filt_err = np.sqrt(self.data[filt_+".err"]**2  + self.sigma_int[filt_]**2)
            data_mc[filt_] = np.random.normal(loc=self.data[filt_], scale=filt_err, size=ndraw)
            data_mc[filt_+".err"] = filt_err
            #Save of the measurement origins
            data_mc[filt_][0] = self.data[filt_]
            
        mcdata = pandas.DataFrame(data_mc)
        for k in self.data.keys():
            if k not in mcdata.columns:
                mcdata[k] = self.data[k]
                
        self._mcdata = mcdata


    # ------- #
    # LOADER  #
    # ------- #
    def load_lephare(self, configfile, dirout=None, **kwargs):
        """
        Load and set a LePhare object as attribute based on the Monte Carlo data..
        
        Parameters
        ----------
        configfile : [string or None]
            Configuration file directory.
            If None, use the default one ($LEPHAREWORK/pylephare/config/default.config).
        
        dirout : [string or None]
            Setting the output folder directory.
            If None, create a default directory based on date and time in $LEPHAREWORK/pylephare/.
            Default is None.
        
        **kwargs
        inhz : [bool]
            Set to True if the fluxes (and flux uncertainties) are given in erg/s/cm2/Hz ; False means in erg/s/cm2/AA.
            Default is False.
        
        verbose : [bool]
            Print informations.
            Default is True.
        
        
        Returns
        -------
        Void
        """
        from . import lephare
        self._lephare = lephare.LePhare(self.mcdata, configfile, dirout=dirout, **kwargs)
        
    def _load_results_(self, updatefilters=False, lbda_range=[1000,10000], verbose=True):
        """
        Load the LePhare results and spectra.
        If some MC draws have failed the SED fitting, redraw them and run again the SED fitting.
        
        Options
        -------
        updatefilters : [bool]
            If True, reset the filters from the 'FILTER_LIST' parameter's value in the LePhare configuration file.
            Default is False.
        
        lbda_range : [list(float or None) or None]
            Wavelength range (in AA) to load for the spectra.
            If None, the full spectrum is recovered. You can also set None to one of the limit to let it be the extremum.
            Default is [1000,10000].
        
        verbose : [bool]
            Print informations.
            Default is True.
        
        
        Returns
        -------
        Void
        """
        if verbose:
            print("loading results...")
            
        from . import configparser, spectrum, lephare, tools
        self._config = configparser.ConfigParser(self._lephare_out["config"])
        if updatefilters:
            self.set_filters(self._config.get_filters(name=True))
            
        self._catin = pandas.read_csv(self._lephare_out["catin"], sep=" ", index_col=0, names=self._config.get_catin_columns())
        # mcdata in AA while _catin in hz
        if not self.has_mcdata():
            self._mcdata = pandas.DataFrame({k:tools.flux_hz_to_aa( self._catin[k], self.filter_bandpasses[k.replace(".err","")].wave_eff)
                                             if k.replace(".err","") in self.filters else self._catin[k] for k in self._catin.columns})
        if not self.has_data():
            self._data = self._mcdata.iloc[0].copy()
            self._data.drop("context", inplace=True)
            self._data.rename({"z-spec":"Z-SPEC", "string":"STRING"}, inplace=True)
        
        self._catout = lephare.read_catout(self._lephare_out["catout"], self.filters)
        try:
            self._spectra = spectrum.LePhareSpectrumCollection.read_files(self._lephare_out["spec"], lbda_range=lbda_range)
        except ValueError:
            self._replace_failed_mc_()
    
    def _replace_failed_mc_(self):
        """
        Detect SED fitting failures by checking 'Z_BEST' values in LePhare output catalog.
        Redraw Monte Carlo data for these failures and run again the SED fitting.
        Automatically run again _load_results_.
        
        
        Returns
        -------
        Void
        """
        data_mc = self._mcdata.copy()
        dirout = self._lephare_out["config"].split("/config")[0]
        list_idx_fails = [ii for ii, k in enumerate(self.catout.Z_BEST.values) if k == "-99.0000"]
        for ii in list_idx_fails:
            for filt_ in self.filters:
                data_mc.loc[ii, filt_] = np.random.normal(loc=self.data[filt_], scale=self.data[filt_+".err"])
        self._mcdata = data_mc
        self.load_lephare(configfile=self._lephare_out["config"], dirout=dirout)
        self.run(verbose=False, update=False)
        self._load_results_(verbose=False, updatefilters=False)
        
            

    # ------- #
    # PLOTTER #
    # ------- #
    def show(self, ax=None, figsize=[7,3.5], ax_rect=[0.1,0.2,0.8,0.7], nmc=10, add_spec=True, inhz=False):
        """
        Plot the synthesized photometry and/or spectrum.
        Return a dictionary containing the figure and the axes ({"fig":fig, "ax":ax}).
        
        Options
        -------
        ax : [plt.Axes or None]
            If an axes is given, draw the spectrum on it.
            Default is None.
        
        figsize : [list(float)]
            Two values list to fix the figure size (in inches).
            Default is [7,3.5].
        
        ax_rect : [list(float)]
            The dimensions [left, bottom, width, height] of the new axes.
            All quantities are in fractions of figure width and height.
            Default is [0.1,0.2,0.8,0.7].
        
        nmc : [int]
            Number of MC draws to overplot (spectra and/or synthesized photometry points).
            Default is 10.
        
        add_spec : [bool]
            If True, plot the spectrum.
            Default is True.
        
        inhz : [bool]
            Should the flux be ploted in erg/s/cm2/Hz (True) or erg/s/cm2/AA (False).
            Default is False (thus in erg/s/cm2/AA).
        
        
        Returns
        -------
        dict
        """
        if inhz:
            raise NotImplementedError("inhz plotting option is not implemented yet.")
            
        import matplotlib.pyplot as mpl
        if ax is None:
            fig = mpl.figure(figsize=figsize)
            ax = fig.add_axes(ax_rect)
        else:
            fig = ax.figure
        
        prop = dict(ls="None", marker="o", ecolor="0.7", mfc="C0", mec="w", ms=10, zorder=5)
        propmc = dict(ls="None", marker="o", mfc="0.7", mec="None", ms=5, alpha=0.3, zorder=3)
        idmc = np.random.choice(np.arange(self.ndraw), replace=False, size=nmc)
        for f in self.filters:
            if self.has_data(): # not single data is loaded from results
                ax.errorbar(self.filter_bandpasses[f].wave_eff, self.data[f], 
                            yerr=self.data[f+".err"], **prop)
            
            ax.plot([self.filter_bandpasses[f].wave_eff]*nmc, self.mcdata[f][idmc],  **propmc)
            
        if add_spec:
            _ = [self.spectra.spectra[i_].show(ax=ax, showdata=False, showmagmodel=False,
                                               set_label=False, lw=1, alpha=0.3)
                 for i_ in idmc]
        
        ax.set_xlabel(r"Wavelentgh [$\AA$]", fontsize="large")
        ax.set_ylabel(r"flux [$\mathrm{erg\,s^{-1}\,cm^{-2}\,%s}$]"%("Hz^{-1}" if inhz else "\AA^{-1}"),  fontsize="large")
        
        return {"fig":fig, "ax":ax}

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
        """ Test that data and not empty data are loaded """
        return self.data is not None and len(self.data)>0

    @property
    def sigma_int(self):
        """ Intrinsic magnitude scatter to be applied to the data when drawing """
        if not hasattr(self, "_sigma_int"):
            self._sigma_int = None
        return self._sigma_int
        
    @property
    def mcdata(self):
        """ DataFrame containing the Monte Carlo drawed data """
        if not hasattr(self,"_mcdata"):
            self._mcdata = None
        return self._mcdata

    def has_mcdata(self):
        """ Test that MC data and not empty MC data are loaded """
        return self.mcdata is not None and len(self.mcdata)>0

    @property
    def ndraw(self):
        """ Monte Carlo draws number """
        return len(self.mcdata) if self.has_mcdata() else None
        
    # - LePhare
    @property
    def lephare(self):
        """ LePhare object """
        if not self.has_lephare():
            raise AttributeError("LePhare is not loaded, see self.load_lephare()")
        return self._lephare
    
    def has_lephare(self):
        """ Test that a LePhare object is loaded """
        return hasattr(self, "_lephare") and self._lephare is not None

    def _did_run_(self):
        """ test if you ran self.run()"""
        return hasattr(self, "_lephare_out")

    # - lephare.run() results
    @property
    def catout(self):
        """ Output catalog """
        if not hasattr(self, "_catout"):
            if self._did_run_():
                self._load_results_()
            else:
                raise AttributeError("you did not run self.run()")
        return self._catout
    
    @property
    def config(self):
        """ Configuration file used to run() ; ConfigParser object """
        if not hasattr(self, "_config"):
            if self._did_run_():
                self._load_results_()
            else:
                raise AttributeError("you did not run self.run()")
        return self._config

    @property
    def spectra(self):
        """ spectra estimated from run() ; LePhareSpectrumCollection object """
        if not hasattr(self, "_spectra"):
            if self._did_run_():
                self._load_results_()
            else:
                raise AttributeError("you did not run self.run()")
        return self._spectra
