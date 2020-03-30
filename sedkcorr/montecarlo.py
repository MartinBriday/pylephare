""" Handles the target SEDfitting Monte Carlo """

import pandas
import numpy as np
from . import io, base

class MCSED( base._FilterHolder_ ):
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
    # SETTER  #
    # ------- #
    def set_data(self, serie):
        """ """
        if not type(serie) == pandas.core.series.Series:
            raise TypeError("insut serie must be pandas.Serie.")
        self._data = serie
        filters = io.keys_to_filters(serie.keys())
        self.set_filters(filters) 
        
    def draw_mc(self, ndraw=500):
        """ """
        data_mc = {k:np.random.normal(loc=self.data[k], scale= self.data[k+".err"], size=ndraw)
                       for k in self.filters}
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
            
        from sedkcorr import configparser, spectrum, lephare
        self._config = configparser.ConfigParser(self._lephare_out["config"])
        self._catin = pandas.read_csv(self._lephare_out["catin"], sep=" ")
        self._catout = lephare.read_catout(self._lephare_out["catout"], self.filters)
        self._spectra = [spectrum.LePhareSpectrum(spec_, lbda_range=[3000,10000]) for spec_ in self._lephare_out["spec"]]
        
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
            _ = [self._spectra[i_].show(ax=ax, showdata=False, showmagmodel=False, set_label=False)
                for i_ in idmc]
            
        
        ax.set_xlabel(r"Wavelentgh [$\AA$]", fontsize="large")
        ax.set_ylabel(r"flux [$\mathrm{erg\,s^{-1}\,cm^{-2}\,%s}$]"%("Hz^{-1}" if inhz else "\AA^{-1}") 
                      if influx else "magnitude",  fontsize="large")

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
