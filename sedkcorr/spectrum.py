""" LePhare Spectrum files """

import pandas
import numpy as np
from .utils import tools

class LePhareSpectrum( object ):
    """ """
    def __init__(self, filename=None):
        """ """
        if filename is not None:
            self.load(filename)
    
    def load(self, filename, onlyfirstrow=True):
        """ """
        datain = open(filename).read().splitlines()
        id_, zspec, zphot = datain[1].split()
        self._id = int(id_)
        self._zspec = float(zspec)
        self._zphot = float(zphot)
        self._nfilters = int(datain[3].split()[1])
        
        # Model 
        col_model = datain[6].split()[1:]
        indexes_model = [i for i,d in enumerate(datain) if np.any([d.startswith(v) for v in ["GAL","STAR","QSO"]])]
        self._dfmodel = pandas.DataFrame([d.split() for d in np.asarray(datain)[indexes_model]], 
                                         columns=col_model).set_index("Type")
        
        # Magnitudes
        col_mag = datain[2].split()[1:]
        self._dfmag = pandas.DataFrame([d.split() for i,d in enumerate(datain) if len(d.split())==len(col_mag)], 
                                         columns=col_mag)
        
        i_base = len(self._dfmag)+len(self._dfmodel)
        istart = np.min([i_base+i for i,l in enumerate(datain[i_base:]) if len(l.split())==2])
        
        self.data = pandas.read_csv(filename,
                                        skiprows=istart, names=["lbda", "mag"], sep="  ",
                                        engine="python")
    # -------- #
    #  SETTER  #
    # -------- #    
    def set_lbda_range(self, lbdamin=None, lbdamax=None):
        """ """
        self._flagin = (self.data.lbda>lbdamin) & (self.data.lbda<lbdamax)
        
        
    # -------- #
    #  GETTER  #
    # -------- # 
    def get_data(self, model="gal", apply_lbdarange=True):
        """ """
        if model in ["gal"]:
            model = "GAL-1"
        index_start, index_end = self._get_model_nlines_(model)
        datas = self.data.iloc[index_start:index_end]
        return datas if not apply_lbdarange else datas[self.flagin[index_start:index_end]]
        
    def get_spectral_data(self, model="gal", influx=True, unit="AA"):
        """ """
        data_ = self.get_data(model)
        lbda_, mag_ = data_.lbda, data_.mag
        if influx:
            return lbda_, tools.mag_to_flux(mag_, band=lbda_, flux_unit=unit)[0]
        return lbda_, mag_
    
    def get_input_data(self, influx=True, unit="AA"):
        """ """
        lbda_, mag_, emag_ = np.asarray(spec.resultmags[["Lbd_mean", "Mag", "emag"]].values, dtype="float").T
        if influx:
            return lbda_, tools.mag_to_flux(mag_, mag_err=emag_, band=lbda_, flux_unit=unit)
        return lbda_, [mag_, emag_]
    
    def get_model_data(self, influx=True, unit="AA"):
        """ """
        lbda_, mag_ = np.asarray(spec.resultmags[["Lbd_mean", "Mag_gal"]].values, dtype="float").T
        if influx:
            return lbda_, tools.mag_to_flux(mag_, band=lbda_, flux_unit=unit)[0]
        return lbda_, mag_
    
    def _get_model_nlines_(self, model):
        """ """
        rangelist_ = [0]+list(np.cumsum(self.resultmodels["Nline"].astype("int")))
        rangelist = [[rl_,rm_] for rl_,rm_ in zip(rangelist_[:-1],rangelist_[1:])]
        return rangelist[self.resultmodels.index.get_loc(model)]

    # -------- #
    # PLOTTER  #
    # -------- #        
    def show(self, ax=None, influx=True, unit="AA", model="gal", 
             scprop={}, **kwargs):
        """ """
        import matplotlib.pyplot as mpl
        if ax is None:
            fig = mpl.figure(figsize=[6,4])
            ax = fig.add_axes([0.12,0.12,0.8,0.8])
        else:
            fig = ax.figure
            
        propfunc = dict(influx=influx, unit=unit)
            
        lbda_, fluxmag_ = self.get_spectral_data(model=model, **propfunc)
        ax.plot(lbda_, fluxmag_*100, **kwargs)
        #if showdata:
        #    lbda_, [data, err] = self.get_input_data(**propfunc)
        #    prop_ = dict(ls="None", marker="o", mfc="C0", mec="C0", ecolor="0.7")
        #    ax.errorbar(lbda_, data, yerr=err, **{**prop_,**scprop})
        #if showmagmodel:
        #    lbda_, data = self.get_model_data(**propfunc)
        #    prop_ = dict(ls="None", marker="o", mfc="None", mec="C0", ecolor="0.7")
        #    ax.errorbar(lbda_, data, yerr=err, **{**prop_,**scprop})   
        return fig
    
    # ============== #
    #   Properties   #
    # ============== #
    @property
    def mag(self):
        """ """
        return self.data.mag[self.flagin]
    
    @property
    def flagin(self):
        """ """
        if not hasattr(self, "_flagin"):
            self.set_lbda_range(0,1e8)
        return self._flagin
    
    #
    # - Results
    #
    @property
    def resultmodels(self):
        """ """
        return self._dfmodel
    
    @property
    def resultmags(self):
        """ """
        return self._dfmag
