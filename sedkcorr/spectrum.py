""" LePhare Spectrum files """

import pandas
import numpy as np
from . import tools, io


def synthesize_photometry(lbda, flux, filter_lbda, filter_trans,
                          normed=True):
    """ Get Photometry from the given spectral information through the given filter.

    This function converts the flux into photons since the transmission provides the
    fraction of photons that goes though.


    Parameters
    -----------
    lbda, flux: [array]
        Wavelength and flux of the spectrum from which you want to synthetize photometry
        
    filter_lbda, filter_trans: [array]
        Wavelength and transmission of the filter.

    normed: [bool] -optional-
        Shall the fitler transmission be normalized?

    Returns
    -------
    Float (photometric point)
    """
    # ---------
    # The Tool
    def integrate_photons(lbda, flux, step, flbda, fthroughput):
        """ """
        filter_interp = np.interp(lbda, flbda, fthroughput)
        dphotons = (filter_interp * flux) * lbda * 5.006909561e7
        return np.trapz(dphotons,lbda) if step is None else np.sum(dphotons*step)
    
    # ---------
    # The Code
    normband = 1. if not normed else \
      integrate_photons(lbda, np.ones(len(lbda)),None,filter_lbda,filter_trans)
      
    return integrate_photons(lbda,flux,None,filter_lbda,filter_trans)/normband

# ----------- #
# Redshifting #
# ----------- #
def deredshift(lbda, flux, z, variance=None, exp=3):
    """Deredshift spectrum from z to 0, and apply a (1+z)**exp
    flux-correction.

    exp=3 is for erg/s/cm2/A spectra to be latter corrected using proper
    (comoving) distance but *not* luminosity distance.

    Parameters
    ----------
    lbda, flux : [array]
        wavelength and flux or the spectrum. 
        y is expected in erg/s/cm2/A (so x in A)
        
    z: [float]
        Cosmological redshift

    variance: [array] -optional-
        spectral variance if any (square of flux units.)

    exp: [float] -optional-
        exposant for the redshift flux dilution 
        [exp=3 is for erg/s/cm2/A spectra]

    Returns
    -------
    arrays (x, y [and variance if any] corresponding to the deredshifted spectrum)

    """
    zp1 = 1 + z
    lbda_rest = lbda/zp1           # Wavelength correction
    flux_rest = flux*(zp1**exp)      # Flux correction
    if variance is not None:
        variance_rest = variance*(zp1**exp)**2

    return (lbda_rest, flux_rest , variance_rest) if variance is not None else (lbda_rest, flux_rest)


#####################
#                   #
#  Spectra Classes  #
#                   #
#####################


class LePhareSpectrum( object ):
    """ """
    def __init__(self, filename=None, lbda_range=None):
        """ """
        if filename is not None:
            self.load(filename)
            if lbda_range is not None:
                self.set_lbda_range(*lbda_range)
    
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

    def get_redshift(self, spectro=True):
        """ """
        return self.zspec if spectro else self.zphoto
    
    def get_spectral_data(self, model="gal", influx=True, inhz=False, restframe=False):
        """ """
        data_ = self.get_data(model)
        lbda_, mag_ = data_.lbda, data_.mag
        flux_aa = tools.mag_to_flux(mag_, None, wavelength=lbda_, inhz=False)[0]
        
        if restframe:
            lbda_, flux_aa = deredshift(lbda_, flux_aa, self.get_redshift(spectro=True))
            
        if influx:
            return lbda_, flux_aa if not inhz else tools.flux_aa_to_hz(flux_aa, lbda_)
        
        mag_ = tools.flux_to_mag(flux_aa, None, wavelength=lbda_, inhz=False)[0]
        return lbda_, mag_
    
    def get_input_data(self, influx=True, inhz=False):
        """ """
        lbda_, mag_, emag_ = np.asarray(self.resultmags[["Lbd_mean", "Mag", "emag"]].values, dtype="float").T
        if influx:
            return lbda_, tools.mag_to_flux(mag_, magerr=emag_, wavelength=lbda_, inhz=inhz)
        return lbda_, [mag_, emag_]
    
    def get_model_data(self, influx=True, inhz=False):
        """ """
        lbda_, mag_ = np.asarray(self.resultmags[["Lbd_mean", "Mag_gal"]].values, dtype="float").T
        if influx:
            return lbda_, tools.mag_to_flux(mag_, None, wavelength=lbda_, inhz=inhz)[0]
        return lbda_, mag_
    
    def _get_model_nlines_(self, model):
        """ """
        rangelist_ = [0]+list(np.cumsum(self.resultmodels["Nline"].astype("int")))
        rangelist = [[rl_,rm_] for rl_,rm_ in zip(rangelist_[:-1],rangelist_[1:])]
        return rangelist[self.resultmodels.index.get_loc(model)]

    # -------- #
    #  Main    #
    # -------- #
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
        from sncosmo import bandpasses
        
        # - Get the corresponding bandpass
        if type(filter_) == str:
            bp = io.get_filter_bandpass(filter_)
        elif type(filter_) = bandpasses.Bandpass:
            bp = filter_
        elif len(filter_) == 2:
            bp = bandpasses.Bandpass(*filter_)
        else:
            raise TypeError(" filter_ must either be a filter name like filter_='sdss.u', a sncosmo.BandPass or a 2D array filter_=[wave, transmission]")

        # - Synthesize through bandpass
        sflux_aa = self.synthesize_photometry(bp.wave, bp.trans, restframe=restframe)
        slbda_ = bp.wave_eff/(1+self.get_redshift()) if restframe else bp.wave_eff

        if influx:
            return slbda_, (sflux_aa if not inhz else tools.flux_aa_to_hz(sflux_aa, slbda_))
        
        mag_ = tools.flux_to_mag(sflux_aa, None, wavelength=slbda_)[0]
        return slbda_, mag_
    
    def synthesize_through_filter(self, filtername,  restframe=False, influx=True, inhz=False):
        """ """
        return self.get_synthetic_photometry(filtername,  restframe=False, influx=True, inhz=False)
    
    def synthesize_photometry(self, filter_lbda, filter_trans, model="gal", restframe=False):
        """ get the synthetic flux in AA """
        lbda_, flux_aa = self.get_spectral_data(model="gal", influx=True, inhz=False, restframe=restframe)
        return synthesize_photometry(lbda_, flux_aa, filter_lbda, filter_trans, normed=True)
    

    # -------- #
    # PLOTTER  #
    # -------- #        
    def show(self, ax=None, influx=True, inhz=False, model="gal", 
             showdata=True, showmagmodel=True, restframe=False,
             scprop={}, colormodel=".7", set_label=True, **kwargs):
        """ """
        import matplotlib.pyplot as mpl
        if ax is None:
            fig = mpl.figure(figsize=[7,3.5])
            ax = fig.add_axes([0.1,0.2,0.8,0.7])
        else:
            fig = ax.figure
            
        propfunc = dict(influx=influx, inhz=inhz)
            
        lbda_, fluxmag_ = self.get_spectral_data(model=model,restframe=restframe, **propfunc)
        ax.plot(lbda_, fluxmag_, color=colormodel, **kwargs)
        if showdata:
            lbda_, [data, err] = self.get_input_data(**propfunc)
            prop_ = dict(ls="None", marker="o", mfc="C0", mec="C0", ecolor="0.7")
            ax.errorbar(lbda_, data, yerr=err, **{**prop_,**scprop})
        if showmagmodel:
            lbda_, data = self.get_model_data(**propfunc)
            prop_ = dict(ls="None", marker="o", mfc="w", mec=colormodel, mew=1)
            ax.plot(lbda_, data, **{**prop_,**scprop})   

        if set_label:
            ax.set_xlabel(r"Wavelentgh [$\AA$]", fontsize="large")
            ax.set_ylabel(r"flux [$\mathrm{erg\,s^{-1}\,cm^{-2}\,%s}$]"%("Hz^{-1}" if inhz else "\AA^{-1}") 
                      if influx else "magnitude",  fontsize="large")
            
        return fig
    
    # ============== #
    #   Properties   #
    # ============== #
    @property
    def zphoto(self):
        """ """
        return self._zphot

    @property
    def zspec(self):
        """ """
        return self._zspec
    
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



class LePhareSpectrumCollection( object ):
    def __init__(self, lepharespectra):
        """ """
        self._spectra = lepharespectra

    @classmethod
    def read_files(cls, filenames, lbda_range=None):
        """ """
        return cls([LePhareSpectrum(file_, lbda_range) for file_ in filenames])
    
    # -------- #
    #  SETTER  #
    # -------- #
    def set_lbda_range(self, lbdamin=None, lbdamax=None):
        """ """
        _ = [spec.set_lbda_range(lbdamin=None, lbdamax=None) for spec in self._spectra]

    # -------- #
    #  GETTER  #
    # -------- #
    def get_data(self, model="gal", apply_lbdarange=True):
        """ """
        return [spec.get_data(model=model, apply_lbdarange=apply_lbdarange) for spec in self._spectra]

    def get_redshift(self, spectro=True):
        """ """
        return [spec.get_redshift(spectro=spectro) for  spec in self._spectra]
    
    def get_spectral_data(self, model="gal", influx=True, inhz=False, restframe=False):
        """ """
        return [spec.get_spectral_data(model=model, influx=influx, inhz=inhz, restframe=restframe) for  spec in self._spectra]
    
    def get_input_data(self, influx=True, inhz=False):
        """ """
        return [spec.get_input_data(influx=influx, inhz=inhz) for  spec in self._spectra]
    
    def get_model_data(self, influx=True, inhz=False):
        """ """
        return [spec.get_model_data(influx=influx, inhz=inhz) for  spec in self._spectra]

    # -------- #
    #  Main    #
    # -------- #
    def get_synthetic_photometry(self, filter_,  restframe=False, influx=True, inhz=False):
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
        return [spec.synthesize_through_filter(filter_, restframe=restframe, influx=influx, inhz=inhz) for  spec in self._spectra]
    
    def synthesize_photometry(self, filter_lbda, filter_trans, model="gal", restframe=False):
        """ get the synthetic flux in AA """
        return [spec.synthesize_photometry(filter_lbda, filter_trans, model=model, restframe=restframe) for  spec in self._spectra]

    # ============== #
    #   Properties   #
    # ============== #
    @property
    def spectra(self):
        """ """
        return self._spectra
