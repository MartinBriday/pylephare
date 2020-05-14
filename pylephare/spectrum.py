""" LePhare Spectrum files """

import pandas
import numpy as np
from . import tools, io


def synthesize_photometry(lbda, flux, filter_lbda, filter_trans, normed=True):
    """
    Return photometry from the given spectral information through the given filter.

    This function converts the flux into photons since the transmission provides the fraction of photons that goes though.

    Parameters
    -----------
    lbda, flux : [array]
        Wavelength and flux of the spectrum from which you want to synthetize photometry.
        
    filter_lbda, filter_trans : [array]
        Wavelength and transmission of the filter.
    
    Options
    -------
    normed : [bool]
        Shall the fitler transmission be normalized (True) or not (False).
        Default is True.

    Returns
    -------
    float (photometric point)
    """
    # ---------
    # The Tool
    def integrate_photons(lbda, flux, step, flbda, fthroughput):
        filter_interp = np.interp(lbda, flbda, fthroughput)
        dphotons = (filter_interp * flux) * lbda * 5.006909561e7
        return np.trapz(dphotons,lbda) if step is None else np.sum(dphotons*step)
    
    # ---------
    # The Code
    normband = 1. if not normed else integrate_photons(lbda, np.ones(len(lbda)),None,filter_lbda,filter_trans)
      
    return integrate_photons(lbda,flux,None,filter_lbda,filter_trans)/normband

# ----------- #
# Redshifting #
# ----------- #
def deredshift(lbda, flux, z, variance=None, exp=3):
    """
    Deredshift spectrum from z to 0, and apply a (1+z)**exp flux-correction.
    
    exp=3 is for erg/s/cm2/A spectra to be later corrected using proper (comoving) distance but *not* luminosity distance.
    
    Return the restframed wavelength, flux and, if any, variance.

    Parameters
    ----------
    lbda, flux : [array]
        Wavelength and flux or the spectrum.
        Flux is expected to be in erg/s/cm2/AA.
        
    z : [float]
        Cosmological redshift

    variance : [array or None]
        Spectral variance if any (square of flux units).
        Default is None.

    exp : [float]
        Exposant for the redshift flux dilution (exp=3 is for erg/s/cm2/AA spectra).
        Default is 3.

    Returns
    -------
    np.array, np.array, np.array if variance is not None
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
    """ Class handling the resulting spectrum from LePhare SED fitting. """
    
    def __init__(self, filename=None, lbda_range=None):
        """
        Class builder.
        
        Parameters
        ----------
        filename : [string or None]
            // Ignored if None //
            Spectrum file directory resulting from LePhare SED fitting.
        
        Options
        -------
        lbda_range : [list(float) or None]
            // Ignored if None or if filename is None //
            Wavelength range (in AA) to load for the spectra.
        
        
        Returns
        -------
        Void
        """
        if filename is not None:
            self.load(filename)
            if lbda_range is not None:
                self.set_lbda_range(*lbda_range)
    
    def load(self, filename):
        """
        Load the spectrum from the given file directory.
        
        Parameters
        ----------
        filename : [string or None]
            Spectrum file directory resulting from LePhare SED fitting.
        
        
        Returns
        -------
        Void
        """
        datain = open(filename).read().splitlines()
        id_, zspec, zphot = datain[1].split()
        self._id = int(id_)
        self._zspec = float(zspec)
        self._zphot = float(zphot)
        self._nfilters = int(datain[3].split()[1])
        
        # Model 
        col_model = datain[6].split()[1:]
        indexes_model = [i for i,d in enumerate(datain) if np.any([d.startswith(v) for v in ["GAL","STAR","QSO"]])]
        self._dfmodel = pandas.DataFrame([d.split() for d in np.asarray(datain)[indexes_model]], columns=col_model).set_index("Type")
        
        # Magnitudes
        col_mag = datain[2].split()[1:]
        self._dfmag = pandas.DataFrame([d.split() for i,d in enumerate(datain) if len(d.split())==len(col_mag)], columns=col_mag)
        
        i_base = len(self._dfmag)+len(self._dfmodel)
        istart = np.min([i_base+i for i,l in enumerate(datain[i_base:]) if len(l.split())==2])
        
        self.data = pandas.read_csv(filename, skiprows=istart, names=["lbda", "mag"], sep="  ", engine="python")
    
    # -------- #
    #  SETTER  #
    # -------- #    
    def set_lbda_range(self, lbdamin=None, lbdamax=None):
        """
        Set the spectrum wavelength range as attribute.
        
        Parameters
        ----------
        lbdamin : [float or None]
            Wavelength lower range limit.
            Default is None.
        
        lbdamax : [float or None]
            Wavelength upper range limit.
            Default is None.
        
        
        Returns
        -------
        Void
        """
        self._flagin = (self.data.lbda>lbdamin) & (self.data.lbda<lbdamax)
        
    # -------- #
    #  GETTER  #
    # -------- #        
    def get_data(self, model="gal", apply_lbdarange=True):
        """
        Return the spectrum (wavelengths, fluxes in erg/s/cm2/AA) as a pandas.DataFrame.
        
        Parameters
        ----------
        model : [string]
            Choice between the available LePhare model: gal, star, qso.
            If the SED fitting hasn't been done on the chosen model, the returned data will be empty.
            Default is "gal".
        
        Options
        -------
        apply_lbdarange : [bool]
            If True, apply the loaded wavelength range to the resulting spectrum.
            Default is True.
        
        
        Returns
        -------
        pandas.DataFrame
        """
        if model in ["gal"]:
            model = "GAL-1"
        index_start, index_end = self._get_model_nlines_(model)
        datas = self.data.iloc[index_start:index_end]
        return datas if not apply_lbdarange else datas[self.flagin[index_start:index_end]]

    def get_redshift(self, spectro=True):
        """
        Return the target redshift.
        
        Options
        -------
        spectro : [bool]
            If True, return the spectrometric redshift.
            Else, return the photometric one.
            Default is True.
        
        
        Returns
        -------
        float
        """
        return self.zspec if spectro else self.zphoto
    
    def get_spectral_data(self, model="gal", apply_lbdarange=True, influx=True, inhz=False, restframe=False):
        """
        Return the spectrum (wavelengths, fluxes or magnitudes) as a tuple.
        
        Parameters
        ----------
        model : [string]
            Choice between the available LePhare model: gal, star, qso.
            If the SED fitting hasn't been done on the chosen model, the returned data will be empty.
            Default is "gal".
        
        Options
        -------
        apply_lbdarange : [bool]
            If True, apply the loaded wavelength range to the resulting spectrum.
            Default is True.
        
        influx : [bool]
            If True, return the spectrum as flux ; False means as magnitude.
            Default is True.
        
        inhz : [bool]
            // Ignored if influx is False //
            Shall the fluxes be returned in erg/s/cm2/Hz (True) or in erg/s/cm2/AA (False).
            Default is False.
        
        restframe : [bool]
            If True, return the restframed spectrum.
            Default is False.
        
        
        Returns
        -------
        np.array, np.array
        """
        data_ = self.get_data(model=model, apply_lbdarange=apply_lbdarange)
        lbda_, mag_ = data_.lbda, data_.mag
        flux_aa = tools.mag_to_flux(mag_, None, wavelength=lbda_, inhz=False)[0]
        
        if restframe:
            lbda_, flux_aa = deredshift(lbda_, flux_aa, self.get_redshift(spectro=True))
            
        if influx:
            return lbda_, flux_aa if not inhz else tools.flux_aa_to_hz(flux_aa, lbda_)
        
        mag_ = tools.flux_to_mag(flux_aa, None, wavelength=lbda_, inhz=False)[0]
        return lbda_, mag_
    
    def get_input_data(self, influx=True, inhz=False):
        """
        Return the input data (used for the SED fitting resulting to this spectrum).
        
        Options
        -------
        influx : [bool]
            If True, return the spectrum as flux ; False means as magnitude.
            Default is True.
        
        inhz : [bool]
            // Ignored if influx = False //
            Shall the fluxes be returned in erg/s/cm2/Hz (True) or in erg/s/cm2/AA (False).
            Default is False.
        
        
        Returns
        -------
        np.array, [np.array, np.array]
        """
        lbda_, mag_, emag_ = np.asarray(self.resultmags[["Lbd_mean", "Mag", "emag"]].values, dtype="float").T
        if influx:
            return lbda_, tools.mag_to_flux(mag_, magerr=emag_, wavelength=lbda_, inhz=inhz)
        return lbda_, [mag_, emag_]
    
    def get_model_data(self, influx=True, inhz=False):
        """
        Return the model data (resulting from the SED fitting).
        
        Options
        -------
        influx : [bool]
            If True, return the spectrum as flux ; False means as magnitude.
            Default is True.
        
        inhz : [bool]
            // Ignored if influx = False //
            Shall the fluxes be returned in erg/s/cm2/Hz (True) or in erg/s/cm2/AA (False).
            Default is False.
        
        
        Returns
        -------
        np.array, np.array
        """
        lbda_, mag_ = np.asarray(self.resultmags[["Lbd_mean", "Mag_gal"]].values, dtype="float").T
        if influx:
            return lbda_, tools.mag_to_flux(mag_, None, wavelength=lbda_, inhz=inhz)[0]
        return lbda_, mag_
    
    def _get_model_nlines_(self, model):
        """
        Return the line range limiting the spectrum given a model.
        
        Parameters
        ----------
        model : [string]
            Choice between the available LePhare model: gal, star, qso.
            If the SED fitting hasn't been done on the chosen model, the returned range will be empty.
        
        
        Returns
        -------
        [int, int]
        """
        rangelist_ = [0]+list(np.cumsum(self.resultmodels["Nline"].astype("int")))
        rangelist = [[rl_,rm_] for rl_,rm_ in zip(rangelist_[:-1],rangelist_[1:])]
        return rangelist[self.resultmodels.index.get_loc(model)]

    # -------- #
    #  Main    #
    # -------- #
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
            If True, the spectrum is first deredshifted before doing the synthetic photometry.
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
        from sncosmo import bandpasses
        
        # - Get the corresponding bandpass
        if type(filter_) == str:
            bp = io.get_filter_bandpass(filter_)
        elif type(filter_) == bandpasses.Bandpass:
            bp = filter_
        elif len(filter_) == 2:
            bp = bandpasses.Bandpass(*filter_)
        else:
            raise TypeError("'filter_' must either be a filter name like filter_='sdss.u',"+
                            "a sncosmo.BandPass or a 2D array filter_=[wave, transmission].")

        # - Synthesize through bandpass
        sflux_aa = self.synthesize_photometry(bp.wave, bp.trans, restframe=restframe)
        slbda_ = bp.wave_eff/(1+self.get_redshift()) if restframe else bp.wave_eff

        if influx:
            return slbda_, (sflux_aa if not inhz else tools.flux_aa_to_hz(sflux_aa, slbda_))
        
        mag_ = tools.flux_to_mag(sflux_aa, None, wavelength=slbda_)[0]
        return slbda_, mag_
    
    def synthesize_through_filter(self, filtername,  restframe=False, influx=True, inhz=False):
        """
        Return photometry synthesized through the given filter/bandpass.
        The returned data are (effective wavelength, synthesize flux/mag) in an array with same size as for the number of draws.

        Parameters
        ----------
        filter_ : [string, sncosmo.BandPass, 2D array]
            The filter through which the spectrum will be synthesized.
            Must be a known filter (instrument.band), the actual bandpass will be grabbed using io.get_filter_bandpass(filter_).
        
        Options
        -------
        restframe : [bool]
            If True, the spectrum is first deredshifted before doing the synthetic photometry.
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
        return self.get_synthetic_photometry(filtername,  restframe=restframe, influx=influx, inhz=inhz)
    
    def synthesize_photometry(self, filter_lbda, filter_trans, model="gal", restframe=False):
        """
        Return the synthetic flux in AA.
        
        Parameters
        ----------
        filter_lbda, filter_trans : [array, array]
            Wavelength and transmission of the filter.
        
        Options
        -------
        model : [string]
            Choice between the available LePhare model: gal, star, qso.
            Default is "gal".
        
        restframe : [bool]
            If True, the spectrum is first deredshifted before doing the synthetic photometry.
            Default is False.
        
        
        Returns
        -------
        float
        """
        lbda_, flux_aa = self.get_spectral_data(model="gal", influx=True, inhz=False, restframe=restframe)
        return synthesize_photometry(lbda_, flux_aa, filter_lbda, filter_trans, normed=True)
    

    # -------- #
    # PLOTTER  #
    # -------- #        
    def show(self, ax=None, figsize=[7,3.5], ax_rect=[0.1,0.2,0.8,0.7],
             influx=True, inhz=False, model="gal",
             showdata=True, showmagmodel=True, restframe=False,
             scprop={}, colormodel=".7", set_label=True, **kwargs):
        """
        Plot the spectrum.
        
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
        
        influx : [bool]
            Shall the spectrum be in plotted in flux (True) or in magnitude (False).
            Default is True.
        
        inhz : [bool]
            // ignored if influx = False //
            Should the flux be returned in erg/s/cm2/Hz (True) or erg/s/cm2/AA (False).
            Default is False (thus in erg/s/cm2/AA).
        
        model : [string]
            Choice between the available LePhare model: gal, star, qso.
            Default is "gal".
        
        showdata : [bool]
            If True, overplot the input data.
            Default is True.
        
        showmagmodel : [bool]
            If True, overplot the model data.
        
        restframe : [bool]
            If True, the spectrum is first deredshifted before doing the synthetic photometry.
            Default is False.
        
        scprop : [dict]
            pyplot.scatter kwargs.
            Default is {}.
        
        colormodel : [pyplot.color]
            Color for the model data points.
            Default is "0.7".
        
        set_label : [bool]
            If True, set the axe labels.
            Default is True.
        
        **kwargs
        pyplot.plot kwargs
        
        
        Returns
        -------
        
        """
        import matplotlib.pyplot as mpl
        if ax is None:
            fig = mpl.figure(figsize=figsize)
            ax = fig.add_axes(ax_rect)
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
            
        return {"fig":fig, "ax":ax}
    
    # ============== #
    #   Properties   #
    # ============== #
    @property
    def zphoto(self):
        """ Fitted photometric redshift """
        return self._zphot

    @property
    def zspec(self):
        """ Spectrometric redshift """
        return self._zspec
    
    @property
    def mag(self):
        """ Magnitudes from spectrum ranged by 'flagin' """
        return self.data.mag[self.flagin]
    
    @property
    def flagin(self):
        """ Wavelength range """
        if not hasattr(self, "_flagin"):
            self.set_lbda_range(0,1e8)
        return self._flagin
    
    #
    # - Results
    #
    @property
    def resultmodels(self):
        """ Output results from LePhare SED fitting """
        return self._dfmodel
    
    @property
    def resultmags(self):
        """ Magnitudes """
        return self._dfmag







class LePhareSpectrumCollection( object ):
    """  """
    
    def __init__(self, lepharespectra):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        self._spectra = lepharespectra

    @classmethod
    def read_files(cls, filenames, lbda_range=None):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        return cls([LePhareSpectrum(file_, lbda_range) for file_ in filenames])
    
    # -------- #
    #  SETTER  #
    # -------- #
    def set_lbda_range(self, lbdamin=None, lbdamax=None):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        _ = [spec.set_lbda_range(lbdamin=None, lbdamax=None) for spec in self._spectra]

    # -------- #
    #  GETTER  #
    # -------- #
    def get_data(self, model="gal", apply_lbdarange=True):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        return [spec.get_data(model=model, apply_lbdarange=apply_lbdarange) for spec in self._spectra]

    def get_redshift(self, spectro=True):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        return [spec.get_redshift(spectro=spectro) for  spec in self._spectra]
    
    def get_spectral_data(self, model="gal", influx=True, inhz=False, restframe=False):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        return [spec.get_spectral_data(model=model, influx=influx, inhz=inhz, restframe=restframe) for  spec in self._spectra]
    
    def get_input_data(self, influx=True, inhz=False):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        return [spec.get_input_data(influx=influx, inhz=inhz) for  spec in self._spectra]
    
    def get_model_data(self, influx=True, inhz=False):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
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
        return [spec.get_synthetic_photometry(filter_, restframe=restframe, influx=influx, inhz=inhz) for  spec in self._spectra]
    
    def synthesize_photometry(self, filter_lbda, filter_trans, model="gal", restframe=False):
        """ get the synthetic flux in AA
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        return [spec.synthesize_photometry(filter_lbda, filter_trans, model=model, restframe=restframe) for  spec in self._spectra]

    # ============== #
    #   Properties   #
    # ============== #
    @property
    def spectra(self):
        """ """
        return self._spectra
