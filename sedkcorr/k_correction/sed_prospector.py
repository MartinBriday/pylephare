
import numpy as np
import pandas

from prospect.io.read_results import results_from, get_sps
from prospect.io.read_results import traceplot, subcorner
from prospect.utils.plotting import quantile


from ..sed_fitting import prospector
from . import basesed




class SED_prospector( basesed.SED ):
    """
    This class read prospector SED fits and apply k-correction on it.
    Child of the class 'basesed.SED'.
    """
    
    PROPERTIES         = ["p_res", "p_obs", "p_mod", "p_sps"]
    SIDE_PROPERTIES    = ["p_run_params"]
    DERIVED_PROPERTIES = ["post_pcts", "sed_stack", "phot_stack"]

    def read_fit_results(self, filename=None, sps=None, **extras):
        """
        Read the prospector output file and initiate the specicific parameters (obs, model, sps).
        
        Parameters
        ----------
        filename : [string]
            Path of the prospector output file.
        
        Options
        -------
        sps : [prospector.sps]
            Give a sps to earn the loading time.
        
        
        Returns
        -------
        Void
        """
        res, obs, mod = results_from(filename, dangerous=False)
        self._properties["p_res"] = res
        self._properties["p_obs"] = obs
        
        buf = prospector.ProspectorSEDFitter()
        buf.set_run_params(auto_add=True, **self.p_run_params)
        buf.load_obs(obs=self.p_obs)
        buf.load_model()
        buf.load_sps(sps)
        self._properties["p_mod"] = buf.model
        self._properties["p_sps"] = buf.sps
    
    def set_data_sed(self, filename=None, fit_errors=True, **extras):
        """
        Extract the fitted SED from the prospector output file and initiate SED data using the parent class (basesed.SED) 'set_data_sed' method.
        
        Parameters
        ----------
        ilename : [string]
            Path of the prospector output file.
        
        Options
        -------
        fit_errors : [bool]
            If True, generate a stack of spectrum to calcultate SED errors.
        
        # **extras #
        
        ### read_fit_results ###
        sps : [prospector.sps]
            Give a sps to earn the loading time.
        
        ### get_sed_error ###
        nb_walkers_points : [int]
            Number of spectrum to generate in the stack using 'set_sed_stack' method.
        
        
        Returns
        -------
        Void
        """
        self.read_fit_results(filename, **extras)
        
        imax = np.argmax(self.p_res["lnprobability"])
        if self.p_run_params["mcmc"] == "emcee":
            i, j = np.unravel_index(imax, self.p_res["lnprobability"].shape)
            theta_max = self.p_res["chain"][i, j, :].copy()
        else:
            theta_max = self.p_res["chain"][imax, :]
        
        if fit_errors:
            self.set_sed_stack(**extras)
        data_sed = {"lbda":self.get_sed_wavelength(),
                    "flux":self.get_sed_flux(theta_max),
                    "flux.err":np.zeros(self.nb_spec_points) if not fit_errors else self.get_sed_error(**extras)}
        _ = super(SED_prospector, self).set_data_sed( pandas.DataFrame(data_sed) )
    
    def context_filters(self, context):
        """
        Return a list of the concerned filter bands relative to the given context.
        
        Parameters
        ----------
        context : [int]
            LePhare type context, it defines the used filter bands for the SED fitting.
            Knowing the bands in 'self.list_bands', the number we want to set is : sum(2**[band_nb]).
            For example, bands = ["u", "g", "r", "i", "z"] (band_nb = [0, 1, 2, 3, 4]) :
            - context = 31 --> ["u", "g", "r", "i", "z"]
            - context = 30 --> ["g", "r", "i", "z"]
            - context = 15 --> ["u", "g", "r", "i"]
            - context = 25 --> ["u", "i", "z"]
            - etc.
        
        
        Returns
        -------
        list(string)
        """
        idx = []
        for ii in range(len(basesed.LIST_BANDS)-1,-1,-1):
            if (context - 2**ii) >= 0:
                context = context - 2**ii
                idx.append(ii)
        return [band for band in basesed.LIST_BANDS if basesed.FILTER_BANDS[band]["context_id"] in idx]
    
    def set_data_meas(self, data_meas=None, z=None, col_syntax=["mag_band", "mag_band_err"], list_bands=None, **extras):
        """
        Set the host redshift and the measured magnitudes for every filter bands used in SED fitting.
        
        Parameters
        ----------
        data_meas : [table like]
            Table like (eg : DataFrame line) of the measurements.
        
        z : [float or None]
            Redshift of the SNeIa host.
            If None, the redshift is supposed to be in the data_meas table under the name "Z-SPEC".
            Default is None.
        
        col_syntax : [list[string]]
            Syntax of measurements and errors column names in the data_meas table.
            Replace the filter band in the column names with the word "band" (eg: ["mag_band", "mag_band_err"]).
        
        Options
        -------
        list_bands : [list[string] or None]
            List of the filter bands used in SED fitting.
            If None, the LePhare context set list_bands. The context is supposed to be in the data_meas table under the name "CONTEXT".
        
        
        Returns
        -------
        Void
        """
        
        ############# Add an option to read phot data from prospector results ###############
        
        self._side_properties["list_bands"] = self.context_filters(data_meas["CONTEXT"]) if list_bands is None else list_bands
        z = z if z is not None else data_meas["Z-SPEC"]
        
        data_meas = {band:{"mag":data_meas[col_syntax[0].replace("band",band)],
                           "mag.err":data_meas[col_syntax[1].replace("band",band)]}
                     for band in self.list_bands}
        _ = super(SED_prospector, self).set_data_meas(data_meas=data_meas, z=z)
    
    def get_sed_flux(self, theta):
        """
        Return a complete spectrum, given by a set of fitted parameters.
        Output unit is erg/cm2/s/Hz.
        
        Parameters
        ----------
        theta : [list(float)]
            List of parameter values fitted in prospector.
        
        
        Returns
        -------
        list(float)
        """
        mspec, mphot, mextra = self.p_mod.mean_model(theta, self.p_obs, sps=self.p_sps)
        mspec = basesed.convert_flux_unit(mspec, lbda=self.get_sed_wavelength(), unit_in="mgy", unit_out="Hz")
        return mspec
    
    def get_sed_wavelength(self):
        """
        Recover the wavelength table from the prospector output file.
        If no 'wavelength' was given in the 'obs' dict during prospector fitting, it recovers the fitted one and redshifts it (as it is the restframe).
        
        
        Returns
        -------
        list(float)
        """
        if self.p_obs["wavelength"] is None:
            # *restframe* spectral wavelengths, since obs["wavelength"] is None
            a = 1.0 + self.p_obs.get("zspec", 0.0)
            wspec = self.p_sps.wavelengths.copy()
            wspec *= a #redshift them
        else:
            wspec = self.p_obs["wavelength"]
        return wspec
    
    def set_post_pcts(self, weights=False, chain_frac=None):
        """
        Set up as attribute the quantiles on fitted parameters.
        
        Options
        -------
        weights : [bool]
            If True, take in account the weights on the chains, given by prospector output file.
            Default is False (same weight for the whole sample).
        
        chain_frac : [float or None]
            Between 0 and 1, fraction of the chains to calculate the quantiles, starting from the end.
        
        
        Returns
        -------
        Void
        """
        chain_iter = 0 if chain_frac is None else -(int(len(self.p_res["chain"])*chain_frac))
        post_pcts = [quantile(self.p_res["chain"][chain_iter:, ii], percents=[16, 50, 84],
                              weights=(self.p_res["weights"][chain_iter:] if weights else None))
                     for ii, param in enumerate(self.p_res["theta_labels"])]
        self._derived_properties["post_pcts"] = post_pcts
    
    def set_sed_stack(self, nb_walkers_points=500, opt_random=True, **extras):
        """
        Set up a stack of spectrum in order to get the SED and the fitted photometry errors.
        
        Parameters
        ----------
        nb_walkers_points : [int]
            Number of spectrum to generate in the stack.
        
        Options
        -------
        opt_random : [bool]
            If True, take *nb_walkers_points* randomly chosen fitted parameters in their chains to generate the spectra.
            Else take the last *nb_walkers_points* fitted parameters from the chains.
        
        
        Returns
        -------
        Void
        """
        randint = np.random.randint
        theta = []
        if self.p_run_params["mcmc"] == "emcee":
            rand_chain = randint(self.p_run_params["nwalkers"], size=nb_walkers_points)
            rand_iter = randint(self.p_run_params["niter"], size=nb_walkers_points)
            theta = self.p_res["chain"][rand_chain, rand_iter, :]
        elif self.p_run_params["mcmc"] == "dynesty":
            rand_iter = randint(len(self.p_res["chain"]), size=nb_walkers_points)
            theta = self.p_res["chain"][rand_iter, :] if opt_random else self.p_res["chain"][-nb_walkers_points:, :]

        mspec = np.empty((nb_walkers_points, self.nb_spec_points))
        mphot = np.empty((nb_walkers_points, self.nb_phot_points))
        for ii in np.arange(nb_walkers_points):
            mspec[ii], mphot[ii], _ = self.p_mod.mean_model(theta[ii], self.p_obs, sps=self.p_sps)
        self._derived_properties["sed_stack"] = mspec
        self._derived_properties["phot_stack"] = mphot
        
    def get_sed_error(self, nb_walkers_points=None, **extras):
        """
        Return the SED error, calculated from the prospector fit result.
        
        Parameters
        ----------
        nb_walkers_points : [int]
            Number of spectrum to generate in the stack using 'set_sed_stack' method.
        
        
        Returns
        -------
        np.array
        """
        if nb_walkers_points is not None and len(self.sed_stack) != nb_walkers_points:
            self.set_sed_stack(nb_walkers_points)
        buf_q = np.quantile(self.sed_stack, [0.16, 0.84], axis=0)
        return basesed.convert_flux_unit((buf_q[1]-buf_q[0])/2, lbda=self.get_sed_wavelength(), unit_in="mgy", unit_out="Hz")
    
    def get_kcorr_error(self, nb_walkers_points=None, **extras):
        """
        Return k-correction error, calculated from the prospector fit result.
        
        Parameters
        ----------
        nb_walkers_points : [int]
            Number of spectrum to generate in the stack using 'set_sed_stack' method.
        
        
        Returns
        -------
        dict
        """
        if nb_walkers_points is not None and len(self.sed_stack) != nb_walkers_points:
            self.set_sed_stack(nb_walkers_points)
        kcorr_stack = []
        data_lbda = self.get_sed_wavelength()
        for data_sed in self.sed_stack:
            data_sed = basesed.convert_flux_unit(data_sed, lbda=data_lbda, unit_in="mgy", unit_out="Hz")
            data_sed = {"lbda":data_lbda, "flux":data_sed}
            data_kcorr = self.get_phot(data_sed=data_sed, bands=None)
            kcorr_stack.append([data_kcorr[band] for band in basesed.LIST_BANDS])
        kcorr_errors = {band:np.quantile(kcorr_stack, [0.16, 0.84], axis=0).T[ii] for ii, band in enumerate(basesed.LIST_BANDS)}
        for band in basesed.LIST_BANDS:
            kcorr_errors[band] = (kcorr_errors[band][1]-kcorr_errors[band][0])/2
        return kcorr_errors
    
    def k_correction(self, **extras):
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
        _ = super(SED_prospector, self).k_correction(kcorr_flux_error=self.get_kcorr_error(**extras))

    def show_walkers(self, figsize=(9,7), **kwargs):
        """
        Run 'traceplot' method of 'prospect.io.read_results' with the results loaded with 'set_data_sed'.
        
        Options
        -------
        figsize : [tuple]
            Size of the figure (matplotlib.pyplot option).
        
        showpars : [list(string)]
            List of the fitted parameters we want to see the chains.
        
        
        Returns
        -------
        fig
        """
        tracefig = traceplot(self.p_res, figsize=figsize, **kwargs)
        return tracefig

    
    





    #-------------------#
    #   Properties      #
    #-------------------#
    @property
    def p_res(self):
        """ prospector 'res' """
        return self._properties["p_res"]

    @property
    def p_obs(self):
        """ prospector 'obs' """
        return self._properties["p_obs"]

    @property
    def p_mod(self):
        """ prospector 'model' """
        return self._properties["p_mod"]
    
    @property
    def p_sps(self):
        """ prospector 'sps' """
        return self._properties["p_sps"]

    @property
    def p_run_params(self):
        """ prospector 'run_params' dictionnary """
        if self._side_properties["p_run_params"] is None:
            self._side_properties["p_run_params"] = self.p_res["run_params"]
        return self._side_properties["p_run_params"]
            
    @property
    def post_pcts(self):
        """ Quantiles of the fitted parameters """
        if self._derived_properties["post_pcts"] is None:
            self.set_post_pcts()
        return self._derived_properties["post_pcts"]

    @property
    def sed_stack(self):
        """ List of SED spectrum """
        return self._derived_properties["sed_stack"]
    
    @property
    def phot_stack(self):
        """ List of fitted photometry """
        return self._derived_properties["phot_stack"]
    
    @property
    def nb_spec_points(self):
        """ Length of the spectrum """
        return len(self.get_sed_wavelength())
    
    @property
    def nb_phot_points(self):
        """ Number of photometric points """
        return len(self.p_obs["maggies"])
