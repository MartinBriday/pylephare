
import numpy as np
import pandas
import time
import sys
import os

#Parameters initialization
from prospect.models import priors, SedModel
from prospect.models.templates import TemplateLibrary
from prospect.models.transforms import zfrac_to_masses, stellar_logzsol
from prospect.sources import CSPSpecBasis, FastStepBasis
from sedpy.observate import load_filters

#SED fitting
from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log
from dynesty.dynamicsampler import stopping_function, weight_function
from dynesty.utils import *


from ..k_correction import basesed
from propobject import BaseObject



class ProspectorSEDFitter( BaseObject ):
    """
    This class is a python SED fitter using Prospector.
    """
    
    PROPERTIES         = ["run_params", "obs", "sps", "model"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = ["mcmc_res"]
    
    RUN_PARAMS = {"verbose":True,
                  "debug":False,
                  "outfile":"test_snf",
                  "output_pickles": False,
                  "model_params": "parametric_sfh",
                  "mcmc": "dynesty",
                  # Optimization parameters
                  "do_powell": False,
                  "ftol":0.5e-5,
                  "maxfev": 5000,
                  "do_levenberg": True,
                  "nmin": 10,
                  # emcee fitting parameters
                  "nwalkers":128,
                  "nburn": [16, 32, 64],
                  "niter": 512,
                  "interval": 0.25,
                  "initial_disp": 0.1,
                  # dynesty Fitter parameters
                  "nested_bound": "multi", # bounding method
                  "nested_sample": "unif", # sampling method
                  "nested_nlive_init": 100,
                  "nested_nlive_batch": 100,
                  "nested_bootstrap": 0,
                  "nested_dlogz_init": 0.05,
                  "nested_weight_kwargs": {"pfrac": 1.0},
                  "nested_stop_kwargs": {"post_thresh": 0.1},
                  # Model parameters
                  "add_neb": False,
                  "add_dust": False,
                  # SPS parameters
                  "zcontinuous": 1,
                  # Fit parameters
                  "noise_model":False
                  }
    
    def __init__(self, **kwargs):
        """
        The class constructor can automatically execute 'set_data'.
            
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if kwargs != {}:
            self.set_data(**kwargs)

    def set_data(self, **kwargs):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        self.load_obs(**kwargs)
        self.load_sps(**kwargs)
        self.load_model(**kwargs)
        return
    
    def set_run_params(self, auto_add=False, **kwargs):
        """
        This method change parameter values.
        
        Parameters
        ----------
        kwargs : [dict]
            Dictionnary containing the parameters we want to change value as keys, and their new value as [dict] value.
        
        Options
        -------
        auto_add : [bool]
            If True and if the given parameters don't exist in the 'run_params' attribute, they will be created in 'run_params'.
        
        
        Returns
        -------
        Void
        """
        for key, value in kwargs.items():
            if key not in self.run_params:
                if auto_add:
                    self.add_run_params(**{key:value})
                else:
                    print("{} is not an existing parameter. If you really want to add it, execute add_run_params().".format(key))
            else:
                self.run_params[key] = value
    
    def add_run_params(self, **kwargs):
        """
        Create new parameters in the 'run_params' attribute.
        
        Parameters
        ----------
        kwargs : [dict]
            Dictionnary containing the new parameters we want to create.
            If the parameters already exist, pass.
        
        
        Returns
        -------
        Void
        """
        for key, value in kwargs.items():
            if key in self.run_params:
                print("{} is already existing. If you really want to change it, execute set_run_params().".format(key))
            else:
                self.run_params[key] = value
    
    def context_filters(self, context):
        """
        Return a list of the concerned filter bands relative to the given context.
        
        Parameters
        ----------
        context : [int]
            LePhare type context, it defines the used filter bands for the SED fitting.
        
        
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
    
    def load_obs(self, data_phot=None, data_spec=None, mask_spec=None, input_flux_unit="Hz", zspec=None, context_filters=None,
                 col_syntax={"mag":"mag_band", "mag.err":"mag_band_err", "lbda":"lbda", "flux":"flux", "flux.err":"flux.err"},
                 obs=None, **extras):
        """
        Prospector function.
        Set the 'obs' dictionnary (dealing with measurements) in a 'prospector' format :
        dict('list_bands', 'filter', 'zspec', 'maggies', 'maggies_unc', 'wavelength', 'spectrum', 'unc', 'mask').
        
        Parameters
        ----------
        data_phot : [dict or table or pandas.DataFrame or None]
            Table of the measured magnitudes and their error.
            Can be None if you want the algorithm to fit the spectrum instead.
            Or you can give both to fit both in the same time.
        
        data_spec : [dict or table or pandas.DataFrame or None]
            Table of the measured flux spectrum and its error.
            Can be None if you want the algorithm to fit the photometry instead.
            Or you can give both to fit both in the same time.
        
        zspec : [list or np.array or pandas.Series]
            List of the redshift we want to impose to each Type Ia Supernovae.
        
        Options
        -------
        mask_spec : [np.array or list or None]
            If a spectrum is given (through 'data_spec'), this mask will be applied to the spectrum.
            If None, the full spectrum is fitted.
        
        input_flux_unit : [string]
            Define the unit of the given spectrum flux :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
            The flux is converted in 'mgy' to be compatible with Prospector.
        
        context_filters : [int or None]
            If photometric measurements are given (through 'data_phot'), this will define the bands to take in account.
            Knowing the bands in 'data_phot' table, the number we want to set is : sum(2**[band_nb]).
            For example, bands = ["u", "g", "r", "i", "z"] (band_nb = [0, 1, 2, 3, 4]) :
            - context_filters = 31 --> ["u", "g", "r", "i", "z"]
            - context_filters = 30 --> ["g", "r", "i", "z"]
            - context_filters = 15 --> ["u", "g", "r", "i"]
            - context_filters = 25 --> ["u", "i", "z"]
            - etc.
        
        col_syntax : [dict]
            Define the syntax for your table columns concerning 'data_phot' and 'data_spec' input.
            Keys :
            - "mag", "mag.err" : put "band" where the filter is specified in the column names.
                                 (for example "mag_band" means your columns are "mag_u", "mag_g", etc. and "mag_band_err" gives "mag_u_err", etc.)
            - "flux", "flux.err" : flux and flux error syntax (don't forget "lbda").
            - "lbda" : wavelenght syntax (if the spectrum is given, so don't forget "flux" and "flux.err").
        
        obs : [dict or None]
            If not None, this will define the 'obs' attribute. Must be prospector compatible.
        
        
        Returns
        -------
        Void
        """
        if obs is not None:
            self._properties["obs"] = obs
            return
        
        self.obs["list_bands"] = self.context_filters(context_filters) if data_phot is not None else None
        self.obs["filters"] = load_filters([basesed.FILTER_BANDS[band]["prospector_name"] for band in self.obs["list_bands"]]) \
                              if data_phot is not None else None
        self.obs["zspec"] = np.asarray(zspec)

        ##### Photometric input #####
        if data_phot is not None:
            data_mag = {band:{"mag":data_phot[col_syntax["mag"].replace("band",band)],
                              "mag.err":data_phot[col_syntax["mag.err"].replace("band",band)]}
                        for band in self.obs["list_bands"]}
            data_phot = [np.asarray([basesed.mag_to_flux(data_mag[band]["mag"], data_mag[band]["mag.err"], band=band, flux_unit="mgy")[0]
                                     for band in self.obs["list_bands"]]),
                         np.asarray([basesed.mag_to_flux(data_mag[band]["mag"], data_mag[band]["mag.err"], band=band, flux_unit="mgy")[1]
                                     for band in self.obs["list_bands"]])]
        
        ##### Spectroscopic input #####
        elif data_spec is not None:
            data_spec = [data_spec[col_syntax["lbda"]],
                         basesed.convert_flux_unit((data_spec[col_syntax["flux"]], data_spec[col_syntax["flux.err"]]),
                                                   lbda=data_spec[col_syntax["lbda"]], unit_in=input_flux_unit, unit_out="mgy")]

        ##### Obs dictionnary #####
        self.obs["maggies"] = None if data_phot is None else data_phot[0]
        self.obs["maggies_unc"] = None if data_phot is None else data_phot[1]
        self.obs["wavelength"] = None if data_spec is None else data_spec[0]
        self.obs["spectrum"] = None if data_spec is None else data_spec[1][0]
        self.obs["unc"] = None if data_spec is None else data_spec[1][1]
        self.obs["mask"] = mask_spec
    
    def load_sps(self, sps=None, **extras):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if self.run_params["model_params"] == "parametric_sfh":
            self._properties["sps"] = CSPSpecBasis(zcontinuous=self.run_params["zcontinuous"]) if sps is None else sps
        else:
            self._properties["sps"] = FastStepBasis(zcontinuous=self.run_params["zcontinuous"]) if sps is None else sps
    
    def load_model(self, imposed_priors=None, **extras):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        model_params = TemplateLibrary[self.run_params["model_params"]]
        
        # Adjust model initial values (only important for emcee)
        model_params["dust2"]["init"] = 0.1
        model_params["logzsol"]["init"] = -0.3
        if self.run_params["model_params"] == "parametric_sfh":
            model_params["tage"]["init"] = 13.
            model_params["mass"]["init"] = 1e8
            model_params["tau"]["init"] = 1.
        
        # If we are going to be using emcee, it is useful to provide an
        # initial scale for the cloud of walkers (the default is 0.1)
        # For dynesty these can be skipped
        if self.run_params["model_params"] == "parametric_sfh":
            model_params["mass"]["init_disp"] = 1e7
            model_params["tau"]["init_disp"] = 3.0
            model_params["tage"]["init_disp"] = 5.0
            model_params["tage"]["disp_floor"] = 2.0
            model_params["dust2"]["disp_floor"] = 0.1
        
        #Priors
        model_params["logzsol"]["prior"] = priors.TopHat(mini=-2., maxi=2.)
        model_params["dust2"]["prior"] = priors.TopHat(mini=0., maxi=5.)
        if self.run_params["model_params"] == "parametric_sfh":
            model_params["tau"]["prior"] = priors.LogUniform(mini=1e-1, maxi=1e2)
            model_params["mass"]["prior"] = priors.LogUniform(mini=1e5, maxi=1e11)
            model_params["tage"]["prior"] = priors.TopHat(mini=0.001, maxi=20.8)
        if "logmass" in model_params.keys():
            model_params["logmass"]["prior"] = priors.TopHat(mini=5., maxi=12.)
        
        if imposed_priors is not None:
            for k, v in imposed_priors.items():
                if k == "dust2":
                    model_params[k]["prior"] = priors.ClippedNormal(mean=v["mean"], sigma=v["sigma"], mini=0., maxi=10.)
                elif k == "logmass":
                    model_params[k]["prior"] = priors.ClippedNormal(mean=v["mean"], sigma=v["sigma"], mini=0., maxi=20.)
                else:
                    try:
                        model_params[k]["prior"] = priors.Normal(mean=v["mean"], sigma=v["sigma"])
                    except:
                        print("'{}' is not a valid model parameter.".format(k))
                        continue
        
        # make sure zred is fixed
        model_params["zred"]["isfree"] = False
        # And set the value to the object_redshift keyword
        model_params["zred"]["init"] = self.obs["zspec"]
    
        if self.run_params["add_dust"]:
            # Add dust emission (with fixed dust SED parameters)
            model_params.update(TemplateLibrary["dust_emission"])
    
        if self.run_params["add_neb"]:
            # Add nebular emission (with fixed parameters)
            model_params.update(TemplateLibrary["nebular"])
    
        # Now instantiate the model using this new dictionary of parameter specifications
        self._properties["model"] = SedModel(model_params)
    
    def lnprobfn(self, theta, model=None, obs=None, sps=None, noise=None, nested=True, residuals=False, verbose=False):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        model = self.model if model is None else model
        obs = self.obs if obs is None else obs
        
        lnp_prior = model.prior_product(theta, nested=nested)
        if not np.isfinite(lnp_prior):
            return -np.infty
        
        # Generate mean model
        t1 = time.time()
        try:
            spec, phot, x = model.mean_model(theta, obs, sps=self.sps)
        except(ValueError):
            return -np.infty
        d1 = time.time() - t1
        
        # Noise modeling
        if self.run_params["noise_model"]:
            spec_noise, phot_noise = model_setup.load_gp(**self.run_params)
            if spec_noise is not None:
                spec_noise.update(**model.params)
            if phot_noise is not None:
                phot_noise.update(**model.params)
        else:
            spec_noise, phot_noise = None, None

        vectors = {'spec': spec,
                   'unc':  obs['unc'],
                   'sed':  model._spec,
                   'cal':  model._speccal,
                   'phot': phot,
                   'maggies_unc': obs['maggies_unc']}

        # Calculate likelihoods
        t2 = time.time()
        lnp_spec = lnlike_spec(spec, obs=obs, spec_noise=spec_noise, **vectors)
        lnp_phot = lnlike_phot(phot, obs=obs, phot_noise=phot_noise, **vectors)
        d2 = time.time() - t2
        
        if verbose:#self.run_params["verbose"]:
            write_log(theta, lnp_prior, lnp_spec, lnp_phot, d1, d2)

        return lnp_phot + lnp_spec + lnp_prior
    
    def prior_transform(self, u):
        #model = self.model if model is None else model
        return self.model.prior_transform(u)
    
    def halt(message, pool=None):
        """
        Exit, closing pool safely.
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        print(message)
        try:
            pool.close()
        except:
            pass
        sys.exit(0)

    def mcmc_dynesty(self, pool, nprocs):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if self.run_params['verbose']:
            print('dynesty sampling...')
        tstart = time.time()  # time it
        dynestyout = fitting.run_dynesty_sampler(self.lnprobfn, self.prior_transform, self.model.ndim,
                                                 pool=pool, queue_size=nprocs,
                                                 stop_function=stopping_function,
                                                 wt_function=weight_function,
                                                 **self.run_params)
        ndur = time.time() - tstart
        print('done dynesty in {0}s'.format(ndur))
        return {"dynestyout":dynestyout, "ndur":ndur}
    
    def init_guess(self, initial_theta, pool):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        from prospect.fitting.fitting import run_minimize
        
        if self.run_params['verbose']:
            print('Starting minimization...')
        
        if not np.isfinite(self.model.prior_product(self.model.initial_theta.copy())):
            halt("Halting: initial parameter position has zero prior probability.")
        
        nmin = self.run_params.get('nmin', 1)
        if pool is not None:
            nmin = max([nmin, pool.size])
        
        if bool(self.run_params.get('do_powell', False)):
            powell_opt = {'ftol': self.run_params['ftol'], 'xtol': 1e-6, 'maxfev': self.run_params['maxfev']}
            guesses, pdur, best = run_minimize(self.obs, self.model, self.sps, noise=None, lnprobfn=self.lnprobfn,
                                               min_method='powell', min_opts={"options": powell_opt},
                                               nmin=nmin, pool=pool)
            initial_center = fitting.reinitialize(guesses[best].x, self.model, edge_trunc=self.run_params.get('edge_trunc', 0.1))
            initial_prob = -guesses[best]['fun']
            if self.run_params['verbose']:
                print('done Powell in {0}s'.format(pdur))
                print('best Powell guess:{0}'.format(initial_center))
        
        elif bool(self.run_params.get('do_levenberg', False)):
            lm_opt = {"xtol": 5e-16, "ftol": 5e-16}
            guesses, pdur, best = run_minimize(self.obs, self.model, self.sps, noise=None, lnprobfn=self.lnprobfn,
                                               min_method='lm', min_opts=lm_opt,
                                               nmin=nmin, pool=pool)
            initial_center = fitting.reinitialize(guesses[best].x, self.model, edge_trunc=self.run_params.get('edge_trunc', 0.1))
            initial_prob = None
            if self.run_params['verbose']:
                print('done L-M in {0}s'.format(pdur))
                print('best L-M guess:{0}'.format(initial_center))
        
        else:
            if self.run_params['verbose']:
                print('No minimization requested.')
            guesses = None
            pdur = 0.0
            initial_center = initial_theta.copy()
            initial_prob = None

        return {"guesses":guesses, "pdur":pdur, "initial_center":initial_center, "initial_prob":initial_prob}

    def mcmc_emcee(self, pool):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        postkwargs = {}
        initial_theta = self.model.rectify_theta(self.model.initial_theta)
        
        try:
            import h5py
            hfile = h5py.File(self.hfilename, "a")
            print("Writing to file {}".format(self.hfilename))
            write_results.write_h5_header(hfile, self.run_params, self.model)
            write_results.write_obs_to_h5(hfile, self.obs)
        except(ImportError):
            hfile = None
        
        init_guess = self.init_guess(initial_theta, pool)
        
        if self.run_params['verbose']:
            print('emcee sampling...')
        tstart = time.time()
        out = fitting.run_emcee_sampler(self.lnprobfn, init_guess["initial_center"], self.model,
                                        postkwargs=postkwargs, prob0=init_guess["initial_prob"],
                                        pool=pool, hdf5=hfile, **self.run_params)
        esampler, burn_p0, burn_prob0 = out
        ndur = time.time() - tstart
        if self.run_params['verbose']:
            print('done emcee in {0}s'.format(ndur))
        
        return {"hfile":hfile, "esampler":esampler, "burn_p0":burn_p0, "burn_prob0":burn_prob0, "init_guess":init_guess, "ndur":ndur}
    
    def run_fit(self, pool=None, nprocs=1, write_res=False):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        try:
            self.run_params['sps_libraries'] = self.sps.ssp.libraries
        except(AttributeError):
            self.run_params['sps_libraries'] = None
        
        if self.run_params.get("debug", False):
            self.halt('stopping for debug')
        
        # Try to set up an HDF5 file and write basic info to it
        odir = os.path.dirname(os.path.abspath(self.run_params['outfile']))
        if (not os.path.exists(odir)):
            badout = 'Target output directory {} does not exist, please make it.'.format(odir)
            halt(badout)
        
        # -------
        # Sample
        # -------
        if self.run_params["mcmc"] == "dynesty":
            mcmc_res = self.mcmc_dynesty(pool, nprocs)
        elif self.run_params["mcmc"] == "emcee":
            mcmc_res = self.mcmc_emcee(pool)
        
        self._derived_properties["mcmc_res"] = mcmc_res

        if write_res:
            self.write_results()

    def write_results(self):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        # -------------------------
        # Output HDF5 (and pickles if asked for)
        # -------------------------
        if self.run_params.get("output_pickles", False):
            # Write the dynesty result object as a pickle
            import pickle
            with open(self.run_params['outfile'] + '_dns.pkl', 'w') as f:
                pickle.dump(self.dynestyout, f)
        
            # Write the model as a pickle
            partext = write_results.paramfile_string(**self.run_params)
            write_results.write_model_pickle(self.run_params['outfile'] + '_model', self.model, powell=None, paramfile_text=partext)
    
        # Write HDF5
        if self.run_params["mcmc"] == "dynesty":
            write_results.write_hdf5(self.hfilename, self.run_params, self.model, self.obs, self.mcmc_res["dynestyout"], None, tsample=self.mcmc_res["ndur"])
        elif self.run_params["mcmc"] == "emcee":
            write_results.write_hdf5(self.mcmc_res["hfile"], self.run_params, self.model, self.obs, self.mcmc_res["esampler"], None, #self.mcmc_res["init_guess"]["guesses"],
                                     toptimize=self.mcmc_res["init_guess"]["pdur"], tsample=self.mcmc_res["ndur"],
                                     sampling_initial_center=self.mcmc_res["init_guess"]["initial_center"],
                                     post_burnin_center=self.mcmc_res["burn_p0"], post_burnin_prob=self.mcmc_res["burn_prob0"])
    
            
    





    #-------------------#
    #   Properties      #
    #-------------------#
    @property
    def run_params(self):
        """  """
        if self._properties["run_params"] is None:
            self._properties["run_params"] = self.RUN_PARAMS
        return self._properties["run_params"]

    @property
    def obs(self):
        """  """
        if self._properties["obs"] is None:
            self._properties["obs"] = {}
        return self._properties["obs"]

    @property
    def sps(self):
        """  """
        return self._properties["sps"]

    @property
    def model(self):
        """  """
        return self._properties["model"]

    @property
    def mcmc_res(self):
        """  """
        return self._derived_properties["mcmc_res"]

    @property
    def hfilename(self):
        """  """
        return self.run_params["outfile"] + "_mcmc_" + self.run_params["mcmc"] + ".h5"

