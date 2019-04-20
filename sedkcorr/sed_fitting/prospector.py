
import numpy as np
import pandas
import time
import sys
import os

#Parameters initialization
from prospect.models import priors, SedModel
from prospect.models.templates import TemplateLibrary
from prospect.sources import CSPSpecBasis
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


# --------------#
#   RUN_PARAMS  #
# --------------#
RUN_PARAMS = {'verbose':True,
              'debug':False,
              'outfile':'test_snf',
              'output_pickles': False,
              # dynesty Fitter parameters
              'nested_bound': 'multi', # bounding method
              'nested_sample': 'unif', # sampling method
              'nested_nlive_init': 100,
              'nested_nlive_batch': 100,
              'nested_bootstrap': 0,
              'nested_dlogz_init': 0.05,
              'nested_weight_kwargs': {"pfrac": 1.0},
              'nested_stop_kwargs': {"post_thresh": 0.1},
              # Model parameters
              'add_neb': False,
              'add_dust': False,
              # SPS parameters
              'zcontinuous': 1,
              # Fit parameters
              'noise_model':False
              }



class ProspectorSEDFitter( BaseObject ):
    """
    
    """
    
    PROPERTIES         = ["run_params", "obs", "sps", "model"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = ["dynestyout"]
    
    def __init__(self, **kwargs):
        """
        
        """
        if kwargs != {}:
            self.set_data(**kwargs)

    def set_data(self, **kwargs):
        """
        
        """
        self.load_obs(**kwargs)
        self.load_sps(**kwargs)
        self.load_model(**kwargs)
        return
    
    def set_run_params(self, **kwargs):
        """
        
        """
        for key, value in kwargs.items():
            if key not in self.run_params:
                print("{} is not an existing parameter. If you really want to add it, execute add_run_params().".format(key))
                continue
            self.run_params[key] = value
    
    def add_run_params(self, **kwargs):
        """
        
        """
        for key, value in kwargs.items():
            if key not in self.run_params:
                print("{} is already existing. If you really want to change it, execute set_run_params().".format(key))
                continue
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
    
    def load_obs(self, data=None, col_syntax=["mag_band", "mag_band_err"], **extras):
        """
        
        """
        self.obs["list_bands"] = self.context_filters(data["CONTEXT"])
        self.obs["filters"] = load_filters([basesed.FILTER_BANDS[band]["prospector_name"] for band in self.obs["list_bands"]])
        
        data_mag = {band:{"mag":data[col_syntax[0].replace("band",band)],
                          "mag.err":data[col_syntax[1].replace("band",band)]}
                    for band in self.obs["list_bands"]}
        
        self.obs["maggies"] = np.asarray([basesed.flux_nu_to_mgy(basesed.band_mag_to_flux(data_mag[band]["mag"], band))
                                          for band in self.obs["list_bands"]])
        self.obs["maggies_unc"] = np.asarray([basesed.flux_nu_to_mgy(basesed.band_mag_to_flux_err(data_mag[band]["mag"], data_mag[band]["mag.err"], band))
                                              for band in self.obs["list_bands"]])
        
        self.obs["zspec"] = data["Z-SPEC"]
        
        self.obs["wavelength"] = None
        self.obs["spectrum"] = None
        self.obs["unc"] = None
        self.obs["mask"] = None
    
    def load_sps(self, zcontinuous=1, sps=None, **extras):
        """
        
        """
        self._properties["sps"] = CSPSpecBasis(zcontinuous=zcontinuous) if sps is None else sps
    
    def load_model(self, add_dust=False, add_neb=False, **extras):
        """
        
        """
        model_params = TemplateLibrary["parametric_sfh"]
        
        # make sure zred is fixed
        model_params["zred"]["isfree"] = False
        # And set the value to the object_redshift keyword
        model_params["zred"]["init"] = self.obs["zspec"]
    
        if add_dust:
            # Add dust emission (with fixed dust SED parameters)
            model_params.update(TemplateLibrary["dust_emission"])
    
        if add_neb:
            # Add nebular emission (with fixed parameters)
            model_params.update(TemplateLibrary["nebular"])
    
        # Now instantiate the model using this new dictionary of parameter specifications
        self._properties["model"] = SedModel(model_params)
    
    def lnprobfn(self, theta, model=None, obs=None, nested=True, verbose=False):
        """
        
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
        """
        print(message)
        try:
            pool.close()
        except:
            pass
        sys.exit(0)
    
    def run_fit(self, pool=None, nprocs=1, write_res=False):
        """
        
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
        self._derived_properties["dynestyout"] = dynestyout

        if write_res:
            self.write_results(ndur)

    def write_results(self, ndur=0.):
        """
        
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
        hfile = self.run_params['outfile'] + '_mcmc.h5'
        write_results.write_hdf5(hfile, self.run_params, self.model, self.obs, self.dynestyout, None, tsample=ndur)
    
            
    





    #-------------------#
    #   Properties      #
    #-------------------#
    @property
    def run_params(self):
        """  """
        if self._properties["run_params"] is None:
            self._properties["run_params"] = RUN_PARAMS
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
    def dynestyout(self):
        """  """
        return self._derived_properties["dynestyout"]

