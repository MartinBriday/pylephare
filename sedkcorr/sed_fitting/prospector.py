
import numpy as np
import pandas

from prospect.models import priors, SedModel
from prospect.models.templates import TemplateLibrary
from prospect.sources import CSPSpecBasis
from sedpy.observate import load_filters

from ..k_correction import basesed
from propobject import BaseObject


# --------------#
#   RUN_PARAMS  #
# --------------#
RUN_PARAMS = {'verbose':True,
              'debug':False,
              'outfile':'test_galphot',
              'output_pickles': False,
              # Optimization parameters
              'do_powell': False,
              'ftol':0.5e-5, 'maxfev': 5000,
              'do_levenberg': True,
              'nmin': 10,
              # emcee fitting parameters
              'nwalkers':128,
              'nburn': [16, 32, 64],
              'niter': 512,
              'interval': 0.25,
              'initial_disp': 0.1,
              # dynesty Fitter parameters
              'nested_bound': 'multi', # bounding method
              'nested_sample': 'unif', # sampling method
              'nested_nlive_init': 100,
              'nested_nlive_batch': 100,
              'nested_bootstrap': 0,
              'nested_dlogz_init': 0.05,
              'nested_weight_kwargs': {"pfrac": 1.0},
              'nested_stop_kwargs': {"post_thresh": 0.1},
              # Obs data parameters
              'objid':0,
              #'phottable': 'demo_photometry.dat',
              'logify_spectrum':False,
              'normalize_spectrum':False,
              'luminosity_distance': 1e-5,  # in Mpc
              # Model parameters
              'add_neb': False,
              'add_dust': False,
              # SPS parameters
              'zcontinuous': 1,
              }



class ProspectorSEDFitter( BaseObject ):
    """
    
    """
    
    PROPERTIES         = ["run_params", "obs", "sps", "model"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
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
        for ii in range(len(basesed.FILTER_BANDS)-1,-1,-1):
            if (context - 2**ii) >= 0:
                context = context - 2**ii
                idx.append(ii)
        return [band for band in basesed.LIST_BANDS if FILTER_BANDS[band]["context_id"] in idx]
    
    def load_obs(self, objid=0, data_photo=None, col_syntax=["mag_band", "mag_band_err"], **kwargs):
        """
        
        """
        list_bands = context_filters(data_photo["CONTEXT"])
        self.obs["filters"] = load_filters([FILTER_BANDS[band]["prospector_name"] for band in list_bands])
        
        data_mag = {band:{"mag":data_photo[col_syntax[0].replace("band",band)],
                          "mag.err":data_photo[col_syntax[1].replace("band",band)]}
                    for band in list_bands}
        
        self.obs["maggies"] = [basesed.flux_nu_to_mgy(basesed.band_mag_to_flux(data_mag[band]["mag"], band))
                               for band in list_bands]
        self.obs["maggies_unc"] = [basesed.flux_nu_to_mgy_err(basesed.band_mag_to_flux(data_mag[band]["mag"], data_mag[band]["mag.err"], band))
                                   for band in list_bands]
        
        self.obs["wavelength"] = None
        self.obs["spectrum"] = None
        self.obs["unc"] = None
        self.obs["mask"] = None
    
    def load_sps(self, zcontinuous=1, **extras):
        """
        
        """
        self._properties["sps"] = CSPSpecBasis(zcontinuous=zcontinuous)
    
    def load_model(self, zspec=None, add_dust=False, add_neb=False, **extras):
        """
        
        """
        model_params = TemplateLibrary["parametric_sfh"]
        
        if zspec is not None:
            # make sure zred is fixed
            model_params["zred"]["isfree"] = False
            # And set the value to the object_redshift keyword
            model_params["zred"]["init"] = zspec
    
        if add_dust:
            # Add dust emission (with fixed dust SED parameters)
            model_params.update(TemplateLibrary["dust_emission"])
    
        if add_neb:
            # Add nebular emission (with fixed parameters)
            model_params.update(TemplateLibrary["nebular"])
    
        # Now instantiate the model using this new dictionary of parameter specifications
        self._properties["model"] = SedModel(model_params)
    
    
    





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
        return self.properties["obs"]

    @property
    def sps(self):
        """  """
        return self._properties["sps"]

    @property
    def model(self):
        """  """
        return self._properties["model"]
