
import numpy as np
from prospect.models import priors, SedModel
from prospect.models.templates import TemplateLibrary
from prospect.sources import CSPSpecBasis
from sedpy.observate import load_filters

from .. import basesed
from propobject import BaseObject


# --------------#
#   RUN_PARAMS  #
# --------------#
RUN_PARAMS = {'verbose':True,
              'debug':False,
              'outfile':'demo_galphot',
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
              'phottable': 'demo_photometry.dat',
              'logify_spectrum':False,
              'normalize_spectrum':False,
              'luminosity_distance': 1e-5,  # in Mpc
              # Model parameters
              'add_neb': False,
              'add_dust': False,
              # SPS parameters
              'zcontinuous': 1,
              }



class ProspectorParameters( BaseObject ):
    """
    
    """
    
    PROPERTIES         = ["run_params"]
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
                    #return [basesed.FILTER_BANDS[ii] for ii in reversed(idx)]
        return [k for k, v in basesed.FILTER_BANDS.items() if v["context_id"] in idx]





    #-------------------#
    #   Properties      #
    #-------------------#
    @property
    def run_params(self):
        """  """
        if self._properties["run_params"] is None:
            self._properties["run_params"] = RUN_PARAMS
        return self._properties["run_params"]

