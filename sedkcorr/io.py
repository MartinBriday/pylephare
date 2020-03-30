""" Module made to manage the io for the project """

import os
import warnings
import configparser # Careful not the .configparser
from datetime import datetime

import numpy as np


PATH_LEPHAREDIR = os.path.expanduser(os.getenv("LEPHAREDIR", default=None))
PATH_LEPHAREWORK = os.path.expanduser(os.getenv("LEPHAREWORK", default=None))



def get_now(format="%d%m%Y_%H%M%S"):
    """ See e.g. https://www.tutorialspoint.com/python/time_strftime.htm for the format """
    return datetime.now().strftime(format)

def get_default_path():
    """ """
    return PATH_LEPHAREWORK+"/pylephare/"+get_now()
    

_DEFAULT_FILTER_CONFIG = {'2mass': {'h': 'H.pb', 'j': 'J.pb', 'ks': 'Ks.pb'},
                          'megacam': {'g': 'gp.pb','r': 'rp.pb', 'u': 'up.pb','z': 'zp.pb','i': 'ip.pb'},
                          'sdss': {'u': 'up.pb','g': 'gp.pb',  'r': 'rp.pb',  'i': 'ip.pb',  'z': 'zp.pb'},
                          'galex': {'fuv': 'FUV.pb','fuvo': 'FUV_old.pb', 'nuv': 'NUV.pb', 'nuvo': 'NUV_old.pb'},
                          'ukidss': {'k': 'K.pb', 'h': 'H.pb', 'h2': 'H2.pb','j': 'J.pb',  'brg': 'Brg.pb','y': 'Y.pb','z': 'Z.pb'}
                          }
FILTER_CONFIGGILE = PATH_LEPHAREWORK+"/filt/config"


def get_filter_bandpass(filtername,**kwargs):
    """ """
    return IO().get_filter_bandpass(filtername,**kwargs)

def keys_to_filters(keys):
    """ """
    return [k for k in keys if k+".err" in keys]
    
class IO( object ):
    """ """
    LEPHAREDIR = PATH_LEPHAREDIR
    LEPHAREWORK = PATH_LEPHAREWORK
    def __init__(self, dirout=None):
        """ """
        self.set_dirout(dirout)
        self.load_filter_config()
    def set_dirout(self, dirout=None):
        """ """
        self._dirout = dirout

    def load_filter_config(self):
        """ """
        self._filtercongig = configparser.ConfigParser(allow_no_value=True)
        if os.path.isfile(FILTER_CONFIGGILE):
            self._filtercongig.read_file( open( FILTER_CONFIGGILE ) )
        else:
            warnings.warn("No %s file yet. This is building a default one"%FILTER_CONFIGGILE)
            self._filtercongig.read_dict(_DEFAULT_FILTER_CONFIG)
            with open(FILTER_CONFIGGILE,"w") as f:
                self._filtercongig.write(f)
        
    def get_config_filterlist(self, filterlist):
        """ get the list of lephare"""
        return ",".join([self.get_filterfile(filt_, fullpath=False) for filt_ in filterlist])

    def get_filterfile(self, filtername, fullpath=True):
        """ """
        source = (self.get_filt_path(from_work=False)+"/") if fullpath else ""
        return source+filtername.split(".")[0]+"/"+self._filtercongig.get(*filtername.split("."))

    def get_filter_bandpass(self, filtername, **kwargs):
        """ returns the sncosmo.Bandpass associated to the given filter """
        from  sncosmo import bandpasses
        wave, trans = np.asarray([l.split() for l in open(self.get_filterfile(filtername)).read().splitlines()
                                      if len(l)>0 and not l.startswith("#")], 
                                     dtype="float").T
        
        return bandpasses.Bandpass(wave, trans, **kwargs)
    # ================ #
    #    Statics       #
    # ================ #
    # // Dirout 
    @staticmethod
    def get_now(format="%d%m%Y_%H%M%S"):
        """ """
        return get_now(format)
    
    @staticmethod
    def get_default_path():
        """ """
        return get_default_path()

    
    # // Filters
    @staticmethod
    def keys_to_filters(keys):
        """ """
        return keys_to_filters(keys)
    
    @staticmethod
    def get_filt_path(from_work=True):
        """ """
        return PATH_LEPHAREWORK+"/filt" if from_work else PATH_LEPHAREDIR+"/filt"

    @staticmethod
    def is_filt_known(filtname):
        """ """
        return os.path.isfile(PATH_LEPHAREWORK+"/filt"+"/"+filtname)

    @staticmethod
    def get_instrument_filters(instrument):
        """ """
        return os.listdir(PATH_LEPHAREDIR+"/filt/"+instrument)
        
    @staticmethod
    def _get_instrument_filters():
        """ dictionary containing, """
        return {l:[l_ for l_ in os.listdir(PATH_LEPHAREDIR+"/filt/"+l) if l_.endswith(".pb")]
                    for l in os.listdir(PATH_LEPHAREDIR+"/filt") if os.path.isdir(PATH_LEPHAREDIR+"/filt/"+l)}
        
    # ================ #
    #    Methods       #
    # ================ #
    def does_dir_exist(self, directory, buildit=False):
        """ """
        if not os.path.isdir(directory):
            if not buildit:
                return False
            os.makedirs(directory)
            
        return True

    # ================ #
    #    Properties    #
    # ================ #
    @property
    def dirout(self):
        """ """
        if not hasattr(self,"_dirout") or self._dirout is None:
            self.set_dirout( self.get_default_path() )
            warnings.warn("Default dirout used: %s"%self._dirout)
        return self._dirout


