""" Module made to manage the io for the project """

import os
import warnings
import configparser # Careful not the .configparser
from datetime import datetime



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
        return ",".join([filt_.split(".")[0]+"/"+self._filtercongig.get(*filt_.split(".")) for filt_ in filterlist])
        
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
    def get_filt_path(from_work=True):
        """ """
        return PATH_LEPHAREWORK+"/filt" if from_work else PATH_LEPHAREDIR+"/filt"

    @staticmethod
    def is_filt_known(filtname):
        """ """
        return os.path.isfile(PATH_LEPHAREWORK+"/filt"+"/"+filtname)

    def get_instrument_filters(instrument):
        """ """
        return os.listdir(PATH_LEPHAREDIR+"/filt/"+l)
        
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


