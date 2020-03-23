""" Module made to manage the io for the project """

import os
import warnings
from datetime import datetime


PATH_LEPHAREDIR = os.path.expanduser(os.getenv("LEPHAREDIR", default=None))
PATH_LEPHAREWORK = os.path.expanduser(os.getenv("LEPHAREWORK", default=None))



def get_now(format="%d%m%Y_%H%M%S"):
    """ See e.g. https://www.tutorialspoint.com/python/time_strftime.htm for the format """
    return datetime.now().strftime(format)

def get_default_path():
    """ """
    return PATH_LEPHAREWORK+"/pylephare/"+get_now()
    

class IO( object ):
    """ """
    LEPHAREDIR = PATH_LEPHAREDIR
    LEPHAREWORK = PATH_LEPHAREWORK
    def __init__(self, dirout=None):
        """ """
        if dirout is None:
            self._dirout = self.get_default_path()
            warnings.warn("Default dirout used: %s"%self._dirout)

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
    def get_filt_path():
        """ """
        return PATH_LEPHAREWORK+"filt"

    @staticmethod
    def is_filt_known(filtname):
        """ """
        return os.path.isfile(PATH_LEPHAREWORK+"/filt/"+filtname)

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
        return self._dirout


