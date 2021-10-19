""" Module made to manage the io for the project """

import os
import warnings
import configparser # Careful not the .configparser
from datetime import datetime

import numpy as np


PATH_LEPHAREDIR = os.path.expanduser(os.getenv("LEPHAREDIR", default=None))
PATH_LEPHAREWORK = os.path.expanduser(os.getenv("LEPHAREWORK", default=None))
PATH_PYLEPHAREWORK = PATH_LEPHAREWORK+"/pylephare"


def get_now(format="%d%m%Y_%H%M%S"):
    """ Return a string filled with the current day and time.
    See e.g. https://www.tutorialspoint.com/python/time_strftime.htm for the format. """
    return datetime.now().strftime(format)

def get_default_path():
    """ Return a folder directory, located in $LEPHAREWORK/pylephare and named by the current day and time. """
    return PATH_PYLEPHAREWORK+"/"+get_now()


_PACKAGE_ROOT = os.path.abspath(os.path.dirname(__file__))+"/"
_DEFAULT_FILTER_CONFIG = {"2mass": {"h": "H.pb", "j": "J.pb", "ks": "Ks.pb"},
                          "megacam": {"g": "gp.pb","r": "rp.pb", "u": "up.pb","z": "zp.pb","i": "ip.pb"},
                          "sdss": {"u": "up.pb","g": "gp.pb",  "r": "rp.pb",  "i": "ip.pb",  "z": "zp.pb"},
                          "galex": {"fuv": "FUV.pb","fuvo": "FUV_old.pb", "nuv": "NUV.pb", "nuvo": "NUV_old.pb"},
                          "ukidss": {"k": "K.pb", "h": "H.pb", "h2": "H2.pb","j": "J.pb",  "brg": "Brg.pb","y": "Y.pb","z": "Z.pb"}, 
                          "ps1": {"g": "g_ps.pb", "r": "r_ps.pb", "i": "i_ps.pb", "z": "z_ps.pb", "y": "y_ps.pb"}
                          }
FILTER_CONFIGFILE = PATH_LEPHAREWORK+"/filt/config"

DEFAULTCONFIG = PATH_PYLEPHAREWORK+"/config/default.config"
DEFAULTCONFIG_OUT = PATH_PYLEPHAREWORK+"/config/default_output.config"


# Create the pylephare environment
if not os.path.isdir(PATH_PYLEPHAREWORK):
    os.mkdir(PATH_PYLEPHAREWORK)
    
if not os.path.isdir(PATH_PYLEPHAREWORK+"/config"):
    os.mkdir(PATH_PYLEPHAREWORK+"/config")
    
if not os.path.isfile(DEFAULTCONFIG):
    from shutil import copyfile
    copyfile(_PACKAGE_ROOT+"/config/default.config", DEFAULTCONFIG)
    
if not os.path.isfile(DEFAULTCONFIG_OUT):
    from shutil import copyfile
    copyfile(_PACKAGE_ROOT+"/config/default_output.config", DEFAULTCONFIG_OUT)

    
def get_filter_bandpass(filtername,**kwargs):
    """
    Return the sncosmo.Bandpass associated to the given filter
    
    Parameters
    ----------
    filtername : [string]
        Filter name to get the bandpass.
        Must be on the format instrument.band (e.g. sdss.u, ps1.g, ...).
    
    **kwargs
        sncosmo.bandpasses.Bandpass(**kwargs)
    
    
    Returns
    -------
    float
    """
    return IO().get_filter_bandpass(filtername,**kwargs)

def filterfile_to_filtername(filterfile):
    """
    Convert LePhare filterfile into filtername(e.g. sdss/up.pb into sdss.u)
    
    Parameters
    ----------
    filterfile : [string]
        Filter file as it appears in the "FILTER_LIST" LePhare configuration file parameter.
    
    
    Returns
    -------
    string
    """
    _io = IO()
    if len(np.atleast_1d(filterfile))==1:
        return _io.get_filtername(filterfile)
    return [_io.get_filtername(filt) for filt in np.atleast_1d(filterfile)]

def filtername_to_filterfile(filtename, fullpath=True):
    """
    Convert filterfile name into a LePhare filterfile (e.g. sdss.u into sdss/up.pb)
    
    Parameters
    ----------
    filtername : [string]
        Filter name to convert into LePhare filterfile.
        Must be on the format instrument.band (e.g. sdss.u, ps1.g, ...).
    
    Options
    -------
    fullpath : [bool]
        If True, return the full path of the corresponding filterfile.
        Default is False.
    
    
    Returns
    -------
    string
    """
    _io = IO()
    if len(np.atleast_1d(filtename))==1:
        return _io.get_filterfile(filtename, fullpath=fullpath)
    return [_io.get_filterfile(filt,fullpath=fullpath) for filt in np.atleast_1d(filtename)]
    
def keys_to_filters(keys):
    """
    Return a list filters included in the given list of keys.
    To be detected, "{filter_name}.err" must appears in the given list of keys.
    
    Parameters
    ----------
    keys : [list(string)]
        List of keys from which to extract the filters.
    
    
    Returns
    -------
    list(string)
    """
    return [k for k in keys if k+".err" in keys]
    
class IO( object ):
    """ Class handling IO for the object """
    
    LEPHAREDIR = PATH_LEPHAREDIR
    LEPHAREWORK = PATH_LEPHAREWORK
    
    def __init__(self, dirout=None):
        """
        Class builder.
        Set the output directory as attribute and load the filter configuration file (create one if needed ; $LEPHAREWORK/filt/config).
        
        Parameters
        ----------
        dirout : [string or None]
            Output directory.
            If None, create a default one (based on current day and time).
            Default is None.
        
        
        Returns
        -------
        Void
        """
        self.set_dirout(dirout)
        self.load_filter_config()
    
    def set_dirout(self, dirout=None):
        """
        Set the output directory as attribute.
        
        Parameters
        ----------
        dirout : [string or None]
            Output directory.
            If None, create a default one (based on current day and time).
            Default is None.
        
        
        Returns
        -------
        Void
        """
        self._dirout = dirout

    def load_filter_config(self):
        """
        Load the filter configuration file (create one if needed).
        ($LEPHAREWORK/filt/config)
        
        Returns
        -------
        Void
        """
        self._filtercongig = configparser.ConfigParser(allow_no_value=True)
        if os.path.isfile(FILTER_CONFIGFILE):
            self._filtercongig.read_file( open( FILTER_CONFIGFILE ) )
        else:
            warnings.warn("No %s file yet. This is building a default one"%FILTER_CONFIGFILE)
            self._filtercongig.read_dict(_DEFAULT_FILTER_CONFIG)
            with open(FILTER_CONFIGFILE,"w") as f:
                self._filtercongig.write(f)
        
    def get_config_filterlist(self, filterlist):
        """
        Return the value to write in the "FILTER_LIST" parameter into the configuration file.
        
        Parameters
        ----------
        filterlist : [list(string)]
            List of filters to convert into the "FILTER_LIST" parameter's value format.
            Filters must be on the format instrument.band (e.g. sdss.u, ps1.g, ...).
        
        
        Returns
        -------
        string
        """
        return ",".join([self.get_filterfile(filt_, fullpath=False) for filt_ in filterlist])

    def get_filterfile(self, filtername, fullpath=True):
        """
        Convert filterfile name into a LePhare filterfile (e.g. sdss.u into sdss/up.pb)
        
        Parameters
        ----------
        filtername : [string]
            Filter name to convert into LePhare filterfile.
            Must be on the format instrument.band (e.g. sdss.u, ps1.g, ...).
        
        Options
        -------
        fullpath : [bool]
            If True, return the full path of the corresponding filterfile.
            Default is False.
        
        
        Returns
        -------
        string
        """
        source = (self.get_filt_path(from_work=False)+"/") if fullpath else ""
        return source+filtername.split(".")[0]+"/"+self._filtercongig.get(*filtername.split("."))

    def get_filtername(self, filterfile):
        """
        Convert LePhare filterfile into filtername(e.g. sdss/up.pb into sdss.u)
        
        Parameters
        ----------
        filterfile : [string]
            Filter file as it appears in the "FILTER_LIST" LePhare configuration file parameter.
        
        
        Returns
        -------
        string
        """
        inst,bandin = filterfile.split("/")[-2:]
        band = [k for k,v in self._filtercongig[inst].items() if v==bandin]
        if len(band)==0:
            raise ValueError("Cannot parse filt, unknown band: %s"%bandin)
        return inst+".%s"%band[0]
    
    def get_filter_bandpass(self, filtername, **kwargs):
        """
        Return the sncosmo.Bandpass associated to the given filter
        
        Parameters
        ----------
        filtername : [string]
            Filter name to get the bandpass.
            Must be on the format instrument.band (e.g. sdss.u, ps1.g, ...).
        
        **kwargs
            sncosmo.bandpasses.Bandpass(**kwargs)
        
        
        Returns
        -------
        float
        """
        from sncosmo import bandpasses
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
        """
        Return a string filled with the current day and time.
        See e.g. https://www.tutorialspoint.com/python/time_strftime.htm for the format.
        
        Parameters
        ----------
        format : [string]
            Current day and time format (https://www.tutorialspoint.com/python/time_strftime.htm).
            Default is "%d%m%Y_%H%M%S".
        
        
        Returns
        -------
        string
        """
        return get_now(format)
    
    @staticmethod
    def get_default_path():
        """
        Return a folder directory, located in $LEPHAREWORK/pylephare and named by the current day and time.
        
        
        Returns
        -------
        Void
        """
        return get_default_path()

    
    # // Filters
    @staticmethod
    def keys_to_filters(keys):
        """
        Return a list filters included in the given list of keys.
        To be detected, "{filter_name}.err" must appears in the given list of keys.
        
        Parameters
        ----------
        keys : [list(string)]
            List of keys from which to extract the filters.
        
        
        Returns
        -------
        list(string)
        """
        return keys_to_filters(keys)
    
    @staticmethod
    def get_filt_path(from_work=True):
        """
        Return the the directory for the filt/ folder.
        
        Options
        -------
        from_work : [bool]
            If True, return $LEPHAREWORK/filt.
            If False, return $LEPHAREDIR/filt.
            Default is True.
        
        
        Returns
        -------
        string
        """
        return PATH_LEPHAREWORK+"/filt" if from_work else PATH_LEPHAREDIR+"/filt"

    @staticmethod
    def is_filt_known(filtname):
        """
        Test that the given filt file exists.
        
        Parameters
        ----------
        filtname : [string]
            Filt file name.
        
        
        Returns
        -------
        bool
        """
        return os.path.isfile(PATH_LEPHAREWORK+"/filt"+"/"+filtname)

    @staticmethod
    def get_instrument_filters(instrument):
        """
        Return the list of LePhare filterfile directories available in $LEPHAREDIR/filt for the given instrument.
        
        Parameters
        ----------
        instrument : [string]
            Instrument from which to get the list of available filterfile directories.
        
        
        Returns
        -------
        list(string)
        """
        return os.listdir(PATH_LEPHAREDIR+"/filt/"+instrument.lower())
        
    @staticmethod
    def _get_instrument_filters():
        """
        Return a dictionary containing for every available instrument in $LEPHAREDIR/filt a list of available LePhare filterfile directories.
        
        
        Returns
        -------
        dict
        """
        return {l:[l_ for l_ in os.listdir(PATH_LEPHAREDIR+"/filt/"+l) if l_.endswith(".pb")]
                for l in os.listdir(PATH_LEPHAREDIR+"/filt") if os.path.isdir(PATH_LEPHAREDIR+"/filt/"+l)}
        
    # ================ #
    #    Methods       #
    # ================ #
    def does_dir_exist(self, directory, buildit=False):
        """
        Test if a given directory exist or not (if not, can build it).
        
        Parameters
        ----------
        directory : [string]
            Directory to check on the existence.
        
        Options
        -------
        buildit : [bool]
            If True, make the given directory if it doesn't exist.
            Default is False.
        
        
        Returns
        -------
        bool
        """
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
        """ Output directory """
        if not hasattr(self,"_dirout") or self._dirout is None:
            self.set_dirout( self.get_default_path() )
            warnings.warn("Default dirout used: %s"%self._dirout)
        return self._dirout


