

import os
import warnings
import subprocess

import pandas
import numpy as np


from . import tools
from .io import PATH_LEPHAREDIR, PATH_LEPHAREWORK


class BC03Installer( object ):
    """ Test if the proper installation is made and enables to run it 

    USAGE:
    simply do: BC03Installer.install()

    """
    LEPHAREDIR = PATH_LEPHAREDIR  # os.getenv('SEDM_USER',default="auto")
    
    # ----------- #
    # Initialize  #
    # ----------- #
    def __init__(self):
        """ """
    
    @classmethod
    def install(cls, compiler="gfortran", reinstall=False, verbose=False):
        """ """
        this = cls()
        if not reinstall and this.has_ised_files():
            if verbose:
                print("No need to install ised files. set reinstall to True to force reinstallation")
            return None
        
        this.run_ised_installer(compiler="gfortran", verbose=verbose)
        return
    
    # =========== #
    #  Methods    #
    # =========== #
    def get_list_sedmodel(self, which="*"):
        """ """
        return [l for l in os.listdir(self.BC03_CHAB) if l.startswith("bc2003")
               and ((which in ["*", "all", "both"]) or 
                    ((which.lower() == "ascii") and l.endswith("ASCII")) or
                    ((which.lower() == "ised") and l.endswith("ised"))
                   )]

    def build_bc03_installer(self, compiler="gfortran"):
        """ """
        os.system("%s  -O5 -Wall -ffixed-line-length-132 %s -o %s"%(compiler,
                                                                    self.bc03_installer+".f",
                                                                    self.bc03_installer))
        
    def run_ised_installer(self, compiler="gfortran", verbose=False):
        """ """
        if not self.has_bc03_installer():
            if verbose:
                print("No bc03_installer. Creating one, using the '%s' compiler"%compiler)
            self.build_bc03_installer(compiler)
            if not self.has_bc03_installer():
                raise AttributeError("build_bc03_installer() ran but still no has_bc03_installer()... sorry")
        
        if verbose:
            print("Creating the .ised files from the .ised_ASCII ones")
                  
        for file in self.get_list_sedmodel('ASCII'):
            report = os.system("%s %s"%(self.bc03_installer, self.BC03_CHAB+"/"+file))
     
    # ----------- #
    # Has tests   #
    # ----------- #
    def has_lephare_env(self):
        """ test that the global LEPHAREDIR is defined"""
        return self.LEPHAREDIR is not None

    def has_bc03_data(self):
        """ """
        if not self.has_lephare_env():
            raise IOError("$LEPHAREDIR is not defined ")
        
        return os.path.isdir( self.BC03_CHAB )
    
    def has_bc03_installer(self):
        """ """
        if not self.has_bc03_data():
            raise IOError("BC03_CHAB lephare gal library not downloaded."+
                          "\n See http://www.cfht.hawaii.edu/~arnouts/LEPHARE/install.html")
            
        return os.path.isfile( self.bc03_installer)
    
    def has_ised_files(self):
        """ """
        return len(self.get_list_sedmodel("ised"))>0
    
    # =========== #
    #  Properties #
    # =========== #
    @property
    def BC03_CHAB(self):
        """ """
        return self.LEPHAREDIR+ "/sed/GAL/BC03_CHAB"
    
    @property
    def bc03_installer(self):
        """ """
        return self.BC03_CHAB+ "/bin_ised"


######################
#                    #
#  LEPHARE           #
#                    #
######################
# ================== #
#                    #
#  Virtual LEPHARE   #
#                    #
# ================== #
class _LePhareBase_( object ):
    """ 
    This Basic Class contains the data - filter - config interactions
    """
    def __init__(self, data=None,  configfile=None, dirout=None, dataunits="AA"):
        """ """
        from .io import IO
        self.io = IO(dirout)

        if data is not None:
            self.set_data(data=data, unit=dataunits)
            
        if configfile is not None:
            self.set_config(configfile)

        
    @classmethod
    def read_csv(cls, catin, configfile=None, filters=None, dirout=None, sep=",", **kwargs):
        """ """
        return cls( pandas.read_csv(catin, sep=sep, **kwargs),
                    configfile=configfile, filters=filters, dirout=dirout,
                  )
    # =============== #
    #  Methods        #
    # =============== #
    # -------- #
    #  SETTER  #
    # -------- #
    def set_data(self, data=None, inhz=False, **kwargs):
        """
        Set up the file paths about the data, the config files (input and output) and the results path.
        
        Parameters
        ----------
        data : [string or pandas.DataFrame or dict]
            Path of the data file or a DataFrame/dict, both of them in a format readable by LePhare fitter.
            
        filters : [list(string) or None]
            List of filters of the given measurements, the 'FILTER_LIST' parameter in the configuration 
            file will be changed to the corresponding list of LePhare file names.
            The filter syntax must be like "project.band" (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).
            If None, the filter list will be based on the configuration file.
        
        inhz:
            set to true if the flux (and flux) are given in erg/s/cm2/Hz
            (False means in erg/s/cm2/AA)

        
        Returns
        -------
        Void
        """
        
        # // data format
        if type(data) == dict:
            data = pandas.DataFrame(data)
            
        elif type(data) != pandas.DataFrame:
            raise TypeError("Data must be dict or pandas.DataFrame")

        data_ = data.copy()

        # // filters
        filters = [k for k in data.columns if k+".err" in data.columns]
        self.set_filters(filters) 
        
        # // flux units
        if not inhz:
            for filter_ in self.filters:
                data_[filter_]        = tools.flux_aa_to_hz(data_[filter_], self.io.get_filter_bandpass(filter_).wave_eff)
                data_[filter_+".err"] = tools.flux_aa_to_hz(data_[filter_+".err"], self.io.get_filter_bandpass(filter_).wave_eff)
        # // Setting
        self._data = data_
        
    def set_config(self, configfile):
        """ """
        from .configparser import ConfigParser
        self._config = ConfigParser.read(configfile)
        self.config.set_fileout(self.io.dirout+"/config")
        if self.has_filters():
            self.config.set_value("FILTER_LIST", self.io.get_config_filterlist(self._filters))
            self.config.set_value("FILTER_FILE", "_".join([l.replace(".","") for l in self.filters])+".filt")
        
    def set_filters(self, filters):
        """ """
        self._filters = filters
        if self.has_config():
            self.config.set_value("FILTER_LIST", self.io.get_config_filterlist(self._filters))
            self.config.set_value("FILTER_FILE", "_".join([l.replace(".","") for l in self.filters])+".filt")
        
        
    def set_redshift(self, redshift, index=None):
        """ """
        self.data.at[index, "zspec"] = redshift

    def set_context(self, context_alue, index=None):
        """ 
        The 'context' number must be : sum(2**[band_nb]). 
            For example, bands = ["u", "g", "r", "i", "z"] (bands_nb = [0, 1, 2, 3, 4]) :

        - context = 31 --> ["u", "g", "r", "i", "z"]
        - context = 30 --> ["g", "r", "i", "z"]
        - context = 15 --> ["u", "g", "r", "i"]
        - context = 25 --> ["u", "i", "z"]

        """
        self.data.at[index, "context"] = context_value

    def set_dirout(self, dirout):
        """ """
        self.io.set_dirout(dirout)
        
    # -------- #
    #  GETTER  #
    # -------- #
    def get_filters_context(self, filters=None):
        """ get the 'CONTEXT' value in 'data_meas' given the list filters to fit on.
        
        The 'context' number must be : sum(2**[band_nb]). 
            For example, bands = ["u", "g", "r", "i", "z"] (bands_nb = [0, 1, 2, 3, 4]) :

        - context = 31 --> ["u", "g", "r", "i", "z"]
        - context = 30 --> ["g", "r", "i", "z"]
        - context = 15 --> ["u", "g", "r", "i"]
        - context = 25 --> ["u", "i", "z"]

        Parameters
        ----------
        id : [int or list(int) or None]
            Index(es) of the line(s) you want to change the context by the given filters.
            If None, the context is changed for every line.
        
        filters : [list(string)]
            List of filters to fit on.
        
        
        Returns
        -------
        Void
        """
        if filters is None:
            filters = self.filters
            
        idx_cross = np.argwhere(np.in1d(self.filters, filters)).ravel()
        context = np.sum([2**ii for ii in idx_cross])
        return self.get_filters_context()

    def get_configfile(self, original=False, update=False):
        """ """
        return self.config._filename if original else self.config.get_fileout(buildit=True, update=update)

    def get_datafile(self, catfile=None):
        """ """
        if catfile is None:
            catfile = self.io.dirout+"/data.csv"

        # Strange LePhare format
        self.data.to_csv(catfile, sep=" ", header=False)
        self.config.set_value("CAT_IN", catfile)
        return catfile
    
    # ============== #
    #  Properties    #
    # ============== #
    # - Data
    @property
    def data(self):
        """ DataFrame containing the data """
        if not hasattr(self,"_data"):
            self._data = None

        return self._data

    def has_data(self):
        """ """
        return self.data is not None and len(self.data)>0
    
    @property
    def ndata(self):
        """ """
        if not self.has_data():
            raise AttributeError("No data set yet.")
        return len(self.data)

    # - Filters
    @property
    def filters(self):
        """ """
        if not hasattr(self,"_filters"):
            self._filters = None
        return self._filters

    def has_filters(self):
        """ """
        return self.filters is not None and len(self.filters)>0
    
    @property
    def nfilters(self):
        """ """
        if not self.has_filters():
            raise AttributeError("No filter set yet.")
        return len(self.filters)
    
    # - Config
    @property
    def config(self):
        """ ConfigParser object to interact with the configuration file """
        if not hasattr(self,"_config"):
            self._config = None
        return self._config

    def has_config(self):
        """ """
        return self.config is not None


# ================== #
#                    #
#  MAIN LEPHARE      #
#                    #
# ================== #
class LePhare( _LePhareBase_ ):
    """ _LePhareBase + command line tools """
    
    # -------- #
    #  RUNNER  #
    # -------- #
    #
    # - Main
    #
    def run(self, update=False, filters=None, contextid=None, dirout=None,
                configfile=None, catinfile=None, originalconfig=False):
        """
        Then execute "$LEPHAREDIR/source/zphota -c [...].para" in the shell.
        
        Options
        -------
        filters : [list or None]
            If a list of filters is given, the context will be changed for every line to the corresponding one.
            The filter syntax must be like 'project.band' (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).
        
        dirout : [string or None]
            Where the data should be stored.
            Default is 'None', which is the path set during the class construction or an execution of 'set_data'.
                        
        Returns
        -------
        Void
        """
        fileout = {}
        if filters is not None:
            self.set_context(self.get_filters_context(filters), index=contextid)

        if dirout is None:
            dirout = self.io.dirout

        if configfile is None:
            if catinfile is None:
                catinfile = self.get_datafile()
            else:
                self.config.set_value("CAT_IN", catinfile)
            configfile = self.get_configfile(originalconfig, update=True)
        else:
            catinfile = "see given configfile"

        # - Initialize
        self.run_init(update=update, configfile=configfile)

        # - Run the fit        
        cmd = "{}/source/zphota -c {} -CAT_OUT {}".format(self.io.LEPHAREDIR, configfile, dirout+"/catout")
        try:
            subprocess.run(cmd.split())
        except:
            raise ValueError("LePhareError : unable to run 'zphota'.")

        fileout["config"] = configfile
        fileout["catin"]  = catinfile
        fileout["catout"] = dirout+"/catout"
        
        # relocate the output spectra if any.
        if self.config.get_value("SPEC_OUT") in ["True", True, "yes", "YES","Yes"]:
            specfiles = [l for l in os.listdir(".") if l.startswith('Id') and l.endswith(".spec")]
            new_location = [dirout+"/"+l for l in specfiles]
            print(specfiles)
            print("moved to")
            _ = [os.rename(l, nl) for l, nl in zip(specfiles,new_location)]
            fileout["spec"] =  new_location            
        else:
            fileout["spec"] = None

        return fileout
        
    #
    # - Secondary
    #
    def build_filter_files(self, update=False, configfile=None, updateconfig=True):
        """ """
        filter_file = self.config.get_value("FILTER_FILE")
        if not self.io.is_filt_known(filter_file) or update:
            cmd = "{}/source/filter -c {}".format(self.io.LEPHAREDIR, self.get_configfile(update=updateconfig) if configfile is None else configfile)
            try:
                subprocess.run(cmd.split())
            except:
                raise ValueError("LePhareError : unable to run 'filter'.")

    def run_sedtolib(self, update=False, configfile=None, updateconfig=True):
        """
        Execute "$LEPHAREDIR/source/sedtolib -t [S,Q,G] -c [...].para" in the shell.
        Exception : if the "$LEPAHAREWORK/lib_bin/[...].bin" files already exist, the command is not executed, unless 'update' is True.
        
        Options
        -------
        
        update : [bool]
            Set to True if you want to force the execution of the command.
        
        Returns
        -------
        Void
        """
    
        lib_s_v = self.config.get_value("STAR_LIB")
        lib_q_v = self.config.get_value("QSO_LIB")
        lib_g_v = self.config.get_value("GAL_LIB")
        
        if not np.prod([os.path.isfile(self.io.LEPHAREWORK+"/lib_bin/"+elt+".bin")
                        for elt in [lib_s_v, lib_q_v, lib_g_v]]) or update:
            for elt in ["S", "Q", "G"]:
                cmd = "{}/source/sedtolib -t {} -c {}".format(self.io.LEPHAREDIR, elt, self.get_configfile(update=updateconfig) if configfile is None else configfile)
                try:
                    subprocess.run(cmd.split())
                except:
                    raise ValueError("LePhareError : unable to run 'sedtolib' for '{}'.".format(elt))

    def run_mag_star(self, update=False, configfile=None, updateconfig=True):
        """
        Execute "$LEPHAREDIR/source/mag_star -c [...].para" in the shell.
        Exception : if the "$LEPAHAREWORK/lib_mag/[...].bin" file already exists, 
                    the command is not executed, unless 'update' is True.
        
        Options
        -------
        update : [bool]
            Set to True if you want to force the execution of the command.
        
        Returns
        -------
        Void
        """
        
        lib_s_v = self.config.get_value("STAR_LIB_OUT")
        
        if not os.path.isfile(self.io.LEPHAREWORK+"/lib_mag/"+lib_s_v+".bin") or update:
            cmd = "{}/source/mag_star -c {}".format(self.io.LEPHAREDIR, self.get_configfile(update=updateconfig) if configfile is None else configfile)
            try:
                subprocess.run(cmd.split())
            except:
                raise ValueError("LePhareError : unable to run 'mag_star'.")

    def run_mag_gal(self, update=False, configfile=None, updateconfig=True):
        """
        Execute "$LEPHAREDIR/source/mag_gal -t [Q,G] -c [...].para" in the shell.
        Exception : if the "$LEPAHAREWORK/lib_mag/[...].bin" files already exist, 
                    the command is not executed, unless 'update' is True.
        
        Options
        -------
        update : [bool]
            Set to True if you want to force the execution of the command.
                        
        Returns
        -------
        Void
        """
        lib_q_v = self.config.get_value("QSO_LIB_OUT")
        lib_g_v = self.config.get_value("GAL_LIB_OUT")
        
        if not np.prod([os.path.isfile(self.io.LEPHAREWORK+"/lib_mag/"+elt+".bin") for elt in [lib_q_v, lib_g_v]]) or update:
            for elt in ["Q", "G"]:
                cmd = "{}/source/mag_gal -t {} -c {}".format(self.io.LEPHAREDIR, elt, self.get_configfile(update=updateconfig) if configfile is None else configfile)
                try:
                    subprocess.run(cmd.split())
                except:
                    raise ValueError("LePhareError : unable to run 'mag_gal' for '{}'.".format(elt))

    def run_init(self, update=False, configfile=None, updateconfig=True):
        """
        Run shell commands to initialize LePhare fitting.
        
        Options
        -------
        update : [bool]
            Set to True if you want to force the execution of the command.
        
        Returns
        -------
        Void
        """
        prop = dict(update=update, configfile=configfile, updateconfig=updateconfig)
        self.build_filter_files(**prop)
        self.run_sedtolib(**prop)
        self.run_mag_star(**prop)
        self.run_mag_gal(**prop)
        
# ================== #
#                    #
#  LEPHARE RESULTS   #
#                    #
# ================== #
class LePhareResults( _LePhareBase_ ):
    """ """
    def __init__(self, data=None, results=None,
                     spectrafiles=None, configfile=None, **kwargs):
        """ """
        
        return super().__init__(data=data, configfile=configfile, **kwargs)

    @classmethod
    def read_csv(cls, datafile, resultsfile=None, spectrafiles=None,
                     configfile=None, sep=",", **kwargs):
        """ """
        return cls( data = pandas.read_csv(datafile, sep=sep, **kwargs),
                    results = pandas.read_csv(resultsfile, sep=sep, **kwargs),
                    spectrafiles=spectrafiles,
                    configfile=configfile, **kwargs
                  )


    # ============== #
    #  Methods       #
    # ============== #
    # ------- #
    # SETTER  #
    # ------- #    
    def set_results(self, results):
        """ """
        self._results = results

    def set_spectra(self, spectralfiles):
        """ """
        self._spectra = spectralfiles

    # ------- #
    # GETTER  #
    # ------- #    

    # ------- #
    # PLOTTER #
    # ------- #    
    def show(self, ax=None):
        """ """
        import matplotlib.pyplot as mpl

    # ============== #
    #  Properties    #
    # ============== #
    @property
    def results(self):
        """ """
        if not hasattr(self,"_results"):
            self._results = None
        return self._results

    def has_results(self):
        """ """
        return self.results is not None

    @property
    def spectra(self):
        """ """
        if not hasattr(self,"_spectra"):
            self._spectra = None
        return self._spectra

    def has_spectra(self):
        """ """
        return self.spectra is not None
