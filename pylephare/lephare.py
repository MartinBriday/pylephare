

import os
import shutil
import warnings
import subprocess

import pandas
import numpy as np


from . import tools
from .io import PATH_LEPHAREDIR, PATH_LEPHAREWORK, _DEFAULT_FILTER_CONFIG
from .base import _FilterHolder_

class BC03Installer( object ):
    """ Test if the proper installation is made and enables to run it 

    USAGE:
    simply do: BC03Installer.install()

    """
    LEPHAREDIR = PATH_LEPHAREDIR  # os.getenv('SEDM_USER',default="auto")
    
    # =========== #
    #  Installer  #
    # =========== #
    @classmethod
    def install(cls, compiler="gfortran", reinstall=False, verbose=False):
        """
        Install, if needed, the BC03 libraries in LePhare package.
        
        Parameters
        ----------
        compiler : [string]
            Compiler to use to build the BC03 libraries installer.
            Default is "gfortran".
        
        Options
        -------
        reinstall : [bool]
            Run the installation, even if the libraries already are installed.
            Default is False.
        
        verbose : [bool]
            Print informations.
            Default is False.
        
        
        Returns
        -------
        None
        """
        this = cls()
        if not reinstall and this.has_ised_files():
            if verbose:
                print("No need to install ised files. Set 'reinstall' to True to force reinstallation")
            return None
        
        this.run_ised_installer(compiler=compiler, verbose=verbose)
        return
    
    # =========== #
    #  Methods    #
    # =========== #
    def get_list_sedmodel(self, which="*"):
        """
        Return the list of files in the chosen format.
        
        Parameters
        ----------
        which : [string]
            File format choice between "ASCII", "ised" or both ("*") to rturn the list of files.
        
        
        Returns
        -------
        list(string)
        """
        return [l for l in os.listdir(self.BC03_CHAB) if l.startswith("bc2003")
                and ((which in ["*", "all", "both"]) or
                     ((which.lower() == "ascii") and l.endswith("ASCII")) or
                     ((which.lower() == "ised") and l.endswith("ised"))
                    )]

    def build_bc03_installer(self, compiler="gfortran"):
        """
        Execute the bash command to build the BC03 libraries installer.
        
        Parameters
        ----------
        compiler : [string]
            Compiler to use to build the BC03 libraries installer.
            Default is "gfortran".
        
        
        Returns
        -------
        Void
        """
        os.system("%s -O5 -Wall -ffixed-line-length-132 %s -o %s"%(compiler,
                                                                   self.bc03_installer+".f",
                                                                   self.bc03_installer))
        
    def run_ised_installer(self, compiler="gfortran", verbose=False):
        """
        Run the BC03 libraries installer on each file.
        
        Parameters
        ----------
        compiler : [string]
            Compiler to use to build the BC03 libraries installer.
            Default is "gfortran".
        
        Options
        -------
        verbose : [bool]
            Print informations.
            Default is False.
        
        
        Returns
        -------
        Void
        """
        if not self.has_bc03_installer():
            if verbose:
                print("No bc03_installer. Creating one, using the '%s' compiler"%compiler)
            self.build_bc03_installer(compiler)
            if not self.has_bc03_installer():
                raise AttributeError("build_bc03_installer() ran but still no has_bc03_installer()... sorry")
        
        if verbose:
            print("Creating the .ised files from the .ised_ASCII ones")
                  
        for file in self.get_list_sedmodel("ASCII"):
            report = os.system("%s %s"%(self.bc03_installer, self.BC03_CHAB+"/"+file))
     
    # ----------- #
    # Has tests   #
    # ----------- #
    def has_lephare_env(self):
        """ Test that the global LEPHAREDIR is defined """
        return self.LEPHAREDIR is not None

    def has_bc03_data(self):
        """ Test that the BC03 directory exists """
        if not self.has_lephare_env():
            raise IOError("$LEPHAREDIR is not defined ")
        
        return os.path.isdir( self.BC03_CHAB )
    
    def has_bc03_installer(self):
        """ Test that the fortran code is compiled and thus taht the BC03 installer exists """
        if not self.has_bc03_data():
            raise IOError("BC03_CHAB lephare gal library not downloaded."+
                          "\nSee http://www.cfht.hawaii.edu/~arnouts/LEPHARE/install.html")
            
        return os.path.isfile( self.bc03_installer)
    
    def has_ised_files(self):
        """ Test that .ised files exist """
        return len(self.get_list_sedmodel("ised"))>0
    
    # =========== #
    #  Properties #
    # =========== #
    @property
    def BC03_CHAB(self):
        """ BC03 directory in the LePhare package """
        return self.LEPHAREDIR+ "/sed/GAL/BC03_CHAB"
    
    @property
    def bc03_installer(self):
        """ BC03 installer directory """
        return self.BC03_CHAB+ "/bin_ised"


######################
#                    #
#  LEPHARE           #
#                    #
######################
def read_catout(filein, filternames=None):
    """
    Read a LePhare output catalog and return a pandas.DataFrame.
    
    Parameters
    ----------
    filein : [string]
        LePhare output catalog directory.
    
    Options
    -------
    filternames : [list(string)]
        List of filter names to use in column names.
        If None, automatically look into the output catalog to find them.
        Default is None.
    
    
    Returns
    -------
    pandas.DataFrame
    """
    def _get_filternames_(lines):
        i_filter_file = [ii for ii, _line in enumerate(lines) if "FILTER_FILE" in _line][0]
        _filter_file = lines[i_filter_file].split(" ")[-1]
        _lp_filters = [ii for ii in lines[i_filter_file+1].split(" ") if ii not in ["", "#", "AB", "VEGA"]]
        _filternames = [ii for ii in _filter_file.split(".filt")[0].split("_")]
        filternames = [jj+"."+ii.split(jj)[-1] for ii in _filternames for jj in _DEFAULT_FILTER_CONFIG.keys() if jj in ii]
        return filternames if len(filternames)==len(_lp_filters) else _lp_filters
    
    with open(filein) as f_:
        d = f_.read().splitlines()
        i_col_names = [ii for ii, line in enumerate(d) if "Output format" in line][0]
        i_header = [i for i, line in enumerate(d) if line.startswith("#")]
        if filternames is None:
            filternames = _get_filternames_(d)
        colname_, colid = np.concatenate([[l_.replace("#","").split() for l_ in l.split(",") if len(l_)>1] for l in d[i_col_names+1:i_header[-1]]],
                                         axis=0).T
        colname = np.concatenate([[l.replace("()","_%s"%f_) for f_ in filternames] if "()" in l else [l] for l in colname_])
        datain = np.asarray([l.split() for l in np.asarray(d)[i_header[-1]+1:]])
        return pandas.DataFrame(datain, columns=colname).set_index("IDENT")
    
# ================== #
#                    #
#  Virtual LEPHARE   #
#                    #
# ================== #
class _LePhareBase_( _FilterHolder_ ):
    """ 
    This Basic Class contains the data - filter - config interactions
    """
    def __init__(self, data=None,  configfile=None, dirout=None, inhz=False):
        """
        Class builder.
        Deal with the input data, the configuration file and the output directory.
        
        Parameters
        ----------
        data : [pandas.DataFrame or dict]
            Input data in a compatible LePhare format. Here are the columns style to adopt:
                filtername0, filtername0.err, filtername1, filtername1.err, ... filternameN, filternameN.err, CONTEXT, Z-SPEC, STRING
            where 'filtername{}' is a known filter by LePhare with the format instrument.band, for instance sdss.u or ps1.z, given in flux.
            The 'CONTEXT' sets the filters to use among the provided ones in input data for the SED fitting in each row.
            It is the sum of 2 to the {i}, {i} being the filtername numbers to use in each row.
            For example, using the filters = ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] in data input:
                * ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3 + 2^4 = 31
                * ["sdss.g", "sdss.r", "sdss.i", "sdss.z"]           --> CONTEXT = 2^1 + 2^2 + 2^3 + 2^4       = 30
                * ["sdss.u", "sdss.g", "sdss.r", "sdss.i"]           --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3       = 15
                * ["sdss.u", "sdss.i", "sdss.z"]                     --> CONTEXT = 2^0 + 2^3 + 2^4             = 25
            If 'CONTEXT' column doesn't exist, automatically create one with the maximum CONTEXT value.
            'Z-SPEC' is the spectroscopic redshift, if available.
            'STRING' and all following columns are not used for the SED fitting. Allow the user to store additional informations.
        
        configfile : [string or None]
            Configuration file directory.
            If None, use the default one ($LEPHAREWORK/pylephare/config/default.config).
            Default is None.
        
        dirout : [string or None]
            Setting the output folder directory.
            If None, create a default directory based on date and time in $LEPHAREWORK/pylephare/.
            Default is None.
        
        Options
        -------
        inhz : [bool]
            Set to True if the fluxes (and flux uncertainties) are given in erg/s/cm2/Hz ; False means in erg/s/cm2/AA.
            Default is False.
        
        
        Returns
        -------
        Void
        """
        from .io import IO
        self.io = IO(dirout)

        if data is not None:
            self.set_data(data=data, inhz=inhz)

        self.set_config(configfile) 

        
    @classmethod
    def read_csv(cls, catin, configfile=None, dirout=None, inhz=False, sep=",", **kwargs):
        """
        Build an instance giving a data directory.
        
        Parameters
        ----------
        catin : [string]
            Input data directory.
        
        configfile : [string or None]
            Configuration file directory.
            If None, use the default one ($LEPHAREWORK/pylephare/config/default.config).
            Default is None.
        
        dirout : [string or None]
            Setting the output folder directory.
            If None, create a default directory based on date and time in $LEPHAREWORK/pylephare/.
            Default is None.
        
        Options
        -------
        inhz : [bool]
            Set to True if the flux (and flux uncertainties) are given in erg/s/cm2/Hz ; False means in erg/s/cm2/AA.
            Default is False.
        
        sep : [string]
            pandas.read_csv parameter setting the separator in data file.
            Default is ",".
        
        **kwargs
            pandas.read_csv kwargs.
        
        
        Returns
        -------
        _LePhareBase_
        """
        return cls( pandas.read_csv(catin, sep=sep, **kwargs),
                    configfile=configfile, inhz=inhz, dirout=dirout )
                  
    # =============== #
    #  Methods        #
    # =============== #
    # -------- #
    #  SETTER  #
    # -------- #
    def set_data(self, data=None, inhz=False, **extras):
        """
        Set the input data as attribute, ready for LePhare SED fitting.
        
        Parameters
        ----------
        data : [pandas.DataFrame or dict]
            Input data in a compatible LePhare format. Here are the columns style to adopt:
                filtername0, filtername0.err, filtername1, filtername1.err, ... filternameN, filternameN.err, CONTEXT, Z-SPEC, STRING
            where 'filtername{}' is a known filter by LePhare with the format instrument.band, for instance sdss.u or ps1.z, given in flux.
            The 'CONTEXT' sets the filters to use among the provided ones in input data for the SED fitting in each row.
            It is the sum of 2 to the {i}, {i} being the filtername numbers to use in each row.
            For example, using the filters = ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] in data input:
                * ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3 + 2^4 = 31
                * ["sdss.g", "sdss.r", "sdss.i", "sdss.z"]           --> CONTEXT = 2^1 + 2^2 + 2^3 + 2^4       = 30
                * ["sdss.u", "sdss.g", "sdss.r", "sdss.i"]           --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3       = 15
                * ["sdss.u", "sdss.i", "sdss.z"]                     --> CONTEXT = 2^0 + 2^3 + 2^4             = 25
            If 'CONTEXT' column doesn't exist, automatically create one with the maximum CONTEXT value.
            'Z-SPEC' is the spectroscopic redshift, if available.
            'STRING' and all following columns are not used for the SED fitting. Allow the user to store additional informations.
        
        Options
        -------
        inhz : [bool]
            Set to True if the fluxes (and flux uncertainties) are given in erg/s/cm2/Hz ; False means in erg/s/cm2/AA.
            Default is False.

        
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
        filters = self.io.keys_to_filters(data_.columns)
        self.set_filters(filters) 
        if "context" not in data_.columns and "CONTEXT" not in data_.columns:
            warnings.warn("No context key in given data. Adding default context in the dataframe. ")
            data_.insert(self.nfilters*2, "context", self.get_filters_context())
        # // flux units
        if not inhz:
            for filter_ in self.filters:
                data_[filter_]        = tools.flux_aa_to_hz(data_[filter_], self.filter_bandpasses[filter_].wave_eff)
                data_[filter_+".err"] = tools.flux_aa_to_hz(data_[filter_+".err"], self.filter_bandpasses[filter_].wave_eff)
        # // Setting
        self._data = data_
        
    def set_config(self, configfile):
        """
        Handle the configuration file.
        If filters have already been set, handle the filter related parameters in the configuration file.
        
        Parameters
        ----------
        configfile : [string or None]
            Configuration file directory.
            If None, use the default one ($LEPHAREWORK/pylephare/config/default.config).
            Default is None.
        
        
        Returns
        -------
        Void
        """
        from .configparser import ConfigParser
        self._config = ConfigParser.read(configfile)
        self.config.set_fileout(self.io.dirout+"/config")
        
        if self.has_filters():
            self.config.set_value("FILTER_LIST", self.io.get_config_filterlist(self._filters))
            self.config.set_filter_suffix(self._filter_labels)
        
    def set_filters(self, filters):
        """
        Handle with the given filters: set them as attribute and modify the configuration file related parameters.
        
        Parameters
        ----------
        filters : [list(string)]
            List of filters provided in input data.
            Must be known filters by LePhare with the format instrument.band, for instance sdss.u or ps1.z.
        
        
        Returns
        -------
        Void
        """
        _ = super().set_filters(filters)
        if self.has_config():
            self.config.set_value("FILTER_LIST", self.io.get_config_filterlist(self._filters))
            self.config.set_filter_suffix(self._filter_labels)
            
    def set_redshift(self, redshift, index=None):
        """
        Set the redshift in input data either for one row or the whole data.
        
        Parameters
        ----------
        redshift : [float]
            Redshift value to set.
        
        index : [int or None]
            Row index(es) to set the redshift.
            If None, set the given redshift for the whole data table.
            Default is None.
        
        
        Returns
        -------
        Void
        """
        if "zspec" not in self.data.columns and "Z-SPEC" not in self.data.columns:
            self.data.insert(self.nfilters*2 + 1, "zsepc", 0.)
        if index is None:
            self.data["zspec"] = np.array([redshift]*self.ndata)
        else:
            for ii in np.atleast_1d(index):
                self.data.at[ii, "zspec"] = redshift

    def set_context(self, context_value, index=None):
        """ 
        The context defines the filters to use among the provided ones in input data for the SED fitting in each row.
        It is the sum of 2 to the {i}, {i} being the filtername numbers to use in each row.
        For example, using the filters = ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] in data input:
            * ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3 + 2^4 = 31
            * ["sdss.g", "sdss.r", "sdss.i", "sdss.z"]           --> CONTEXT = 2^1 + 2^2 + 2^3 + 2^4       = 30
            * ["sdss.u", "sdss.g", "sdss.r", "sdss.i"]           --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3       = 15
            * ["sdss.u", "sdss.i", "sdss.z"]                     --> CONTEXT = 2^0 + 2^3 + 2^4             = 25
        
        Parameters
        ----------
        context_value : [float]
            Context value to set.
        
        index : [int or list(int) or None]
            Row index to set the context.
            If None, set the given context value for the whole data table.
            Default is None.
        
        
        Returns
        -------
        Void
        """
        if index is None:
            self.data["context"] = np.array([context_value] * self.ndata)
        else:
            for ii in np.atleast_1d(index):
                self.data.at[ii, "context"] = context_value

    def set_photolib_prop(self, gal=True, star=False, qso=False, gallib="BC03", verbose=True):
        """
        Set which kind(s) of SED fitting is(are) done.
        Handle the related parameters in the configuration file.
        
        Parameters
        ----------
        gal : [bool]
            If True, run LePhare on galaxy SED templates.
            Default is True.
        
        star : [bool]
            If True, run LePhare on star SED templates.
            Default is False.
        
        qso : [bool]
            If True, run LePhare on QSO SED templates.
            Default is False.
        
        gallib : [bool]
            If 'gal' is True, define the used galaxy library among the available ones ($LEPHAREDIR/sed/GAL/).
            Default is "BC03".
        
        Options
        -------
        verbose : [bool]
            Print informations.
            Default is True.
        
        
        Returns
        -------
        Void
        """
        if not self.has_config():
            raise AttributeError("No config set.")
        
        self._photolib_set = True
        if verbose:
            print("gal: ",gal, "star: ",star, "qso: ",qso, "; gallib: ",gallib)
        self.config.set_zphotlib(gal=gal, star=star, qso=qso, gallib=gallib)

    def set_dirout(self, dirout):
        """
        Set the output folder directory as attribute.
        
        Parameters
        ----------
        dirout : [string or None]
            Setting the output folder directory.
            If None, create a default directory based on date and time in $LEPHAREWORK/pylephare/.
            Default is None.
        
        
        Returns
        -------
        Void
        """
        self.io.set_dirout(dirout)
        
    # -------- #
    #  GETTER  #
    # -------- #
    def get_filters_context(self, filters=None):
        """
        Return the 'CONTEXT' value in 'data_meas' given the list filters to fit on.
        
        The context defines the filters to use among the provided ones in input data for the SED fitting in each row.
        It is the sum of 2 to the {i}, {i} being the filtername numbers to use in each row.
        For example, using the filters = ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] in data input:
            * ["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"] --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3 + 2^4 = 31
            * ["sdss.g", "sdss.r", "sdss.i", "sdss.z"]           --> CONTEXT = 2^1 + 2^2 + 2^3 + 2^4       = 30
            * ["sdss.u", "sdss.g", "sdss.r", "sdss.i"]           --> CONTEXT = 2^0 + 2^1 + 2^2 + 2^3       = 15
            * ["sdss.u", "sdss.i", "sdss.z"]                     --> CONTEXT = 2^0 + 2^3 + 2^4             = 25

        Parameters
        ----------
        filters : [list(string)]
            List of filters to fit on.
            Must be known filters by LePhare with the format instrument.band, for instance sdss.u or ps1.z.
        
        
        Returns
        -------
        Void
        """
        if filters is None:
            filters = self.filters
            
        idx_cross = np.argwhere(np.in1d(self.filters, filters)).ravel()
        context = np.sum([2**ii for ii in idx_cross])
        return context

    def get_configfile(self, original=False, update=False):
        """
        Return a configuration file directory.
        Either return the original one, or the one at this point of modifications saved in the output folder.
        
        Parameters
        ----------
        original : [bool]
            Set True to return the original configuration file directory ; False to return the one in the output directory.
            Default is False.
        
        Options
        -------
        update : [bool]
            If True, update the configuration file in the output folder with this object configuration parameters.
            Default is False.
        
        
        Returns
        -------
        string
        """
        return self.config._filename if original else self.config.get_fileout(buildit=True, update=update)

    def get_datafile(self, catfile=None):
        """
        Return the input data catalog directory.
        Also set it in the configuration file.
        
        Parameters
        ----------
        catfile : [string or None]
            Input catalog directory.
            Set this in the configuration file.
            If None, return the directory of the data catalog saved in the output directory.
            Default is None.
        
        
        Returns
        -------
        string
        """
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
        """ Test that not void data are loaded """
        return self.data is not None and len(self.data)>0
    
    @property
    def ndata(self):
        """ Number of rows in data """
        if not self.has_data():
            raise AttributeError("No data set yet.")
        return len(self.data)

    # - Config
    @property
    def config(self):
        """ ConfigParser object to interact with the configuration file """
        if not hasattr(self,"_config"):
            self._config = None
        return self._config

    def has_config(self):
        """ Test that a config is loaded """
        return self.config is not None


# ================== #
#                    #
#  MAIN LEPHARE      #
#                    #
# ================== #
class LePhare( _LePhareBase_ ):
    """ _LePhareBase_ + command line tools """
    
    # -------- #
    #  RUNNER  #
    # -------- #
    #
    # - Main
    #
    def run(self, filters=None, contextid=None, dirout=None,
            configfile=None, catinfile=None, originalconfig=False,
            onwhat=["gal","star","qso"], gallib="BC03", update_init=False, verbose=True):
        """
        Execute the LePhare SED fitting ($LEPHAREDIR/source/zphota).
        Return a dictionary containing the directory for the input data catalog, the configuration file, the output folder and the spectra (if any).
        
        Options
        -------
        filters : [list or None]
            List of filters on which to execute the SED fitting, meaning that the context is modified.
            Must be known filters by LePhare with the format instrument.band, for instance sdss.u or ps1.z.
            Go with 'contextid' parameter.
            If None, the context is not changed from what has been set before this point (you can check with {self}.data).
            Default is None.
        
        contextid : [int or list(int) or None]
            Row index(es) to set the new context.
            Go with 'filters' parameter (setting the new context value).
            If None, set the new context value for the whole data table.
            Default is None.
        
        dirout : [string or None]
            Modify the output folder directory.
            If None, use the one already set for this object.
            Default is None.
        
        configfile : [string or None]
            New configuration file directory.
            If None, use the one already set for this object.
            Default is None.
        
        catinfile : [string or None]
            New input data catalog (automatically change it in the configuration file too).
            Go with 'configfile' parameter (must be None).
            If None, use the one already set for this object.
            Default is None.
        
        originalconfig : [bool]
            If True, use the original configuration file ; False acounts for any brought modifications at this point.
            Go with 'configfile' parameter (must be None).
            Default is False.
        
        onwhat: [list or None]
            // ignored if None //
            Defines on which kind of templates to run the SED fitting.
            List which could contain "gal", "star", "qso" (any or combination of).
            For example: ["gal", "qso"]
            Default is ["gal","star","qso"].
        
        gallib: [string]
            If 'gal' in 'onwhat', define the used galaxy library among the available ones ($LEPHAREDIR/sed/GAL/).
            Default is "BC03".
        
        update_init : [bool]
            If True, run the LePhare initialization, even if the templates already exist.
            Default is False.
        
        verbose : [bool]
            Print informations.
            Default is True.
        

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
        self.run_init(update=update_init, configfile=configfile, onwhat=onwhat, gallib=gallib, verbose=verbose)

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
            specfiles = [l for l in os.listdir(".") if l.startswith("Id") and l.endswith(".spec")]
            new_location = [dirout+"/"+l for l in specfiles]
            if verbose:
                print(specfiles)
                print("moved to")
                print(new_location)
            try:
                _ = [os.rename(l, nl) for l, nl in zip(specfiles,new_location)]
            except OSError:
                _ = [shutil.move(l, nl) for l, nl in zip(specfiles,new_location)]
            fileout["spec"] =  new_location            
        else:
            fileout["spec"] = None

        return fileout
        
    #
    # - Secondary
    #
    def build_filter_file(self, update=False, configfile=None, updateconfig=True, verbose=True):
        """
        Create the filter file ($LEPHAREDIR/source/filter).
        
        Options
        -------
        update : [bool]
            If True, create the filter file even if it already exists.
            Default is False.
        
        configfile : [string or None]
            New configuration file directory.
            If None, use the one already set for this object.
            Default is None.
        
        updateconfig : [bool]
            If True, update the configuration file in the output folder with this object configuration parameters.
            Default is False.
        
        verbose : [bool]
            Print informations.
            Default is True.
        
        
        Returns
        -------
        Void
        """
        filter_file = self.config.get_value("FILTER_FILE")
        if configfile is None:
            configfile = self.get_configfile(update=updateconfig)
        if verbose:
            print("filter_file: ", filter_file)
            print("configfile:", configfile)
        if not self.io.is_filt_known(filter_file) or update:
            cmd = "{}/source/filter -c {} -FILTER_FILE {}".format(self.io.LEPHAREDIR, configfile, filter_file)
            try:
                subprocess.run(cmd.split())
            except:
                raise ValueError("LePhareError : unable to run 'filter'.")

    def run_sedtolib(self, update=False, configfile=None, updateconfig=True, photolibprop=None):
        """
        Create the SED templates ($LEPHAREDIR/source/sedtolib)
        
        Options
        -------
        update : [bool]
            If True, run the command even if the templates already exist.
            Default is False.
        
        configfile : [string or None]
            New configuration file directory.
            If None, use the one already set for this object.
            Default is None.
        
        updateconfig : [bool]
            If True, update the configuration file in the output folder with this object configuration parameters.
            Default is False.
        
        photolibprop : [dict]
            // ignired if None //
            Kwargs for the 'set_photolib_prop' method:
                gal : [bool]
                    If True, run LePhare on galaxy SED templates.
                    Default is True.
                star : [bool]
                    If True, run LePhare on star SED templates.
                    Default is False.
                qso : [bool]
                    If True, run LePhare on QSO SED templates.
                    Default is False.
                gallib : [bool]
                    If 'gal' is True, define the used galaxy library among the available ones ($LEPHAREDIR/sed/GAL/).
                    Default is "BC03".
                verbose : [bool]
                    Print informations.
                    Default is True.
            Default is None.
        
        
        Returns
        -------
        Void
        """
        if photolibprop is not None:
            self.set_photolib_prop(**photolibprop)

        if "STAR_LIB" not in self.config.switched_off_keys:
            lib_s_v = self.config.get_value("STAR_LIB")
        else:
            lib_s_v = None
            
        if "QSO_LIB" not in self.config.switched_off_keys:
            lib_q_v = self.config.get_value("QSO_LIB")
        else:
            lib_q_v = None
            
        if "GAL_LIB" not in self.config.switched_off_keys:
            lib_g_v = self.config.get_value("GAL_LIB")
        else:
            lib_g_v =None
        
        if not np.prod([os.path.isfile(self.io.LEPHAREWORK+"/lib_bin/"+elt+".bin")
                        for elt in [lib_s_v, lib_q_v, lib_g_v] if elt is not None]) or update:
            for elt in ["S", "Q", "G"]:
                cmd = "{}/source/sedtolib -t {} -c {}".format(self.io.LEPHAREDIR, elt, self.get_configfile(update=updateconfig) if configfile is None else configfile)
                try:
                    subprocess.run(cmd.split())
                except:
                    raise ValueError("LePhareError : unable to run 'sedtolib' for '{}'.".format(elt))

    def run_mag_star(self, update=False, configfile=None, updateconfig=True, photolibprop=None):
        """
        Create the star magnitude templates ($LEPHAREDIR/source/mag_star).
        
        Options
        -------
        update : [bool]
            If True, run the command even if the templates already exist.
            Default is False.
        
        configfile : [string or None]
            New configuration file directory.
            If None, use the one already set for this object.
            Default is None.
        
        updateconfig : [bool]
            If True, update the configuration file in the output folder with this object configuration parameters.
            Default is False.
        
        photolibprop : [dict]
            // ignired if None //
            Kwargs for the 'set_photolib_prop' method:
                gal : [bool]
                    If True, run LePhare on galaxy SED templates.
                    Default is True.
                star : [bool]
                    If True, run LePhare on star SED templates.
                    Default is False.
                qso : [bool]
                    If True, run LePhare on QSO SED templates.
                    Default is False.
                gallib : [bool]
                    If 'gal' is True, define the used galaxy library among the available ones ($LEPHAREDIR/sed/GAL/).
                    Default is "BC03".
                verbose : [bool]
                    Print informations.
                    Default is True.
            Default is None.
        
        
        Returns
        -------
        Void
        """
        if photolibprop is not None:
            self.set_photolib_prop(**photolibprop)

        if "STAR_LIB_OUT" in self.config.switched_off_keys:
            return None
        
        lib_s_v = self.config.get_value("STAR_LIB_OUT")
        
        if not os.path.isfile(self.io.LEPHAREWORK+"/lib_mag/"+lib_s_v+".bin") or update:
            cmd = "{}/source/mag_star -c {}".format(self.io.LEPHAREDIR, self.get_configfile(update=updateconfig) if configfile is None else configfile)
            try:
                subprocess.run(cmd.split())
            except:
                raise ValueError("LePhareError : unable to run 'mag_star'.")

    def run_mag_gal(self, update=False, configfile=None, updateconfig=True, photolibprop=None):
        """
        Create the galaxy and QSO magnitude templates ($LEPHAREDIR/source/mag_gal).
        
        Options
        -------
        update : [bool]
            If True, run the command even if the templates already exist.
            Default is False.
        
        configfile : [string or None]
            New configuration file directory.
            If None, use the one already set for this object.
            Default is None.
        
        updateconfig : [bool]
            If True, update the configuration file in the output folder with this object configuration parameters.
            Default is False.
        
        photolibprop : [dict]
            // ignired if None //
            Kwargs for the 'set_photolib_prop' method:
                gal : [bool]
                    If True, run LePhare on galaxy SED templates.
                    Default is True.
                star : [bool]
                    If True, run LePhare on star SED templates.
                    Default is False.
                qso : [bool]
                    If True, run LePhare on QSO SED templates.
                    Default is False.
                gallib : [bool]
                    If 'gal' is True, define the used galaxy library among the available ones ($LEPHAREDIR/sed/GAL/).
                    Default is "BC03".
                verbose : [bool]
                    Print informations.
                    Default is True.
            Default is None.
        
        
        Returns
        -------
        Void
        """
        if photolibprop is not None:
            self.set_photolib_prop(**photolibprop)
         
        if "QSO_LIB_OUT" not in self.config.switched_off_keys:
            lib_q_v = self.config.get_value("QSO_LIB_OUT")
        else:
            lib_q_v = None
            
        if "GAL_LIB_OUT" not in self.config.switched_off_keys:
            lib_g_v = self.config.get_value("GAL_LIB_OUT")
        else:
            lib_g_v =None
        
        if not np.prod([os.path.isfile(self.io.LEPHAREWORK+"/lib_mag/"+elt+".bin") for elt in [lib_q_v, lib_g_v] if elt is not None]) or update:
            for elt in ["Q", "G"]:
                cmd = "{}/source/mag_gal -t {} -c {}".format(self.io.LEPHAREDIR, elt, self.get_configfile(update=updateconfig) if configfile is None else configfile)
                try:
                    subprocess.run(cmd.split())
                except:
                    raise ValueError("LePhareError : unable to run 'mag_gal' for '{}'.".format(elt))

    def run_init(self, update=False, configfile=None, updateconfig=True, onwhat=None, gallib="BC03", verbose=True):
        """
        Run the whole LePhare initialization.
        
        Options
        -------
        update : [bool]
            If True, run the command even if the templates already exist.
            Default is False.
        
        configfile : [string or None]
            New configuration file directory.
            If None, use the one already set for this object.
            Default is None.
        
        updateconfig : [bool]
            If True, update the configuration file in the output folder with this object configuration parameters.
            Default is False.
        
        onwhat: [list or None]
            // ignored if None //
            Defines on which kind of templates to run the SED fitting.
            List which could contain "gal", "star", "qso" (any or combination of).
            For example: ["gal", "qso"]
            Default is ["gal","star","qso"].
        
        gallib: [string]
            If 'gal' in 'onwhat', define the used galaxy library among the available ones ($LEPHAREDIR/sed/GAL/).
            Default is "BC03".
        
        verbose : [bool]
            Print informations.
            Default is True.
        
        
        Returns
        -------
        Void
        """
        if onwhat is not None:
            photolibprop = {**{"gal":False, "qso":False, "star":False},**{k:True for k in np.atleast_1d(onwhat)}}
            self.set_photolib_prop(gallib=gallib, verbose=verbose, **photolibprop)
            
        prop = dict(update=update, configfile=configfile, updateconfig=updateconfig)
        self.build_filter_file(verbose=verbose, **prop)
        self.run_sedtolib(**prop)
        self.run_mag_star(**prop)
        self.run_mag_gal(**prop)
        
# ================== #
#                    #
#  LEPHARE RESULTS   #
#                    #
# ================== #
