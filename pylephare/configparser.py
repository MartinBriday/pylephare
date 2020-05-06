
import os
import numpy as np
import warnings
from . import io

class ConfigParser( object ):
    """ Handle configuration files and their parameters. """
    
    _TYPE = "_configholder_"
    # ================ #
    #    Methods       #
    # ================ #
    def __init__(self, filename=None, fileconfout=None):
        """
        Class builder.
        
        Parameters
        ----------
        filename : [string or None]
            Configuration file (input) directory.
            If None, the configuration file will be the default one ($LEPHAREWORK/pylephare/config/default.config).
            Default is None.
        
        fileconfout : [string or None]
            Output configuration file.
            If None, the output configuration file will be the default one ($LEPHAREWORK/pylephare/config/default_output.config).
            Default is None.
        
        
        Returns
        -------
        Void
        """
        self._init_config_()
        if filename is None:
            filename = io.DEFAULTCONFIG
            warnings.warn("Using default configuration file.")
        if fileconfout is None:
            fileconfout = io.DEFAULTCONFIG_OUT
            
        self.load_config(filename, fileconfout=fileconfout)

    def _init_config_(self):
        """
        Reset the configuration parameters attribute.
        
        
        Returns
        -------
        Void
        """
        self._config= {}
        
    @classmethod
    def read(cls, filename):
        """
        Read a given configuration file and build a ConfigParser object.
        
        Parameters
        ----------
        filename : [string or None]
            Configuration file (input) directory.
        
        
        Returns
        -------
        ConfigParser
        """
        return cls(filename)

    def writeto(self, filename=None):
        """
        Write the configuration parameters into the given configuration file.
        
        Parameters
        ----------
        filename : [string or None]
            Configuration file (input) directory.
            If None, the directory will be the already set 'fileout'.
            Default is None
        
        
        Returns
        -------
        
        """
        if filename is None:
            if self.fileout is not None:
                filename = self.fileout
            else:
                raise IOError("filename is None and no self.fileout set.")
            
        with open(filename, "w") as f:
            f.writelines("\n".join(self.get_config("configfile")))

    # -------- #
    #  LOADDER #
    # -------- #
    def load_config(self, filename, fileconfout=None):
        """
        Handle the original configuration file directories (input and output).
        Load the input configuration parameters from the file into attribute.
        
        Parameters
        ----------
        filename : [string or None]
            Configuration file (input) directory.
        
        fileconfout : [string or None]
            // ignored if None //
            Output configuration file.
        
        
        Returns
        -------
        Void
        """
        self._filename = filename
        self._init_config = open(filename).read().splitlines()
        metaconfig = [l for l in self._init_config
                      if not (l.startswith("#") or l.startswith("\t") or l.startswith(" ")) and len(l)>0]

        for k in metaconfig:
            key, *value = k.split()
            if len(value)==1:
                self.set_value(key, value[0], None)
                
            elif len(value)>1:
                if value[1]=="#":
                    self.set_value(key, value[0], " ".join(value[2:]))
                else:
                    raise IOError("Cannot parse the line %s"%k)
            else:
                raise IOError("cannot parse the line %s"%k)
        if fileconfout is not None:
            self.set_value("PARA_OUT", fileconfout)
            
    def _build_new_config_(self):
        """
        Build and return the
        
        
        Returns
        -------
        np.array
        """
        if hasattr(self,"_init_config"):
            # Has been loaded
            newconfig = []
            set_key = []
            for l in self._init_config:
                if ((l.startswith("#") or l.startswith("\t") or l.startswith(" ")) or len(l)==0):
                    newconfig.append(l)
                else:
                    key = l.split()[0]
                    if key in self.switched_off_keys:
                        key = "# "+key
                    newconfig.append(self.get_config_lines(key))
                    set_key.append(key)

            keys = np.asarray(list(self.config.keys()))
            not_set_key = keys[~np.in1d(keys, set_key)]
            
            if len(not_set_key)>0:
                newconfig.append("# == NEW KEYS == #")
                for key in not_set_key:
                    newconfig.append(self.get_config_lines(key))
            return np.array(newconfig)
        else:
            # new
            return self.get_config("array")


    # -------- #
    #  SETTER  #
    # -------- #
    def set_fileout(self, fileout, builddir=True):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if builddir:
            dir_ = os.path.dirname(fileout)
            if not os.path.isdir(dir_):
                os.makedirs(dir_)
                
        self._fileout = fileout

        
    def set_value(self, key, value, comments=None):
        """
        Mother setter methods. All data set this way will be recorded with using self.writeto()
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        self._config[key.upper()] = {"value":value,
                                     "comments":comments
                                    }            
        
    def switch_off_key(self, key):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if key not in self.switched_off_keys:
            self._switched_off_keys.append(key)
            self._config["# "+key] = self._config.pop(key)
        
    def switch_on_key(self, key):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if key in self.switched_off_keys:
            self._switched_off_keys.pop(self._switched_off_keys.index(key))
            self._config[key] = self._config.pop("# "+key)

    def set_filter_suffix(self, suffix ):
        """ 
        Add the suffix to the {}_LIB, {}_LIB_IN {}_LIB_OUT values
        {} = {STAR, QSO or GAL}

        It also sets FILTER_FILE = suffix.filt
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        print("Running set_filter_suffix.")
        for key in ["STAR","QSO","GAL"]:
            self.set_value("%s_LIB"%key, "LIB_%s_%s"%(key,suffix))
            self.set_value("%s_LIB_IN"%key, "LIB_%s_%s"%(key,suffix))
            self.set_value("%s_LIB_OUT"%key, "%s_%s"%(key,suffix))

        self.set_value("FILTER_FILE", suffix+".filt")
        self.set_zphotlib()
        
    def set_zphotlib(self, gal=True, star=False, qso=False, gallib="BC03"):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        def switch_lib(self, lib, button):
            _button = "on" if button else "off"
            for ii in ["_LIB", "_LIB_IN", "_LIB_OUT"]:
                eval("""self.switch_{}_key("{}")""".format(_button, lib+ii))
            return eval("""self.get_value("{}_LIB_OUT")""".format(lib)) if button else None
            
        print("calling set_zphotlib")
        
        if gal:
            for key in ["GAL_LIB","GAL_LIB_IN","GAL_LIB_OUT"]:
                self.set_value(key, self.get_value(key).replace("GAL",gallib) )
        
        _zphotlib = []
        for _lib in ["STAR", "QSO", "GAL"]:
            _lib_out = switch_lib(self, lib=_lib, button=eval(_lib.lower()))
            if _lib_out is not None:
                _zphotlib.append(_lib_out)

        self.set_value("ZPHOTLIB", ",".join(_zphotlib))
        return self.get_value("ZPHOTLIB")

    def set_intrinsic_error(self, err_scale):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if err_scale is None:
            self.switch_off_key("ERR_SCALE")
            return
        elif len(err_scale) != len(self.get_filters()):
            raise ValueError(" you must provide the number of intrinsic magnitude as you have filters (%d vs. %d)"%(len(err_scale), len(self.get_filters())))
        
        self.switch_on_key("ERR_SCALE") # ignored if already on
        self.set_value("ERR_SCALE",",".join(["{}".format(l) for l in err_scale]))
        
    # -------- #
    #  GETTER  #
    # -------- #
    def get_value(self, key):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if key not in self._config:
            raise ValueError("%s not in self.config"%key)
        
        return self._config[key]["value"]
    
    def get_comments(self, key):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if key not in self._config:
            raise ValueError("%s not in self.config"%key)
        return self._config[key]["comments"]

    def get_config_lines(self, key):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        return " ".join([key, self._config[key]["value"], "#", self._config[key]["comments"]]
                            if self._config[key]["comments"] is not None else
                        [key, self._config[key]["value"]]
                        )

    def get_config(self, astype="dict"):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if astype == "dict":
            return self.config
        elif astype.lower() == "dataframe":
            import pandas
            return pandas.DataFrame(self.config).T
        
        elif astype.lower() == "array":
            return np.array([self.get_config_lines(k_) for k_ in self.config.keys()])
        
        elif astype.lower() in ["file", "configfile"]:
            return self._build_new_config_()
        
        elif astype.lower() in ["original","input","inputfile"]:
            if hasattr(self,"_init_config"):
                return self._init_config
            else:
                raise ValueError("No original config file")
            
        raise ValueError("astype '%s' not supported. Could be dict, dataframe, array, configfile, inputfile"%astype)

    def get_fileout(self, buildit=True, update=False):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        if self.has_fileout():
            if update or (not os.path.isfile(self.fileout) and buildit):
                self.writeto()
                
        return self.fileout

    def get_filters(self, name=False):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        filtfile_list = self.get_value("FILTER_LIST").split(",")
        if not name:
            return filtfile_list
        return [io.filterfile_to_filtername(filt) for filt in filtfile_list]
        
            
    def get_catin_columns(self):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        catformat = self.get_value("CAT_FMT")
        cattype   = self.get_value("CAT_TYPE")
        filters   = self.get_filters(name=True)
        
        if catformat == "MMEE":
            names = list(filters)+[f+".err" for f in filters]
        elif catformat == "MEME":
            names = list(np.concatenate([[f,f+".err"] for f in filters]))
        else:
            raise ValueError("Cannot Parse the CAT_FMT (%s) from the configfile"%catformat)
        
        if cattype == "LONG":
            names += ["context","z-spec", "string"]
        elif cattype == "SHORT":
            names += ["context"]
        else:
            raise ValueError("Cannot Parse the CAT_TYPE (%s) from the configfile"%cattype)
        
        return names
    # ================ #
    #  Properties      #
    # ================ #
    @property
    def switched_off_keys(self):
        """ """
        if not hasattr(self, "_switched_off_keys"):
            self._switched_off_keys = []
            
        return self._switched_off_keys
    @property
    def config(self):
        """ """
        if not hasattr(self,"_config") or self._config is None:
            self._config = {}
            
        return self._config

    @property
    def fileout(self):
        """ """
        if not hasattr(self, "_fileout"):
            self._fileout = None
        return self._fileout

    def has_fileout(self):
        """ """
        return self.fileout is not None
