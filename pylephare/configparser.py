
import os
import numpy as np
import warnings
from . import io

class ConfigParser( object ):
    """ """
    _TYPE = "_configholder_"
    # ================ #
    #    Methods       #
    # ================ #
    def __init__(self, filename=None, fileconfout=None):
        """ """
        self._init_config_()
        if filename is None:
            filename = io.DEFAULTCONFIG
            warnings.warn("using default config file")
        if fileconfout is None:
            fileconfout = io.DEFAULTCONFIG_OUT
            
        self.load_config(filename, fileconfout=fileconfout)

    def _init_config_(self):
        """ """
        self._config= {}
        
    @classmethod
    def read(cls, filename):
        """ """
        return cls(filename)

    def writeto(self, filename=None):
        """ """
        if filename is None:
            if self.fileout is not None:
                filename = self.fileout
            else:
                raise IOError("filename is None and no self.fileout set.")
            
        with open(filename, 'w') as f:
            f.writelines("\n".join(self.get_config("configfile")))

    # -------- #
    #  LOADDER #
    # -------- #
    def load_config(self, filename, fileconfout=None):
        """ """
        self._filename = filename
        self._init_config = open(filename).read().splitlines()
        metaconfig = [l for l in self._init_config
                      if not (l.startswith("#") or l.startswith("\t") or l.startswith(" ")) and len(l)>0]

        for k in metaconfig:
            key, *value = k.split()
            if len(value)==1:
                self.set_value(key, value[0],None)
                
            elif len(value)>1:
                if value[1]=="#":
                    self.set_value(key, value[0], " ".join(value[2:]))
                else:
                    raise IOError("cannot parse the line %s"%k)
            else:
                raise IOError("cannot parse the line %s"%k)
        if fileconfout is not None:
            self.set_value("PARA_OUT", fileconfout)
            
    def _build_new_config_(self):
        """ """
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
            return newconfig
        else:
            # new
            return self.get_config("array")


    # -------- #
    #  SETTER  #
    # -------- #
    def set_fileout(self, fileout, builddir=True):
        """ """
        if builddir:
            dir_ = os.path.dirname(fileout)
            if not os.path.isdir(dir_):
                os.makedirs(dir_)
                
        self._fileout = fileout

        
    def set_value(self, key, value, comments=None):
        """ Mother setter methods. All data set this way will be recorded with using self.writeto() """
        self._config[key.upper()] = {"value":value,
                                     "comments":comments
                                    }            
        
    def switch_off_key(self, key):
        """ """
        if key not in self.switched_off_keys:
            self.switched_off_keys.append(key)
            self._config["# "+key] = self._config.pop(key)
        
    def switch_on_key(self, key):
        """ """
        if key in self.switched_off_keys:
            self.switched_off_keys.pop(self.switched_off_keys.index(key))
            self._config[key] = self._config.pop("# "+key)

    def set_filter_suffix(self, suffix ):
        """ 
        Add the suffix to the {}_LIB, {}_LIB_IN {}_LIB_OUT values
        {} = {STAR, QSO or GAL}

        It also sets FILTER_FILE = suffix.filt
        """
        print('running set_filter_suffix')
        for key in ["STAR","QSO","GAL"]:
            self.set_value("%s_LIB"%key, "LIB_%s_%s"%(key,suffix))
            self.set_value("%s_LIB_IN"%key, "LIB_%s_%s"%(key,suffix))
            self.set_value("%s_LIB_OUT"%key, "%s_%s"%(key,suffix))

        self.set_value("FILTER_FILE", suffix+".filt")
        self.set_zphotlib()
        
    def set_zphotlib(self, gal=True, star=False, qso=False, gallib="BC03"):
        """ """
        print("calling set_zphotlib")
        zphoto=""
        if not star:
            self.switch_off_key("STAR_LIB")
            self.switch_off_key("STAR_LIB_IN")
            self.switch_off_key("STAR_LIB_OUT")
            starkey = []
        else:
            starkey = [self.get_value("STAR_LIB_OUT")]
            
        if not qso:
            self.switch_off_key("QSO_LIB")
            self.switch_off_key("QSO_LIB_IN")
            self.switch_off_key("QSO_LIB_OUT")
            qsokey = []
        else:
            qsokey = [self.get_value("QSO_LIB_OUT")]

        if not gal:
            self.switch_off_key("GAL_LIB")
            self.switch_off_key("GAL_LIB_IN")
            self.switch_off_key("GAL_LIB_OUT")
            galkey = []
        else:
            for key in ["GAL_LIB","GAL_LIB_IN","GAL_LIB_OUT"]:
                self.set_value(key, self.get_value(key).replace("GAL",gallib) )
            galkey = [self.get_value("GAL_LIB_OUT")]

        self.set_value("ZPHOTLIB", ','.join(starkey+qsokey+galkey))
        return self.get_value("ZPHOTLIB")

    def set_intrinsic_error(self, err_scale):
        """ """
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
        """ """
        if key not in self._config:
            raise ValueError("%s not in self.config"%key)
        
        return self._config[key]["value"]
    
    def get_comments(self, key):
        """ """
        if key not in self._config:
            raise ValueError("%s not in self.config"%key)
        return self._config[key]["comments"]

    def get_config_lines(self, key):
        """ """
        return " ".join([key, self._config[key]["value"], "#", self._config[key]["comments"]]
                            if self._config[key]["comments"] is not None else
                        [key, self._config[key]["value"]]
                        )

    def get_config(self, astype="dict"):
        """ """
        if astype == "dict":
            return self.config
        elif astype.lower() == "dataframe":
            import pandas
            return pandas.DataFrame(self.config).T
        
        elif astype.lower() == "array":
            return [self.get_config_lines(k_) for k_ in self.config.keys()]
        
        elif astype.lower() in ["file", "configfile"]:
            return self._build_new_config_()
        
        elif astype.lower() in ["original","input","inputfile"]:
            if hasattr(self,"_init_config"):
                return self._init_config
            else:
                raise ValueError("No original config file")
            
        raise ValueError("astype '%s' not supported. Could be dict, dataframe, array, configfile, inputfile"%astype)

    def get_fileout(self, buildit=True, update=False):
        """ """
        if self.has_fileout():
            if update or (not os.path.isfile(self.fileout) and buildit):
                self.writeto()
                
        return self.fileout

    def get_filters(self, name=False):
        """ """
        filtfile_list = self.get_value("FILTER_LIST").split(",")
        if not name:
            return filtfile_list
        return [io.filterfile_to_filtername(filt) for filt in filtfile_list]
        
            
    def get_catin_columns(self):
        """ """
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
