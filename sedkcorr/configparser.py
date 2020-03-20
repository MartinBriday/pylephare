import numpy as np
import warnings

class ConfigParser( object ):
    """ """
    _TYPE = "_configholder_"
    # ================ #
    #    Methods       #
    # ================ #
    def __init__(self, filename):
        """ """
        self._init_config_()
        self.load_config(filename)

    def _init_config_(self):
        """ """
        self._config= {}
        
    @classmethod
    def read(cls, filename):
        """ """
        return cls(filename)
    # -------- #
    #  LOADDER #
    # -------- #
    def load_config(self, filename):
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
                    newconfig.append(self.get_config_line(key))
                    set_key.append(key)

            keys = np.asarray(list(self.config.keys()))
            not_set_key = keys[~np.in1d(keys, set_key)]
            
            if len(not_set_key)>0:
                newconfig.append("# == NEW KEYS == #")
                for key in not_set_key:
                    newconfig.append(self.get_config_line(key))
            return newconfig
        else:
            # new
            return self.get_config("array")
    # -------- #
    #  SETTER  #
    # -------- #
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
        
    # -------- #
    #  GETTER  #
    # -------- #
    def get_value(self, key):
        """ """
        if key not in self._config:
            raise ValueError("%s not in self.config"%key)
        
        return self._config[key]["value"]
    
    def get_comments(self, key, default=np.NaN):
        """ """
        if key not in self._config:
            raise ValueError("%s not in self.config"%key)
        return self._config[key]["comments"]

    def get_config_line(self, key):
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
            return [self.get_config_line(k_) for k_ in self.config.keys()]
        elif astype.lower() in ["file", "configfile"]:
            return self._build_new_config_()
        elif astype.lower() in ["original","input","inputfile"]:
            if hasattr(self,"_init_config"):
                return self._init_config
            else:
                raise ValueError("No original config file")
            
        raise ValueError("astype not value (%s). Could be dict, dataframe, array, configfile, inputfile"%astype)



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
