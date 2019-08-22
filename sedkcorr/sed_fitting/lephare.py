
import numpy as np
import pandas
import os
import subprocess
import pkg_resources

from propobject import BaseObject


class LePhareSEDFitter( BaseObject ):
    """
    
    """

    PROPERTIES         = ["input_param", "output_param"]
    SIDE_PROPERTIES    = ["path_config", "path_init", "path_results"]
    DERIVED_PROPERTIES = []
    
    INPUT_PARAM = {"STAR_SED":"$LEPHAREDIR/sed/STAR/STAR_MOD.list",
                   "STAR_FSCALE":3.432E-09,
                   "STAR_LIB":"LIB_STAR",
                   "QSO_SED":"$LEPHAREDIR/sed/QSO/QSO_MOD.list",
                   "QSO_FSCALE":1.,
                   "QSO_LIB":"LIB_QSO",
                   "GAL_SED":"$LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_MOD.list",
                   "GAL_FSCALE":1.,
                   "GAL_LIB":"LIB_BC03",
                   "SEL_AGE":"$LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_AGE.list",
                   "AGE_RANGE":(0.,13.55e9)
                   "end":None}
    
    PATH_LEPHAREDIR = os.path.expanduser(os.getenv("LEPHAREDIR"))

    def __init__(self, **kwargs):
        """
        
        """
        self.set_data(**kwargs)

    def set_data(self, **kwargs):
        """
        
        """
        self.set_path({"config":None, "init":None, "results":None})
        for key, value in self.INPUT_PARAM.items():
            self._write_param_(key, value)
    
    def set_path(self, which={"config":None, "init":None, "results":None}):
        """
        
        """
        for key, value in which.items():
            path = (pkg_resources.resource_filename(__name__, key+"/") + "/") if v is None else value
            self._side_properties["path_"+key] = path
    
    def _write_param_(self, line, new_param_value):
        """
        
        """
        if type(line) == int:
            idx_line = line
        elif type(line) == str:
            with open(self.path_config+"lephare_zphot_input.para", "r") as file:
                file_buf = [line for line in file]
            for ii, line in enumerate(file_buf):
                if line.split(" ")[0] == line:
                    idx_line = ii
                    break
        else:
            raise ValueError("'line' must be either integer (line index) or string (first word of the line)")
        
        if type(new_param_value) in [int, float]:
            new_param_value = str(new_param_value)
        elif type(new_param_value) in [tuple, list]:
            new_param_value = ",".join(new_param_value)
        
        line_splited = file_buf[idx_line].split()
        line_splited[1] = new_param_value
        file_buf[idx_line] = " ".join(line_splited) + "\n"

        with open(self.path_config+"lephare_zphot_input.para", "w") as file:
            for line in file_buf:
                file.write(line)
    
    def _set_data_path_(self, data_path):
        """
        
        """
        self._write_para_("CAT_IN", data_path)
    
    def run_sedtolib(self, path_config=None):
        """
        
        """
        if path_config is not None:
            self.set_path({"path_config":path_config})
        
        for elt in ["S", "Q", "G"]:
            cmd = "{}/source/sedtolib -t {} -c {}/lephare_zphot_input.para".format(self.PATH_LEPHAREDIR, elt, self.path_config)
            subprocess.run(cmd.split())
    
    def run_mag_star(self, path_config=None, path_init=None):
        """
        
        """
        if path_config is not None:
            self.set_path({"path_config":path_config})
        if path_init is not None:
            self.set_path({"path_init":path_init})
        
        os.chdir(self.path_init)
        cmd = "{}/source/mag_star -c {}/lephare_zphot_input.para".format(self.PATH_LEPHAREDIR, self.path_config)
        subprocess.run(cmd.split())
    
    def run_mag_gal(self, path_config=None, path_init=None):
        """
        
        """
        if path_config is not None:
            self.set_path({"path_config":path_config})
        if path_init is not None:
            self.set_path({"path_init":path_init})
        
        os.chdir(self.path_init)
        for elt in ["Q", "G"]:
            cmd = "{}/source/mag_gal -t {} -c {}/lephare_zphot_input.para".format(self.PATH_LEPHAREDIR, elt, self.path_config)
            subprocess.run(cmd.split())
    
    def run_zphota(self, path_config=None, path_results=None):
        """
        
        """
        if path_config is not None:
            self.set_path({"path_config":path_config})
        if path_results is not None:
            self.set_path({"path_results":path_results})
        
        os.chdir(self.path_results)
        cmd = "{}/source/zphota -c {}/lephare_zphot_input.para".format(self.PATH_LEPHAREDIR, self.path_config)
        subprocess.run(cmd.split())
        return
    
    

    #-------------------#
    #   Properties      #
    #-------------------#
    @property
    def input_param(self):
        """  """
        return self._properties["input_param"]
    
    @property
    def output_param(self):
        """  """
        return self._properties["output_param"]
    
    @property
    def path_config(self):
        """  """
        return self._side_properties["path_config"]
    
    @property
    def path_init(self):
        """  """
        return self._side_properties["path_init"]
    
    @property
    def path_results(self):
        """  """
        return self._side_properties["path_results"]


