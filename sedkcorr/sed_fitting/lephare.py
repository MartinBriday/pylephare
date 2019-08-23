
import numpy as np
import pandas
import os
import subprocess
import pkg_resources

from propobject import BaseObject


class LePhareSEDFitter( BaseObject ):
    """
    
    """

    PROPERTIES         = ["input_param_file", "output_param_file"]
    SIDE_PROPERTIES    = ["path_results"]
    DERIVED_PROPERTIES = []
    
    PATH_LEPHAREDIR = os.path.expanduser(os.getenv("LEPHAREDIR"))

    def __init__(self, data_path=None, **kwargs):
        """
        
        """
        self.set_data(data_path, **kwargs)

    def set_data(self, data_path=None, input_param_file=None, output_param_file=None, path_results=None, **kwargs):
        """
        
        """
        self._properties["input_param_file"] = pkg_resources.resource_filename(__name__, "config/") + "/lephare_zphot_input.para" \
                                               if input_param_file is None else input_param_file
        self._properties["input_param_file"] = pkg_resources.resource_filename(__name__, "config/") + "/lephare_zphot_output.para" \
                                               if input_param_file is None else input_param_file
        self._side_properties["path_results"] = pkg_resources.resource_filename(__name__, "results/") + "/" \
                                                if path_results is None else path_results
        
        if data_path is not None:
            self.change_input_param("CAT_IN", data_path)

    def _get_idx_line_(self, file, line):
        """
        
        """
        if type(line) == int:
            idx_line = line
        elif type(line) == str:
            flag_no_match = True
            for ii, line in enumerate(file_buf):
                splitted_line = line.split()
                if line in splitted_line[0] or splitted_line[1] == line:
                    idx_line = ii
                    flag_no_match = False
                    break
            if flag_no_match:
                raise ValueError("{} is not a known parameter in the input parameters file.".format(line))
        else:
            raise ValueError("'line' must be either integer (line index) or string (first word of the line)")
        return idx_line
    
    def _get_new_param_value_(self, new_param_value):
        """
        
        """
        if new_param_value in ["YES", "yes", "Yes", True]:
            new_param_value = "YES"
        elif new_param_value in ["NO", "no", "No", False]:
            new_param_value = "NO"
        elif type(new_param_value) in [int, float]:
            new_param_value = str(new_param_value)
        elif type(new_param_value) in [tuple, list]:
            new_param_value = ",".join([str(elt) for elt in new_param_value])
        return new_param_value
    
    def change_input_param(self, line, new_param_value, force_comment=False):
        """
        
        """
        with open(self.input_param_file, "r") as file:
            file_buf = [line for line in file]
        idx_line = self._get_idx_line_(file_buf, line)
        new_param_value = self._get_new_param_value_(new_param_value)
        
        splitted_line = file_buf[idx_line].split()
        idx_param = 2 if splitted_line[0][0] == "#" else 1
        splitted_line[idx_param] = new_param_value
        if force_comment:
            splitted_line.insert(0, "#")
        file_buf[idx_line] = " ".join(splitted_line) + "\n"

        with open(self.input_param_file, "w") as file:
            for line in file_buf:
                file.write(line)
    
    def run_sedtolib(self, input_param_file=None):
        """
        
        """
        if path_config is not None:
            self._properties["input_param_file"] = input_param_file
        
        for elt in ["S", "Q", "G"]:
            cmd = "{}/source/sedtolib -t {} -c {}".format(self.PATH_LEPHAREDIR, elt, self.input_param_file)
            subprocess.run(cmd.split())
    
    def run_mag_star(self, input_param_file=None):
        """
        
        """
        if path_config is not None:
            self._properties["input_param_file"] = input_param_file
        
        cmd = "{}/source/mag_star -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
        subprocess.run(cmd.split())
    
    def run_mag_gal(self, input_param_file=None):
        """
        
        """
        if path_config is not None:
            self._properties["input_param_file"] = input_param_file
        
        for elt in ["Q", "G"]:
            cmd = "{}/source/mag_gal -t {} -c {}".format(self.PATH_LEPHAREDIR, elt, self.input_param_file)
            subprocess.run(cmd.split())
    
    def run_zphota(self, input_param_file=None, path_results=None):
        """
        
        """
        if path_config is not None:
            self._properties["input_param_file"] = input_param_file
        if path_results is not None:
            self._side_properties["path_results"] = path_results
        
        os.chdir(self.path_results)
        cmd = "{}/source/zphota -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
        subprocess.run(cmd.split())
        return
    
    

    #-------------------#
    #   Properties      #
    #-------------------#
    @property
    def input_param_file(self):
        """  """
        return self._properties["input_param_file"]

    @property
    def output_param_file(self):
        """  """
        return self._properties["output_param_file"]
    
    @property
    def path_results(self):
        """  """
        return self._side_properties["path_results"]


