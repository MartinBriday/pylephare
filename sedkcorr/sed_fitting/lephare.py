
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
                raise ValueError("{} is not a known parameter in the parameters file.".format(line))
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
    
    def _get_cleared_line_(self, idx_line):
        """
        
        """
        flag_comment = False
        splitted_line = file_buf[idx_line].split()
        if splitted_line[0] == "#":
            splitted_line.pop(0)
            flag_comment = True
        elif len(splitted_line[0])>1 and splitted_line[0][0] == "#":
            splitted_line[0] = splitted_line[0].replace("#", "")
            flag_comment = True
        return splitted_line, flag_comment
    
    def _get_file_lines(self, file):
        """
        
        """
        with open(file, "r") as file:
            file_buf = [line for line in file]
        return file_buf
    
    def change_input_param(self, line, new_param_value, force_comment=False):
        """
        
        """
        file_buf = self._get_file_lines(self.input_param_file)
        idx_line = self._get_idx_line_(file_buf, line)
        new_param_value = self._get_new_param_value_(new_param_value)
        
        splitted_line, _ = self._get_cleared_line_(idx_line)
        #idx_param = 2 if splitted_line[0][0] == "#" else 1
        splitted_line[1] = new_param_value
        if force_comment:
            splitted_line.insert(0, "#")
        file_buf[idx_line] = " ".join(splitted_line) + "\n"

        with open(self.input_param_file, "w") as file:
            for line in file_buf:
                file.write(line)

    def change_output_param(self, line, force_comment):
        """
        
        """
        file_buf = self._get_file_lines(self.output_param_file)
        idx_line = self._get_idx_line_(file_buf, line)
        splitted_line, _ = self._get_cleared_line_(idx_line)
        
        if force_comment:
            splitted_line.insert(0, "#")
        file_buf[idx_line] = " ".join(splitted_line) + "\n"
        
        with open(self.input_param_file, "w") as file:
            for line in file_buf:
                file.write(line)

    def _print_param_details_(self, param):
        """
        
        """
        if param in self.INPUT_PARAM:
            which = "input"
        elif param in self.OUTPUT_PARAM:
            which = "output"
        else:
            raise ValueError("'{}' neither is in input parameter list nor output parameter list.".format(param))
        file_buf = self._get_file_lines(self._properties[which+"_param_file"])
        idx_line = self._get_idx_line_(file_buf, param)
        splitted_line, flag_comment = self._get_cleared_line_(idx_line)
        
        txt_param = "{} : {} (comented? : {})".format(param, splitted_line[1] if which=="input" else "", flag_comment)
        
        print(txt_param)
        return

    def describe_params(self, which="input"):
        """
        
        """
        if which == "input":
            list_param = self.INPUT_PARAM
        elif which == "output":
            list_param = self.OUTPUT_PARAM
        else:
            raise ValueError("'which' must be either 'input' or 'output'".)
        for _param in list_param:
            self._print_param_details_(_param)
    
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


