
import numpy as np
import pandas
import os
import subprocess
import pkg_resources

from propobject import BaseObject


class LePhareSEDFitter( BaseObject ):
    """
    
    """

    PROPERTIES         = ["data", "input_param_file", "output_param_file"]
    SIDE_PROPERTIES    = ["path_results"]
    DERIVED_PROPERTIES = []
    
    INPUT_PARAM = {#CREATION OF LIBRARIES FROM SEDs List
                   "STAR_SED":"$LEPHAREDIR/sed/STAR/STAR_MOD.list",
                   "STAR_FSCALE":3.432E-09,
                   "STAR_LIB":"LIB_STAR",
                   "QSO_SED":"$LEPHAREDIR/sed/QSO/QSO_MOD.list",
                   "QSO_FSCALE":1.,
                   "QSO_LIB":"LIB_QSO",
                   "GAL_SED":"$LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_MOD.list",
                   "GAL_FSCALE":1.,
                   "GAL_LIB":"LIB_BC03",
                   "SEL_AGE":"$LEPHAREDIR/sed/GAL/BC03_CHAB/BC03_AGE.list",
                   "AGE_RANGE":(0.,13.55e9),
                   #FILTERS
                   "FILTER_LIST":["sdss/up.pb","sdss/gp.pb","sdss/rp.pb","sdss/ip.pb","sdss/zp.pb"],
                   "TRANS_TYPE":0,
                   "FILTER_CALIB":0,
                   "FILTER_FILE":"sdss.filt",
                   #THEORETICAL MAGNITUDES
                   "STAR_LIB_IN":"LIB_STAR",
                   "STAR_LIB_OUT":"STAR",
                   "QSO_LIB_IN":"LIB_QSO",
                   "QSO_LIB_OUT":"QSO",
                   "GAL_LIB_IN":"LIB_BC03",
                   "GAL_LIB_OUT":"BC03",
                   "MAGTYPE":"AB",
                   "Z_STEP":(0.004,0.20,0.1),
                   "COSMOLOGY":(70,0.3,0.7),
                   "MOD_EXTINC":(1,27),
                   "EXTINC_LAW":"LMC_Fitzpatrick.dat",
                   "EB_V":(0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8),
                   "EM_LINES":"YES",
                   "Z_FORM":(8,7,6,5,4,3),
                   "LIB_ASCII":"YES",
                   #PHOTOMETRIC REDSHIFTS
                   "CAT_IN":"",
                   "INP_TYPE":"M",
                   "CAT_MAG":"AB",
                   "CAT_FMT":"MEME",
                   "CAT_LINES":(1,100000),
                   "CAT_TYPE":"LONG",
                   "CAT_OUT":"",
                   "PARA_OUT":pkg_resources.resource_filename(__name__, "config/") + "/lephare_zphot_output.para",
                   "BD_SCALE":0,
                   "GLB_CONTEXT":-1,
                   "FORB_CONTEXT":-1,
                   "ERR_SCALE":(0.052,0.026,0.05,0.02,0.02,0.02,0.03),
                   "ERR_FACTOR":1.,
                   "ZPHOTLIB":("BC03","STAR","QSO"),
                   "ADD_EMLINES":"YES",
                   "FIR_LIB":"NONE",
                   "FIR_LMIN":7.,
                   "FIR_CONT":-1,
                   "FIR_SCALE":-1,
                   "FIR_FREESCALE":"NO",
                   "FIR_SUBSTELLAR":"NO",
                   "PHYS_LIB":"NONE",
                   "PHYS_CONT":-1,
                   "PHYS_SCALE":-1,
                   "PHYS_NMAX":100000,
                   "MASS_SCALE":(6.,16.),
                   "MAG_ABS":(-10.,-26.),
                   "MAG_REF":6,
                   "ZFORM_MIN":(5,5,5,5,5,5,3,1),
                   "Z_RANGE":(0,0.20),
                   "EBV_RANGE":(0.,1.),
                   "NZ_PRIOR":(4,2,4),
                   "ZFIX":"YES",
                   "Z_INTERP":"NO",
                   "DZ_WIN":0.25,
                   "MIN_THRES":0.1,
                   "PROB_INTZ":(0,0.5,0.5,1.,1.,1.5),
                   "MABS_METHOD":1,
                   "MABS_CONTEXT":-1,
                   "MABS_REF":4,
                   "MABS_FILT":(3,4,5,6),
                   "MABS_ZBIN":(0,0.5,1,1.5,2,3,3.5,4),
                   "SPEC_OUT":"YES",
                   "CHI2_OUT":"NO",
                   "PDZ_OUT":"NONE",
                   "PDZ_MABS_FILT":(2,10,14),
                   "FAST_MODE":"NO",
                   "COL_NUM":3,
                   "COL_SIGMA":3,
                   "COL_SEL":"AND",
                   "APPLY_SYSSHIFT":(0.007,-0.001,0.001,-0.003,0.006),
                   "AUTO_ADAPT":"NO",
                   "ADAPT_BAND":(4,2,4),
                   "ADAPT_LIM":(10,22.0),
                   "ADAPT_POLY":1,
                   "ADAPT_METH":1,
                   "ADAPT_CONTEXT":-1,
                   "ADAPT_ZBIN":(0.01,0.1),
                   "ADAPT_MODBIN":(1,1000),
                   "ERROR_ADAPT":"NO"}
    
    PATH_LEPHAREDIR = os.path.expanduser(os.getenv("LEPHAREDIR"))

    def __init__(self, **kwargs):
        """
        
        """
        self.set_data(**kwargs)

    def set_data(self, data=None, input_param_file=None, output_param_file=None, path_results=None, **kwargs):
        """
        
        """
        self._properties["input_param_file"] = pkg_resources.resource_filename(__name__, "config/") + "/lephare_zphot_input.para" \
                                               if input_param_file is None else input_param_file
        self._properties["input_param_file"] = pkg_resources.resource_filename(__name__, "config/") + "/lephare_zphot_output.para" \
                                               if input_param_file is None else input_param_file
        self._side_properties["path_results"] = pkg_resources.resource_filename(__name__, "results/") + "/" \
                                                if path_results is None else path_results
        
        list_commented_param = ["Z_FORM", "FORB_CONTEXT", "ZFORM_MIN", "NZ_PRIOR", "PROB_INTZ", "PDZ_MABS_FILT", "APPLY_SYSSHIFT"]
        for key, value in self.INPUT_PARAM.items():
            self.change_input_param(key, value, True is key in list_commented_param else False)

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
    
    def _set_data_path_(self, data_path):
        """
        
        """
        self._write_para_("CAT_IN", data_path)
    
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


