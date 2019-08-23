
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
                   "AGE_RANGE":(0.,13.55e9),
                   "FILTER_LIST":["galex/FUV.pb","galex/NUV.pb","sdss/up.pb","sdss/gp.pb","sdss/rp.pb","sdss/ip.pb","sdss/zp.pb"],
                   "TRANS_TYPE":0,
                   "FILTER_CALIB":0,
                   "FILTER_FILE":"snf_galex_sdss.filt",
                   "STAR_LIB_IN":"snf_LIB_STAR",
                   "STAR_LIB_OUT":"snf_STAR",
                   "QSO_LIB_IN":"snf_LIB_QSO",
                   "QSO_LIB_OUT":"snf_QSO",
                   "GAL_LIB_IN":"snf_LIB_BC03",
                   "GAL_LIB_OUT":"snf_BC03",
                   "MAGTYPE":"AB",
                   "Z_STEP":(0.004,0.20,0.1),
                   "COSMOLOGY":(70,0.3,0.7),
                   "MOD_EXTINC":(1,27),
                   "EXTINC_LAW":"LMC_Fitzpatrick.dat",
                   "EB_V":(0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8),
                   "EM_LINES":"YES",
                   #"Z_FORM":(8,7,6,5,4,3),
                   "LIB_ASCII":"YES",
                   "CAT_IN":"$LEPHAREDIR/data/SNf/SNf_mag_1kpc_124_lephare.csv",
                   "INP_TYPE":"M",
                   "CAT_MAG":"AB",
                   "CAT_FMT":"MEME",
                   "CAT_LINES":(1,100000),
                   "CAT_TYPE":"LONG",
                   "CAT_OUT":"$LEPHAREDIR/data/SNF/SNf_mag_1kpc_124_lephare.out",
                   "PARA_OUT":"$LEPHAREDIR/config/snf_host_zphot_output.para",
                   "BD_SCALE":0,
                   "GLB_CONTEXT":-1,
                   #"FORB_CONTEXT":-1,
                   "ERR_SCALE":(0.052,0.026,0.05,0.02,0.02,0.02,0.03),
                   "ERR_FACTOR":1.,
                   "ZPHOTLIB":("snf_BC03","snf_STAR","snf_QSO"),
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
                   #"ZFORM_MIN":(5,5,5,5,5,5,3,1),
                   "Z_RANGE":(0,0.20),
                   "EBV_RANGE":(0.,1.),
                   #"NZ_PRIOR":(4,2,4),
                   "ZFIX":"YES",
                   "Z_INTERP":"NO",
                   "DZ_WIN":0.25,
                   "MIN_THRES":0.1,
                   #"PROB_INTZ":(0,0.5,0.5,1.,1.,1.5),
                   "MABS_METHOD":1,
                   "MABS_CONTEXT":-1,
                   "MABS_REF":4,
                   "MABS_FILT":(3,4,5,6),
                   "MABS_ZBIN":(0,0.5,1,1.5,2,3,3.5,4),
                   "SPEC_OUT":"YES",
                   "CHI2_OUT":"NO",
                   "PDZ_OUT":"NONE",
                   #"PDZ_MABS_FILT":(2,10,14),
                   "FAST_MODE":"NO",
                   "COL_NUM":3,
                   "COL_SIGMA":3,
                   "COL_SEL":"AND",
                   #"APPLY_SYSSHIFT":(0.007,-0.001,0.001,-0.003,0.006),
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

    def _get_idx_line_(self, line):
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
    
    def _write_param_(self, line, new_param_value):
        """
        
        """
        idx_line = self._get_idx_line_(line)
        new_param_value = self._get_new_param_value_(new_param_value)
        
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


