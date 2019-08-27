
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
    SIDE_PROPERTIES    = ["results_path"]
    DERIVED_PROPERTIES = []
    
    INPUT_PARAM = ["STAR_SED", "STAR_FSCALE", "STAR_LIB", "QSO_SED", "QSO_FSCALE", "QSO_LIB", "GAL_SED", "GAL_FSCALE", "GAL_LIB",
                   "SEL_AGE", "AGE_RANGE", "FILTER_LIST", "TRANS_TYPE", "FILTER_CALIB", "FILTER_FILE", "STAR_LIB_IN", "STAR_LIB_OUT",
                   "QSO_LIB_IN", "QSO_LIB_OUT", "GAL_LIB_IN", "GAL_LIB_OUT", "MAGTYPE", "Z_STEP", "COSMOLOGY", "MOD_EXTINC",
                   "EXTINC_LAW", "EB_V", "EM_LINES", "Z_FORM", "LIB_ASCII", "CAT_IN", "INP_TYPE", "CAT_MAG", "CAT_FMT", "CAT_LINES",
                   "CAT_TYPE", "CAT_OUT", "PARA_OUT", "BD_SCALE", "GLB_CONTEXT", "FORB_CONTEXT", "ERR_SCALE", "ERR_FACTOR", "ZPHOTLIB",
                   "ADD_EMLINES", "FIR_LIB", "FIR_LMIN", "FIR_CONT", "FIR_SCALE", "FIR_FREESCALE", "FIR_SUBSTELLAR", "PHYS_LIB",
                   "PHYS_CONT", "PHYS_SCALE", "PHYS_NMAX", "MASS_SCALE", "MAG_ABS", "MAG_REF", "ZFORM_MIN", "Z_RANGE", "EBV_RANGE",
                   "NZ_PRIOR", "ZFIX", "Z_INTERP", "DZ_WIN", "MIN_THRES", "PROB_INTZ", "MABS_METHOD", "MABS_CONTEXT", "MABS_REF",
                   "MABS_FILT", "MABS_ZBIN", "SPEC_OUT", "CHI2_OUT", "PDZ_OUT", "PDZ_MABS_FILT", "FAST_MODE", "COL_NUM", "COL_SIGMA",
                   "COL_SEL", "APPLY_SYSSHIFT", "AUTO_ADAPT", "ADAPT_BAND", "ADAPT_LIM", "ADAPT_POLY", "ADAPT_METH", "ADAPT_CONTEXT",
                   "ADAPT_ZBIN", "ADAPT_MODBIN", "ERROR_ADAPT"]
    
    OUTPUT_PARAM = ["IDENT", "Z_BEST", "Z_BEST68_LOW", "Z_BEST68_HIGH", "Z_ML", "Z_ML68_LOW", "Z_ML68_HIGH", "Z_BEST90_LOW",
                   "Z_BEST90_HIGH", "Z_BEST99_LOW", "Z_BEST99_HIGH", "CHI_BEST", "MOD_BEST", "EXTLAW_BEST", "EBV_BEST", "ZF_BEST",
                   "MAG_ABS_BEST", "PDZ_BEST", "SCALE_BEST", "DIST_MOD_BEST", "NBAND_USED", "NBAND_ULIM", "Z_SEC", "CHI_SEC",
                   "MOD_SEC", "AGE_SEC", "EBV_SEC", "ZF_SEC", "MAG_ABS_SEC", "PDZ_SEC", "SCALE_SEC", "Z_QSO", "CHI_QSO", "MOD_QSO",
                   "MAG_ABS_QSO", "DIST_MOD_QSO", "MOD_STAR", "CHI_STAR", "MAG_OBS()", "ERR_MAG_OBS()", "MAG_MOD()", "K_COR()",
                   "MAG_ABS()", "MABS_FILT()", "K_COR_QSO()", "MAG_ABS_QSO()", "PDZ()", "CONTEXT", "ZSPEC", "STRING_INPUT",
                   "LUM_TIR_BEST", "LIB_FIR", "MOD_FIR", "CHI2_FIR", "FSCALE_FIR", "NBAND_FIR", "LUM_TIR_MED", "LUM_TIR_INF",
                   "LUM_TIR_SUP", "MAG_MOD_FIR()", "MAG_ABS_FIR()", "K_COR_FIR()", "AGE_BEST", "AGE_INF", "AGE_MED", "AGE_SUP",
                   "LDUST_BEST", "LDUST_INF", "LDUST_MED", "LDUST_SUP", "MASS_BEST", "MASS_INF", "MASS_MED", "MASS_SUP", "SFR_BEST",
                   "SFR_INF", "SFR_MED", "SFR_SUP", "SSFR_BEST", "SSFR_INF", "SSFR_MED", "SSFR_SUP", "LUM_NUV_BEST", "LUM_R_BEST",
                   "LUM_K_BEST","PHYS_CHI2_BEST", "PHYS_MOD_BEST", "PHYS_MAG_MOD()", "PHYS_MAG_ABS()", "PHYS_K_COR()",
                   "PHYS_PARA{}_BEST", "PHYS_PARA{}_MED", "PHYS_PARA{}_INF", "PHYS_PARA{}_SUP"]
                   
    PHYS_PARAM = {"1":"Lg(Age[yr])",
                  "2":"Lg(Luv[Lo])",
                  "3":"Lg(Lr[Lo])",
                  "4":"Lg(Lk[Lo])",
                  "5":"Lg(Ldust [Lo])",
                  "6":"Lg(Mass[Mo])",
                  "7": "Lg(SFR_1e8[Mo/yr])",
                  "8":"Lg(SFR_1e9[Mo/yr])",
                  "9":"Lg(Gamma[Gyr])",
                  "10":"Metallicity",
                  "11":"Tauv",
                  "12":"mu",
                  "13":"Lg(Age_weightRlight[yr])",
                  "14":"Lg(Age_weightMass[yr])",
                  "15":"Lg(Age_tlastburst[yr])",
                  "16":"Log(Frac_burst_1.e8[yr])",
                  "17":"Lg(Frac_burst_1e9[yr])",
                  "18":"Nb burst",
                  "19":"Lg(SSFR_1e8 [yr-1])",
                  "20":"Lg(SSFR_1e9 [yr-1])",
                  "21":"A(1500) [mag]",
                  "22":"A(2300) [mag]",
                  "23":"A(4500) [mag]",
                  "24":"A(6000) [mag]",
                  "25":"A(8500) [mag]",
                  "26":"A(22000)[mag]",
                  "27":"Lg(Ldust_hot[Lo])"}
    
    PATH_LEPHAREDIR = os.path.expanduser(os.getenv("LEPHAREDIR"))
    PATH_LEPHAREWORK = os.path.expanduser(os.getenv("LEPHAREWORK"))

    def __init__(self, data_path=None, **kwargs):
        """
        
        """
        self.set_data(data_path, **kwargs)

    def set_data(self, data_path=None, input_param_file=None, output_param_file=None, results_path=None, **kwargs):
        """
        
        """
        self._properties["input_param_file"] = pkg_resources.resource_filename(__name__, "config/") + "/lephare_zphot_input.para" \
                                               if input_param_file is None else os.path.abspath(input_param_file)
        self._properties["output_param_file"] = pkg_resources.resource_filename(__name__, "config/") + "/lephare_zphot_output.para" \
                                               if output_param_file is None else os.path.abspath(output_param_file)
        
        self.change_param("PARA_OUT", self.output_param_file, False)
        
        if data_path is not None:
            self.change_param("CAT_IN", os.path.abspath(data_path))
        if results_path is not None:
            self.change_param("CAT_OUT", os.path.abspath(results_path))
        
        if results_path is not None:
            results_path = "/".join(results_path.split("/")[:-1]) + "/"
        self._side_properties["results_path"] = pkg_resources.resource_filename(__name__, "results/") + "/" \
            if results_path is None else os.path.abspath(results_path)

    def _get_idx_line_(self, file, line):
        """
        
        """
        if type(line) == int:
            idx_line = line
        elif type(line) == str:
            flag_no_match = True
            for ii, _line in enumerate(file):
                splitted_line = _line.split()
                if len(splitted_line)>0 and (line in splitted_line[0] or (len(splitted_line)>1 and splitted_line[1] == line)):
                    idx_line = ii
                    flag_no_match = False
                    break
            if flag_no_match:
                raise ValueError("{} is not a known parameter in the parameters file.".format(line))
        else:
            raise ValueError("'line' ( = {}) must be either integer (line index) or string (first word of the line)".format(line))
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
    
    def _get_cleared_line_(self, line):
        """
        
        """
        flag_comment = False
        splitted_line = line.split()
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
    
    def _get_config_(self, param):
        """
        
        """
        if param in self.INPUT_PARAM:
            config = "input"
        elif param in self.OUTPUT_PARAM + ["PHYS_PARA{}_{}".format(ii, jj) for ii in self.PHYS_PARAM.keys() for jj in ["BEST", "MED", "INF", "SUP"]]:
            config = "output"
        else:
            raise ValueError("'{}' neither is in input parameter list nor output parameter list.".format(param))
        return config
    
    def change_param(self, param, new_param_value=None, force_comment=False):
        """
        
        """
        config = self._get_config_(param)
        file_buf = self._get_file_lines(self._properties[config+"_param_file"])
        idx_line = self._get_idx_line_(file_buf, param)
        new_param_value = self._get_new_param_value_(new_param_value)
        
        splitted_line, _ = self._get_cleared_line_(file_buf[idx_line])
        if config == "input":
            splitted_line[1] = new_param_value
        if force_comment:
            splitted_line.insert(0, "#")
        file_buf[idx_line] = " ".join(splitted_line) + "\n"

        with open(self._properties[config+"_param_file"], "w") as file:
            for line in file_buf:
                file.write(line)

    def _get_param_details_(self, param):
        """
        
        """
        config = self._get_config_(param)
        file_buf = self._get_file_lines(self._properties[config+"_param_file"])
        idx_line = self._get_idx_line_(file_buf, param)
        splitted_line, flag_comment = self._get_cleared_line_(file_buf[idx_line])
        
        if "PHYS_PARA" in param:
            phys_param_idx = ["PARA{}".format(ii) for ii in self.PHYS_PARAM.keys()].index(param.split("_")[1]) + 1
            param += " ({})".format(self.PHYS_PARAM[str(phys_param_idx)])
        
        return param, splitted_line[1] if config=="input" else "", flag_comment

    def _print_param_details_(self, param):
        """
        
        """
        param, value, flag_comment = self._get_param_details_(param)
        
        txt_param = "{} : {} (comented? : {})".format(param, value, flag_comment)
        
        print(txt_param)

    def describe_params(self, which_config="input", which_param=None):
        """
        
        """
        if which_config == "input" and which_param in [None, "all", "*"]:
            list_param = self.INPUT_PARAM
        elif which_config == "output" and which_param in [None, "all", "*"]:
            list_param = self.OUTPUT_PARAM
        elif which_param not in [None, "all", "*"] and type(which_param) in [str, list]:
            list_param = which_param if type(which_param) == list else [which_param]
        else:
            raise ValueError("'which_config' must be either 'input' or 'output' ; \n",
                             "'which_param' must be in [None, 'all', '*'], or one parameter or a list of parameters (in those two cases, which_config is not needed.")
        
        for _param in list_param:
            if "{}" in _param:
                for ii in self.PHYS_PARAM.keys():
                    self._print_param_details_(_param.format(ii))
            else:
                self._print_param_details_(_param)
    
    def run_sedtolib(self, input_param_file=None, update=False):
        """
        
        """
        if input_param_file is not None:
            self._properties["input_param_file"] = os.path.abspath(input_param_file)
    
        lib_s_k, lib_s_v, _ = self._get_param_details_("STAR_LIB")
        lib_q_k, lib_q_v, _ = self._get_param_details_("QSO_LIB")
        lib_g_k, lib_g_v, _ = self._get_param_details_("GAL_LIB")
        
        if not np.prod([os.path.isfile(self.PATH_LEPHAREWORK+"/lib_bin/"+elt+".bin") for elt in [lib_s_v, lib_q_v, lib_g_v]]) or update:
            for elt in ["S", "Q", "G"]:
                cmd = "{}/source/sedtolib -t {} -c {}".format(self.PATH_LEPHAREDIR, elt, self.input_param_file)
                subprocess.run(cmd.split())
    
    def run_mag_star(self, input_param_file=None, update=False):
        """
        
        """
        if input_param_file is not None:
            self._properties["input_param_file"] = os.path.abspath(input_param_file)
        
        lib_s_k, lib_s_v, _ = self._get_param_details_("STAR_LIB_OUT")
        
        if not os.path.isfile(self.PATH_LEPHAREWORK+"/lib_mag/"+lib_s_v+".bin") or update:
            cmd = "{}/source/mag_star -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
            subprocess.run(cmd.split())
    
    def run_mag_gal(self, input_param_file=None, update=False):
        """
        
        """
        if input_param_file is not None:
            self._properties["input_param_file"] = os.path.abspath(input_param_file)
        
        lib_q_k, lib_q_v, _ = self._get_param_details_("QSO_LIB_OUT")
        lib_g_k, lib_g_v, _ = self._get_param_details_("GAL_LIB_OUT")
        
        if not np.prod([os.path.isfile(self.PATH_LEPHAREWORK+"/lib_mag/"+elt+".bin") for elt in [lib_q_v, lib_g_v]]) or update:
            for elt in ["Q", "G"]:
                cmd = "{}/source/mag_gal -t {} -c {}".format(self.PATH_LEPHAREDIR, elt, self.input_param_file)
                subprocess.run(cmd.split())
    
    def run_zphota(self, input_param_file=None, results_path=None):
        """
        
        """
        if input_param_file is not None:
            self._properties["input_param_file"] = os.path.abspath(input_param_file)
        if results_path is not None:
            self._side_properties["results_path"] = os.path.abspath(results_path)
        
        os.chdir(self.results_path)
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
    def results_path(self):
        """  """
        return self._side_properties["results_path"]


