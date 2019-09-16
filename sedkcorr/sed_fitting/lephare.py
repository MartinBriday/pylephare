
import numpy as np
import pandas
import os
import subprocess
import pkg_resources
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)

from propobject import BaseObject
from astrobject import photometry
from astrobject.instruments import sdss

from ..k_correction import kcorrection
from ..k_correction.kcorrection import KCorrection


class LePhareSEDFitter( BaseObject ):
    """
    This class is a python "wrapper" (sort of) to run LePhare SED fitting.
    """

    PROPERTIES         = ["input_param_file", "output_param_file", "data_meas"]
    SIDE_PROPERTIES    = ["config_path", "results_path", "header_spec_file", "filt_list"]
    DERIVED_PROPERTIES = ["data_sed", "data_res"]
    
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

    def __init__(self, data=None, **kwargs):
        """
        The class constructor can automatically execute the 'set_data' method.
        
        Options
        -------
        data : [string or pandas.DataFrame or dict]
            Path of the data file or a DataFrame/dict, both of them in a format readable by LePhare fitter.
        
        input_param_file : [string or None]
            Path of the input parameter file.
            If 'None', the default file is imported from the package ('/config').
        
        output_param_file : [string or None]
            Path of the output parameter file.
            If 'None', the default file is imported from the package ('/config').
        
        results_path : [string or None]
            Path for the results of the SED fitter.
            If 'None', the default folder is located in the package ('/results').
        
        
        Returns
        -------
        Void
        """
        if data is not None:
            self.set_data(data, **kwargs)

    def set_data(self, data=None, input_param_file=None, output_param_file=None, results_path=None, flux_unit="Hz", data_filename="data.csv", **kwargs):
        """
        Set up the file paths about the data, the config files (input and output) and the results path.
        
        Parameters
        ----------
        data : [string or pandas.DataFrame or dict]
            Path of the data file or a DataFrame/dict, both of them in a format readable by LePhare fitter.
        
        input_param_file : [string or None]
            Path of the input parameter file.
            If 'None', the default file is imported from the package ('/config/lephare_zphot_input.para').
        
        output_param_file : [string or None]
            Path of the output parameter file.
            If 'None', the default file is imported from the package ('/config/lephare_zphot_output.para').
        
        results_path : [string or None]
            Path for the results of the SED fitter.
            If 'None', the default file is located in the package ('/results/data.out').
        
        flux_unit : [string]
            If 'data' is in flux, you can precise here the unit :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
            As LePhare needs the flux to be in the "Hz" (see above) unit, it will be converted if necessary.
        
        
        Returns
        -------
        Void
        """
        self._properties["input_param_file"] = self.config_path + "/lephare_zphot_input.para" \
                                               if input_param_file is None else os.path.abspath(input_param_file)
        self._properties["output_param_file"] = self.config_path + "/lephare_zphot_output.para" \
                                                if output_param_file is None else os.path.abspath(output_param_file)
        
        if data is not None:
            if type(data) == dict:
                try:
                    data = pandas.DataFrame(data)
                except(ValueError):
                    for k, v in data.items():
                        data[k] = [v]
                    data = pandas.DataFrame(data)
            if type(data) == pandas.DataFrame:
                data_path = data_filename if "/" in data_filename else self.config_path+"/"+data_filename
                data.to_csv(data_path, sep=" ", header=False)
            elif type(data) == str:
                data_path = data
            else:
                raise TypeError("data must be a DataFrame or a string")
            self.change_param("CAT_IN", os.path.abspath(data_path))
        self._set_input_type_()
        self._convert_flux_(flux_unit)

        self.change_param("PARA_OUT", self.output_param_file, False)
        self._set_results_path_(results_path)
            
    def _set_input_type_(self):
        """
        Automatically changes the input configuration file 'INP_TYPE' parameter either the input data or magnitudes or flux.
        
        
        Returns
        -------
        Void
        """
        if self._get_param_details_("INP_TYPE")[1] == "M" and self.data_meas.iloc[0]["mag_g"] < 1:
            self.change_param("INP_TYPE", "F", False)
        elif self._get_param_details_("INP_TYPE")[1] == "F" and self.data_meas.iloc[0]["flux_g"] > 1:
            self.change_param("INP_TYPE", "M", False)

    def _convert_flux_(self, flux_unit):
        """
        
        """
        data_buf = self.data_meas
        if self._get_param_details_("INP_TYPE")[1] == "F" and flux_unit != "Hz":
            for _filt in self.filt_list:
                data_buf["flux_"+_filt] = KCorrection.convert_flux_unit(data_buf["flux_"+_filt], kcorrection.FILTER_BANDS[_filt]["lbda"], flux_unit, "Hz")
                data_buf["flux_"+_filt+".err"] = KCorrection.convert_flux_unit(data_buf["flux_"+_filt+".err"], kcorrection.FILTER_BANDS[_filt]["lbda"], flux_unit, "Hz")
        data_buf.to_csv(self._get_param_details_("CAT_IN")[1], sep=" ", header=False)

    def _set_results_path_(self, results_path=None):
        """
        Set attribute 'results_path' and write the result file path in the config file.
        
        Parameters
        ----------
        results_path : [string or None]
            Path for the results of the SED fitter.
            If 'None', the default folder is located in the package ('/results').
        
        
        Returns
        -------
        Void
        """
        if results_path is None:
            results_path = self.results_path + "/data.out"
        
        self.change_param("CAT_OUT", os.path.abspath(results_path))
        self._side_properties["results_path"] = os.path.abspath("/".join(results_path.split("/")[:-1])) + "/"

    def _get_idx_line_(self, file, line):
        """
        Return the line index from a file, corresponding to a given parameter.
        
        Parameters
        ----------
        file : [list]
            List of the lines included in a file.
        
        line : [string or int]
            Name of the parameter we want to look for in the file.
            You can give an integer (the index of the line, the method will return the same number...
        
        
        Returns
        -------
        int
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
        Return the new parameter value in string in a "LePhare" format.
        
        Parameters
        ----------
        new_param_value : [string or int or float or tuple or list]
            New value we want to give to an input parameter.
        
        
        Returns
        -------
        str
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
        Return the splitted line and a comment flag (False by default).
        If it's recognized as a commented line, the '#' is removed, but the comment flag is set to True.
        
        Parameters
        ----------
        line : [string]
            Line of the file we want to split and clear.
        
        
        Returns
        -------
        list[str], bool
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
        Return a list of the lines included in a given file.
        
        Parameters
        ----------
        file : [string]
            Full path of the file.
        
        
        Returns
        -------
        list[string]
        """
        with open(file, "r") as file:
            file_buf = [line for line in file]
        return file_buf
    
    def _get_config_(self, param):
        """
        Return either 'input' or 'output', depending on the given parameter.
        
        Parameters
        ----------
        param : [string]
            Parameter we want to know either it is in the input parameter file or the output one.
        
        
        Returns
        -------
        str
        """
        if param in self.INPUT_PARAM:
            config = "input"
        elif param in self.OUTPUT_PARAM + ["PHYS_PARA{}_{}".format(ii, jj) for ii in self.PHYS_PARAM.keys() for jj in ["BEST", "MED", "INF", "SUP"]]:
            config = "output"
        else:
            raise ValueError("'{}' neither is in input parameter list nor output parameter list.".format(param))
        return config
    
    def change_param(self, param, new_param_value, force_comment=False):
        """
        This method change a parameter value in the input/output parameter file.
        You also can impose the parameter as a comment or not.
        
        Parameters
        ----------
        param : [string]
            Parameter we want to change the value and/or the comment state.
        
        new_param_value : [string or int or float or tuple or list]
            New value we want to give to an input parameter.
            Output parameters don't have values (only comment state), so the new value is not taken in acount for this case.
        
        force_comment : [bool]
            Impose a comment state.
            If True, the parameter is set as a comment (will write '#' at the beginnig of the corresponding line).
        
        
        Returns
        -------
        Void
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
        
        if param in ["{}_LIB".format(elt) for elt in ["STAR", "QSO", "GAL"]]:
            self.change_param(param+"_IN", new_param_value, False)
        elif param == "CAT_OUT":
            self._side_properties["results_path"] = os.path.abspath("/".join(new_param_value.split("/")[:-1]) + "/")

    def _get_param_details_(self, param):
        """
        Return the details (parameter name, value, comment state) of a given parameter
        
        Parameters
        ----------
        param : [string]
            Parameter we want to get the details.
            If the parameter is an output one, the value is a void string.
        
        
        Returns
        -------
        str, str, bool
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
        Print the details (parameter name, value, comment state) from a given parameter.
        
        Parameters
        ----------
        param : [string]
            Parameter we want to print the details.
        
        
        Returns
        -------
        Void
        """
        param, value, flag_comment = self._get_param_details_(param)
        
        txt_param = "{} : {} (comented? : {})".format(param, value, flag_comment)
        
        print(txt_param)

    def describe_params(self, which="input"):
        """
        This method prints the details (parameter name, value, comment state) of parameters.
        
        Parameters
        ----------
        which : [string]
            Either 'input' or 'output' or one or a list of parameter names inside configuration files.
            Either 'input' or 'output' will print every parameter details from the specific file.
            One or a list of parameters will restrict the print to the/these parameters.
        
        
        Returns
        -------
        Void
        """
        if which == "input":
            list_param = self.INPUT_PARAM
        elif which == "output":
            list_param = self.OUTPUT_PARAM
        elif type(which) in [str, list]:
            list_param = which if type(which) == list else [which]
        else:
            raise ValueError("'which_config' must be either 'input' or 'output' ; \n",
                             "'which_param' must be in [None, 'all', '*'], or one parameter or a list of parameters (in those two cases, which_config is not needed.")
        
        for _param in list_param:
            if "{}" in _param:
                for ii in self.PHYS_PARAM.keys():
                    self._print_param_details_(_param.format(ii))
            else:
                self._print_param_details_(_param)

    def _init_changes_(self, input_param_file=None, results_path=None, change_params=None):
        """
        LePhare shell commands initialization change options.
        
        Parameters
        ----------
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        results_path : [string or None]
            If you want to set a new results path, give the new path here.
            Default is 'None', which is the path set during the class construction or an execution of 'set_data'.
        
        change_params : [dict or None]
        If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, and as values a list containing the new parameter value as first element and the force comment option as second.
        
        
        Returns
        -------
        Void
        """
        if input_param_file is not None:
            self._properties["input_param_file"] = os.path.abspath(input_param_file)
        if results_path is not None:
            self._set_results_path_(results_path)
        if change_params is not None:
            try:
                for k, v in change_params.items():
                    self.change_param(param=k, new_param_value=v[0], force_comment=v[1])
            except(AttributeError):
                raise TypeError("'change_params' must be a dictionary containing the parameters you want to change as keys and for values a list with the new parameter value first and the force comment option second.")

    def run_filter(self, input_param_file=None, update=False, change_params=None):
        """
        Execute "$LEPHAREDIR/source/filter -c [...].para" in the shell.
        Exception : if the "$LEPAHAREWORK/filt/[...].filt" files already exist, the command is not executed, unless 'update' is True.
        
        Options
        -------
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        update : [bool]
            Set to True if you want to execute the command, even if the "$LEPAHAREWORK/filt/[...].filt" files already exist.
        
        change_params : [dict or None]
            If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, and as values a list containing the new parameter value as first element and the force comment option as second.
        
        
        Returns
        -------
        Void
        """
        self._init_changes_(input_param_file=input_param_file, results_path=None, change_params=change_params)
                
        filt_k, filt_v, _ = self._get_param_details_("FILTER_FILE")
        
        if not os.path.isfile(self.PATH_LEPHAREWORK+"/filt/"+filt_v) or update:
            cmd = "{}/source/filter -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
            subprocess.run(cmd.split())
    
    def run_sedtolib(self, input_param_file=None, update=False, change_params=None):
        """
        Execute "$LEPHAREDIR/source/sedtolib -t [S,Q,G] -c [...].para" in the shell.
        Exception : if the "$LEPAHAREWORK/lib_bin/[...].bin" files already exist, the command is not executed, unless 'update' is True.
        
        Options
        -------
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        update : [bool]
            Set to True if you want to execute the command, even if the "$LEPAHAREWORK/lib_bin/[...].bin" files already exist.
            
        change_params : [dict or None]
            If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, and as values a list containing the new parameter value as first element and the force comment option as second.
        
        
        Returns
        -------
        Void
        """
        self._init_changes_(input_param_file=input_param_file, results_path=None, change_params=change_params)
    
        lib_s_k, lib_s_v, _ = self._get_param_details_("STAR_LIB")
        lib_q_k, lib_q_v, _ = self._get_param_details_("QSO_LIB")
        lib_g_k, lib_g_v, _ = self._get_param_details_("GAL_LIB")
        
        if not np.prod([os.path.isfile(self.PATH_LEPHAREWORK+"/lib_bin/"+elt+".bin") for elt in [lib_s_v, lib_q_v, lib_g_v]]) or update:
            for elt in ["S", "Q", "G"]:
                cmd = "{}/source/sedtolib -t {} -c {}".format(self.PATH_LEPHAREDIR, elt, self.input_param_file)
                subprocess.run(cmd.split())
    
    def run_mag_star(self, input_param_file=None, update=False, change_params=None):
        """
        Execute "$LEPHAREDIR/source/mag_star -c [...].para" in the shell.
        Exception : if the "$LEPAHAREWORK/lib_mag/[...].bin" file already exists, the command is not executed, unless 'update' is True.
        
        Options
        -------
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        update : [bool]
            Set to True if you want to execute the command, even if the "$LEPAHAREWORK/lib_mag/[...].bin" file already exists.
            
        change_params : [dict or None]
            If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, and as values a list containing the new parameter value as first element and the force comment option as second.
            
            
        Returns
        -------
        Void
        """
        self._init_changes_(input_param_file=input_param_file, results_path=None, change_params=change_params)
        
        lib_s_k, lib_s_v, _ = self._get_param_details_("STAR_LIB_OUT")
        
        if not os.path.isfile(self.PATH_LEPHAREWORK+"/lib_mag/"+lib_s_v+".bin") or update:
            cmd = "{}/source/mag_star -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
            subprocess.run(cmd.split())
    
    def run_mag_gal(self, input_param_file=None, update=False, change_params=None):
        """
        Execute "$LEPHAREDIR/source/mag_gal -t [Q,G] -c [...].para" in the shell.
        Exception : if the "$LEPAHAREWORK/lib_mag/[...].bin" files already exist, the command is not executed, unless 'update' is True.
        
        Options
        -------
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        update : [bool]
            Set to True if you want to execute the command, even if the "$LEPAHAREWORK/lib_mag/[...].bin" files already exist.
            
        change_params : [dict or None]
            If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, and as values a list containing the new parameter value as first element and the force comment option as second.
            
            
        Returns
        -------
        Void
        """
        self._init_changes_(input_param_file=input_param_file, results_path=None, change_params=change_params)
        
        lib_q_k, lib_q_v, _ = self._get_param_details_("QSO_LIB_OUT")
        lib_g_k, lib_g_v, _ = self._get_param_details_("GAL_LIB_OUT")
        
        if not np.prod([os.path.isfile(self.PATH_LEPHAREWORK+"/lib_mag/"+elt+".bin") for elt in [lib_q_v, lib_g_v]]) or update:
            for elt in ["Q", "G"]:
                cmd = "{}/source/mag_gal -t {} -c {}".format(self.PATH_LEPHAREDIR, elt, self.input_param_file)
                subprocess.run(cmd.split())
    
    def run_zphota(self, input_param_file=None, results_path=None, change_params=None, change_context=None, savefile=None):
        """
        First change current directory to the results path.
        Then execute "$LEPHAREDIR/source/zphota -c [...].para" in the shell.
        
        Options
        -------
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        results_path : [string or None]
            If you want to set a new results path, give the new path here.
            Default is 'None', which is the path set during the class construction or an execution of 'set_data'.
            
        change_params : [dict or None]
            If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, and as values a list containing the new parameter value as first element and the force comment option as second.
        
        change_context : [list(string) or dict or None]
            If a list of filters is given, every context is changed to the corresponding one.
            If you want to change only a few indexes context, you can give a dictionary containing a list of the indexes under the key "id" and the list of filters under the key "filters".
            If None, nothing change from the 'data_meas'.
        
        savefile : [string or None]
            If not None, the 'data_sed' will be saved in the given file path.
            
            
        Returns
        -------
        Void
        """
        if change_context is not None:
            if type(change_context) == list:
                self.set_context_filters(id=None, filters=change_context)
            elif type(change_context) == dict:
                self.set_context_filters(**change_context)
        self._init_changes_(input_param_file=input_param_file, results_path=results_path, change_params=change_params)
        
        os.chdir(self.results_path)
        cmd = "{}/source/zphota -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
        subprocess.run(cmd.split())
        self.set_data_sed()
        self.set_data_res()
        if savefile is not None:
            self.write(savefile, None)

    def run_fit(self, input_param_file=None, results_path=None, update=False, change_params=None, change_context=None, savefile=None):
        """
        Run shell commands to execute LePhare fitting.
        
        Options
        -------
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        results_path : [string or None]
            If you want to set a new results path, give the new path here.
            Default is 'None', which is the path set during the class construction or an execution of 'set_data'.
            
        update : [bool]
            Set to True if you want to update the initialization on filters, sed libraries, etc.
            
        change_params : [dict or None]
            If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, and as values a list containing the new parameter value as first element and the force comment option as second.
        
        change_context : [list(string) or dict or None]
            If a list of filters is given, every context is changed to the corresponding one.
            If you want to change only a few indexes context, you can give a dictionary containing a list of the indexes under the key "id" and the list of filters under the key "filters".
            If None, nothing change from the 'data_meas'.
        
        savefile : [string or None]
            If not None, the 'data_sed' will be saved in the given file path.
        
        
        Returns
        -------
        Void
        """
        self._init_changes_(input_param_file=input_param_file, results_path=results_path, change_params=change_params)
        
        self.run_filter(input_param_file=input_param_file, update=update, change_params=change_params)
        self.run_sedtolib(input_param_file=input_param_file, update=update, change_params=change_params)
        self.run_mag_star(input_param_file=input_param_file, update=update, change_params=change_params)
        self.run_mag_gal(input_param_file=input_param_file, update=update, change_params=change_params)
        self.run_zphota(input_param_file=input_param_file, results_path=results_path, change_params=change_params, change_context=change_context, savefile=savefile)
    
    def _get_filt_list_(self):
        """
        Return the number of filters included in the configuration file ("FILTER_LIST").
        
        
        Returns
        -------
        int
        """
        _, filt_list, _ = self._get_param_details_("FILTER_LIST")
        return [self.lephare_filt_to_filt(kcorrection.FILTER_BANDS, _filt) for _filt in filt_list.split(",")]
    
    def _get_context_filters_(self, data):
        """
        Return a list of the concerned filter bands relative to the given context.
        Knowing the bands in 'self.filt_list', the context number is : sum(2**[band_nb]).
        For example, bands = ["u", "g", "r", "i", "z"] (band_nb = [0, 1, 2, 3, 4]) :
        - context = 31 --> ["u", "g", "r", "i", "z"]
        - context = 30 --> ["g", "r", "i", "z"]
        - context = 15 --> ["u", "g", "r", "i"]
        - context = 25 --> ["u", "i", "z"]
        - etc.
        
        
        Returns
        -------
        list(string)
        """
        context = int(data["CONTEXT"])
        idx_list = []
        for ii in np.arange(len(self.filt_list)-1,-1,-1):
            if (context - 2**ii) >= 0:
                context = context - 2**ii
                idx_list.append(ii)
        return [band for band in self.filt_list if self.filt_list.index(band) in idx_list]

    def set_context_filters(self, id=None, filters=["u", "g", "r", "i", "z"]):
        """
        Automatically changes the 'CONTEXT' value in 'data_meas' given the list filters to fit on.
        
        Parameters
        ----------
        id : [int or list(int) or None]
            Index(es) of the line(s) you want to change the context by the given filters.
        
        filters : [list(string)]
            List of filters to fit on.
        
        
        Returns
        -------
        Void
        """
        data_buf = self.data_meas
        idx_cross = [self.filt_list.index(_filt) for _filt in filters]
        context = np.sum([2**ii for ii in idx_cross])
        if id is None:
            data_buf["CONTEXT"] = np.array([context]*len(self.data_meas))
        else:
            for _id in np.atleast_1d(id):
                data_buf.at[_id, "CONTEXT"] = context
        data_buf.to_csv(self._get_param_details_("CAT_IN")[1], sep=" ", header=False)

    def _get_data_meas_(self):
        """
        Return a DataFrame of the input data.
        
        
        Returns
        -------
        pandas.DataFrame
        """
        data = pandas.read_csv(self._get_param_details_("CAT_IN")[1], sep=" ")
        prefix = "mag" if self._get_param_details_("INP_TYPE")[1] == "M" else "flux"
        names = []
        for _filt in self.filt_list:
            names.append("{}_{}".format(prefix, _filt))
            names.append("{}_{}.err".format(prefix, _filt))
        names += ["CONTEXT", "Z-SPEC"]
        names += ["STRING_{}".format(ii) for ii in np.arange(data.shape[1]-(2*len(self.filt_list)+3))]
        data = pandas.read_csv(self._get_param_details_("CAT_IN")[1], sep=" ", names=names)
        return data
    
    def _get_header_spec_file_(self):
        """
        
        """
        idx_start = 14
        sed_filename = self._get_sed_filename_(0)
        data_sed =  pandas.read_csv(os.path.expanduser(sed_filename),
                                    skiprows=idx_start, sep="  ", engine="python", nrows=1)
        while data_sed.shape[1] != 2:
            idx_start += 1
            data_sed =  pandas.read_csv(os.path.expanduser(sed_filename),
                                        skiprows=idx_start, sep="  ", engine="python", nrows=1)
        return idx_start
    
    def _get_data_sed_(self, sed_filename):
        """
        Return a DataFrame of the fitted SED (wavelength, magnitude).
        
        
        Returns
        -------
        pandas.DataFrame
        """
        data_sed =  pandas.read_csv(os.path.expanduser(sed_filename),
                                    skiprows=self.header_spec_file, names=["lbda", "mag"], sep="  ",
                                    engine="python", nrows=1050)
        return data_sed
                        
    def _get_sed_filename_(self, id_sed):
        """
        Return the SED spectrum file name (full path) given an index.
        
        Parameters
        ----------
        id_sed : [int]
            Index of the SED spectrum you want to get the file name.
        
        
        Returns
        -------
        string
        """
        id_sed = str(id_sed)
        while len(id_sed)<9:
            id_sed = "0"+id_sed
        sed_filename = "/Id"+id_sed+".spec"
        return self.results_path+sed_filename

    def set_data_sed(self):
        """
        Set the LePhare fit results in a dictionary containing the fitted spectrum tables.
        
        
        Returns
        -------
        Void
        """
        self._derived_properties["data_sed"] = {ii:self._get_data_sed_(self._get_sed_filename_(ii)) for ii in np.arange(len(self.data_meas))}
    
    def set_data_res(self):
        """
        Set the LePhare outfile as a pandas table attribute.
        
        
        Returns
        -------
        Void
        """
        self._derived_properties["data_res"] = self.lephare_output_file_reader(filename=self._get_param_details_("CAT_OUT")[1], filter_list=self.filt_list)

    def show(self, ax=None, id_sed=0, y_unit="AA", plot_phot=True, xlim=(None, None), ylim=(None, None), xscale="linear", yscale="linear", savefile=None, **kwargs):
        """
        Plot method.
        Return dict("fig", "ax").
        
        Parameters
        ----------
        id_sed : [int]
            Index of the SED you want to plot, corresponding to index in the 'data_meas' table.
        
        y_unit : [string]
            Choice to plot "mag" or flux with :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        Options
        -------
        ax : [matplotlib.axes]
            Already existing axes you want to add stuff in.
            Else, None.
        
        plot_phot : [bool]
            If True, plot the photometry points with errors, either in flux or magnitude.
        
        xlim : [tuple[float or None]]
            Set the limits on the x axis.
            If (None, None), the figure has free x axis limits.
        
        ylim : [tuple[float or None]]
            Set the limits on the y axis.
            If (None, None), the figure has free y axis limits.
        
        xscale : [string]
            Scale of the x axis : "linear", "log", ...
        
        yscale : [string]
            Scale of the y axis : "linear", "log", ...
        
        savefile : [string or None]
            If None, the figure won't be saved.
            To save it, input a path directory + filename.
        
        **kwargs : [dict]
            pyplot.plot options to apply on the SED spectrum.
        
        
        Returns
        -------
        dict
        """
        if ax is not None:
            fig = ax.figure
        else:
            fig, ax = plt.subplots()
        
        if len(self.data_meas) == 1:
            id_sed = 0
        elif id_sed is None:
            id_sed = -1
        x_sed = self.data_sed[id_sed]["lbda"]
        y_sed = self.data_sed[id_sed]["mag"]
        if y_unit in ["Hz", "AA", "mgy"]:
            y_sed, _ = KCorrection.mag_to_flux(y_sed, np.zeros(len(y_sed)), band=x_sed, flux_unit=y_unit, opt_mAB0=False)
        
        #SED
        opt_sed = {"ls":"-", "marker":"", "color":"0.4"}
        ax.plot(x_sed, y_sed, label="_nolegend_", **{**opt_sed, **kwargs})
        
        #Photometry
        if plot_phot:
            prefix = "mag" if self._get_param_details_("INP_TYPE")[1] == "M" else "flux"
            data_meas = self.data_meas if id_sed is not None else self.data_orig
            id_sed = 0 if id_sed is None else id_sed
            for _filt in self._get_context_filters_(data_meas.iloc[id_sed]):
                x_phot = kcorrection.FILTER_BANDS[_filt]["lbda"]
                y_phot = float(data_meas.iloc[id_sed]["{}_{}".format(prefix, _filt)])
                y_phot_err = float(data_meas.iloc[id_sed]["{}_{}.err".format(prefix, _filt)])
                if y_unit in ["Hz", "AA", "mgy"] and prefix == "mag":
                    y_phot, y_phot_err = KCorrection.mag_to_flux(y_phot, y_phot_err, band=_filt, flux_unit=y_unit, opt_mAB0=True)
                elif y_unit == "mag" and prefix == "flux":
                    y_phot, y_phot_err = KCorrection.flux_to_mag(y_phot, y_phot_err, band=_filt, flux_unit="Hz", opt_mAB0=True)
                elif y_unit in ["Hz", "AA", "mgy"] and prefix == "flux":
                    y_phot, y_phot_err = KCorrection.convert_flux_unit((y_phot, y_phot_err), kcorrection.FILTER_BANDS[_filt]["lbda"], "Hz", y_unit)
                ax.errorbar(x_phot, y_phot, yerr=y_phot_err, ls="", marker="o", color=kcorrection.FILTER_BANDS[_filt]["color"], label=_filt)
        
        #Writings
        ax.set_xlabel(r"$\lambda$ [\AA]", fontsize="large")
        ylabel = "mgy" if y_unit == "mgy" else \
                 r"${{f}}_{{\nu}}$ $[erg.{{s}}^{{-1}}.{{cm}}^{{-2}}.{Hz}^{{-1}}]$" if y_unit == "Hz" else \
                 r"${{f}}_{{\lambda}}$ $[erg.{{s}}^{{-1}}.{{cm}}^{{-2}}.{\AA}^{{-1}}]$" if y_unit == "AA" else \
                 "mag"
        ax.set_ylabel(ylabel, fontsize="large")
        ax.legend(loc="upper right", ncol=1)

        #Fig view
        ax.set_xlim(xlim)
        if ylim == (None, None):
            xmin, xmax = ax.get_xlim()
            mask = (xmin < np.asarray(x_sed)) * (np.asarray(x_sed) < xmax)
            ymin, ymax = np.min(np.asarray(y_sed)[mask]), np.max(np.asarray(y_sed)[mask])
            ylim = (ymin if ymin>0 else 0, ymax)
        ax.set_ylim(ylim)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        
        #Save
        if savefile is not None:
            fig.savefig(savefile)
        
        return {"ax":ax, "fig":fig}

    def write(self, savefile=None, id_sed=0):
        """
        Write 'data_sed' into a file.
        
        Parameters
        ----------
        savefile : [string]
            File path to save the data.
        
        id_sed : [int or None]
            Index of the 'data_sed' dictionary you want to save in a file.
            If None, the whole dictionary will be saved, each column name containing the id of the income.
        
        
        Returns
        -------
        Void
        """
        if type(id_sed) == int:
            self.data_sed[id_sed].to_csv(savefile, sep=" ", index=False)
        elif id_sed is None:
            data_out = self.data_sed[0].copy()
            data_out.set_index("lbda", inplace=True)
            data_out.rename(columns={k:k+"_id_0" for k in data_out.keys()}, inplace=True)
            for k1, v1 in self.data_sed.items():
                if k1 == 0:
                    continue
                pd_buf = v1.set_index("lbda").copy()
                for k2, v2 in pd_buf.items():
                    data_out[k2+"_id_{}".format(k1)] = v2
            data_out.reset_index(inplace=True)
            data_out.to_csv(savefile, sep=" ", index=False)
        else:
            raise ValueError("'id_sed' must be either an integer index contained in the 'data_sed' dictionary, or None if you want to save everything.")




    
    
    @staticmethod
    def lephare_filt_to_filt(dict, value):
        """
        Return the filter name of the input lephare filter transmission file name.
        
        Parameters
        ----------
        dict : [dict]
            Specific dictionnary containing as keys every filter and as values another
            dictionnary referencing the lephare filter transmission file names.
        
        value : [string]
            Lephare filter transmission file name you want to know the filter name.
        
        
        Returns
        -------
        string
        """
        for key1, value1 in dict.items():
            for key2, value2 in value1.items():
                if value2 == value:
                    return key1
        raise ValueError("{} is an unknown lephare filter syntax.".format(value))
    
    @staticmethod
    def lephare_output_file_reader(filename=None, filter_list=None):
        """
        Read a LePhare output file, returning a pandas.DataFrame.

        Parameters
        ----------
        filename : [string or None]
            Path of the LePhare file.
            If None, read the outfile file made at the end of the fit of the corresponding object.

        filter_list : [list(string) or None]
            List of filters used in LePhare.
            Necessary to create the good number of names for columns.


        Returns
        -------
        pandas.DataFrame
        """
        with open(filename, "r") as f1:
            buf_file = f1.readlines()
        buf_delimiter = buf_file[0]

        #Skiprows
        skiprows = [ii for ii, line in enumerate(buf_file) if line==buf_delimiter]
        skiprows = (skiprows[1] - skiprows[0])

        #Column names
        ii_col_names = [ii for ii, line in enumerate(buf_file) if "Output format" in line][0]
        col_names = [name for line in buf_file[ii_col_names+1:skiprows] for name in line.split(" ")]
        for elt in ["", ",", "\n", "#"]+[str(ii) for ii in range(1000)]:
            while elt in col_names: col_names.remove(elt)
        for ii, name in enumerate(col_names):
            if "()" in name:
                col_names.remove(name)
                for jj, filt in enumerate(filter_list):
                    col_names.insert(ii+jj, name.replace("()", "_{}".format(filt)))
        return pandas.read_csv(filename, sep=" ", skipinitialspace=True, skiprows=skiprows+1, names=col_names)
    
    
    #-------------------#
    #   Properties      #
    #-------------------#
    @property
    def input_param_file(self):
        """ Path of the input parameter file. """
        return self._properties["input_param_file"]

    @property
    def output_param_file(self):
        """ Path of the output parameter file. """
        return self._properties["output_param_file"]
    
    @property
    def config_path(self):
        """ Path of the configuration files """
        if self._side_properties["config_path"] is None:
            self._side_properties["config_path"] = pkg_resources.resource_filename(__name__, "config/") + "/"
        return self._side_properties["config_path"]
    
    @property
    def results_path(self):
        """ Path of the result path. """
        if self._side_properties["results_path"] is None:
            self._side_properties["results_path"] = pkg_resources.resource_filename(__name__, "results/") + "/"
        return self._side_properties["results_path"]
    
    @property
    def header_spec_file(self):
        """ Number of lines to skip when reading the spec files. """
        if self._side_properties["header_spec_file"] is None:
            self._side_properties["header_spec_file"] = self._get_header_spec_file_()
        return self._side_properties["header_spec_file"]

    @property
    def filt_list(self):
        """ List of filters set in the input configuration file. """
#        if self._side_properties["filt_list"] is None:
#            self._get_filt_list_()
#        return self._side_properties["filt_list"]
        return self._get_filt_list_()

    @property
    def data_meas(self):
        """ Input data """
#        if self._properties["data_meas"] is None:
#            self._get_data_meas_()
#        return self._properties["data_meas"]
        return self._get_data_meas_()

    @property
    def data_sed(self):
        """ Fitted SED data """
        if self._derived_properties["data_sed"] is None:
            self._derived_properties["data_sed"] = {}
        return self._derived_properties["data_sed"]

    @property
    def data_res(self):
        """ DataFrame of the LePhare results containing the wanted output parameters. """
        return self._derived_properties["data_res"]


class LePhareRand( LePhareSEDFitter ):
    """
    This class is a child of 'LePhareSEDFitter' class. Its goal is to fit SED errors.
    """
    
    PROPERTIES         = ["data_orig", "data_rand"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []

    def __init__(self, data=None, **kwargs):
        """
        The class constructor can automatically execute the 'set_data' method.
        
        Options
        -------
        data : [string or pandas.DataFrame or dict]
            Path of the data file or a DataFrame/dict, both of them in a format readable by LePhare fitter.
            Only accept the data for one by one object.
        
        input_param_file : [string or None]
            Path of the input parameter file.
            If 'None', the default file is imported from the package ('/config/lephare_zphot_input.para').
        
        output_param_file : [string or None]
            Path of the output parameter file.
            If 'None', the default file is imported from the package ('/config/lephare_zphot_output.para').
        
        results_path : [string or None]
            Path for the results of the SED fitter.
            If 'None', the default file is located in the package ('/results/data.out').
        
        flux_unit : [string]
            If 'data' is in flux, you can precise here the unit :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
            As LePhare needs the flux to be in the "Hz" (see above) unit, it will be converted if necessary.
        
        nb_draw : [int]
            Number of draw for the Monte Carlo fitting errors on the SED.
        
        
        Returns
        -------
        Void
        """
        if data is not None:
            self.set_data(data, **kwargs)

    def set_data(self, data=None, input_param_file=None, output_param_file=None, results_path=None, flux_unit="Hz", data_filename="data_rand.csv", nb_draw=100, **kwargs):
        """
        Set up the file paths about the data, the config files (input and output) and the results path.
        
        Parameters
        ----------
        data : [string or pandas.DataFrame or dict]
            Path of the data file or a DataFrame/dict, both of them in a format readable by LePhare fitter.
            Only accept the data for one by one object.
        
        input_param_file : [string or None]
            Path of the input parameter file.
            If 'None', the default file is imported from the package ('/config/lephare_zphot_input.para').
        
        output_param_file : [string or None]
            Path of the output parameter file.
            If 'None', the default file is imported from the package ('/config/lephare_zphot_output.para').
        
        results_path : [string or None]
            Path for the results of the SED fitter.
            If 'None', the default file is located in the package ('/results/data.out').
        
        flux_unit : [string]
            If 'data' is in flux, you can precise here the unit :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
            As LePhare needs the flux to be in the "Hz" (see above) unit, it will be converted if necessary.
        
        nb_draw : [int]
            Number of draw for the Monte Carlo fitting errors on the SED.
        
        
        Returns
        -------
        Void
        """
        _ = super(LePhareRand, self).set_data(data=data, input_param_file=input_param_file, output_param_file=output_param_file,
                                              results_path=results_path, flux_unit=flux_unit, data_filename=data_filename, **kwargs)
        self.set_rand_data(nb_draw=nb_draw, data_filename=data_filename)
    
    def set_rand_data(self, nb_draw=100, data_filename="data_rand.csv"):
        """
        Initialize the Monte Carlo random data to fit through LePhare in order to get SED errors.
        
        Parameters
        ----------
        nb_draw : [int]
            Number of draw for the Monte Carlo fitting errors on the SED.
        
        
        Returns
        -------
        Void
        """
        self._properties["data_orig"] = self.data_meas
        for _filt in self.filt_list:
            photopoint = photometry.PhotoPoint(flux=self.data_meas["flux_{}".format(_filt)],
                                               var=self.data_meas["flux_{}.err".format(_filt)]**2,
                                               lbda=kcorrection.FILTER_BANDS[_filt]["lbda"])
            photopoint.draw_photosamplers(nsamplers=nb_draw, negative_fluxmag=40)
            self.data_rand["flux_{}".format(_filt)] = photopoint.photosamplers.samplers
            self.data_rand["flux_{}.err".format(_filt)] = np.array([self.data_meas["flux_{}.err".format(_filt)][0]]*nb_draw)
        self.data_rand["CONTEXT"] = np.array([self.data_meas["CONTEXT"][0]]*nb_draw)
        self.data_rand["Z-SPEC"] = np.array([self.data_meas["Z-SPEC"][0]]*nb_draw)
        for _col_name in self.data_meas.columns:
            if "STRING" in _col_name:
                self.data_rand[_col_name] = np.array([self.data_meas[_col_name][0]]*nb_draw)
        self._properties["data_rand"] = pandas.DataFrame(self.data_rand)

        data_path = data_filename if "/" in data_filename else self.config_path+"/"+data_filename
        self.data_rand.to_csv(data_path, sep=" ", header=False)
        self.change_param("CAT_IN", os.path.abspath(data_path))
        self._set_input_type_()
    
    def _get_fit_quantiles_(self, quants=[0.16, 0.5, 0.84], y_unit="AA"):
        """
        Return the given quantiles on the Monte Carlo fitted spectra.
        
        Parameters
        ----------
        quants : [list(float)]
            Quantiles you want to get on the spectra. Each quantile must be between 0 and 1.
        
        y_unit : [string]
            Unit of the output spectrum quantiles:
            - "mag" : magnitudes
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        
        Returns
        -------
        Void
        """
        lbda = self.data_sed[0]["lbda"]
        mags = np.quantile([v["mag"] for k, v in self.data_sed.items() if (k != -1 and len(v) != 0)], quants, axis=0)
        return {k:m if y_unit=="mag" else KCorrection.mag_to_flux(m, np.zeros(len(m)), band=lbda, flux_unit=y_unit, opt_mAB0=False)[0] for k, m in zip(quants, mags)}
    
    def set_data_sed(self):
        """
        Set the LePhare fit results in a dictionary.
        Add the Monte Carlo results with a dictionary containing the wavelength, the spectra median and one sigma error below and above the median.
        
        
        Returns
        -------
        Void
        """
        _ = super(LePhareRand, self).set_data_sed()
        quants = self._get_fit_quantiles_(quants=[0.16, 0.5, 0.84], y_unit="mag")
        pd_none = {"lbda":self.data_sed[0]["lbda"], "mag":quants[0.5], "mag.err_low":quants[0.5]-quants[0.16], "mag.err_up":quants[0.84]-quants[0.5]}
        self.data_sed[-1] = pandas.DataFrame(pd_none)

    def show(self, ax=None, id_sed=None, y_unit="AA", plot_phot=True, xlim=(None, None), ylim=(None, None), xscale="linear", yscale="linear", savefile=None, show_sigmas=[1, 2], **kwargs):
        """
        Plot method.
        Return dict("fig", "ax").
        
        Parameters
        ----------
        id_sed : [int or None]
            Index of the SED you want to plot, corresponding to index in the 'data_meas' table.
            If None, plot the Monte Carlo fitted median of the spectra, and sigmas if specified (see 'show_sigmas').
        
        y_unit : [string]
            Choice to plot "mag" or flux with :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        Options
        -------
        ax : [matplotlib.axes]
            Already existing axes you want to add stuff in.
            Else, None.
        
        plot_phot : [bool]
            If True, plot the photometry points with errors, either in flux or magnitude.
        
        xlim : [tuple[float or None]]
            Set the limits on the x axis.
            If (None, None), the figure has free x axis limits.
        
        ylim : [tuple[float or None]]
            Set the limits on the y axis.
            If (None, None), the figure has free y axis limits.
        
        xscale : [string]
            Scale of the x axis : "linear", "log", ...
        
        yscale : [string]
            Scale of the y axis : "linear", "log", ...
        
        savefile : [string or None]
            If None, the figure won't be saved.
            To save it, input a path directory + filename.
        
        show_sigmas : [int or list(int) or None]
            Show 1 and/or 2 (or None) sigmas error of the SED.
        
        **kwargs : [dict]
            pyplot.plot options to apply on the SED spectrum.
        
        
        Returns
        -------
        dict
        """
        kwargs = {**{"zorder":3}, **kwargs}
        dict_fig = super(LePhareRand, self).show(ax=ax, id_sed=id_sed, y_unit=y_unit, plot_phot=plot_phot, xlim=xlim, ylim=ylim, xscale=xscale, yscale=yscale, savefile=savefile, **kwargs)
        
        if id_sed is None:
            lbda = self.data_sed[0]["lbda"]
            nsigmas = len(np.atleast_1d(show_sigmas))
            if 2 in show_sigmas:
                ff = self._get_fit_quantiles_(quants=[0.05, 0.95], y_unit=y_unit)
                dict_fig["ax"].fill_between(lbda, ff[0.05], ff[0.95], alpha=0.3/nsigmas, color="C0", lw=0, zorder=1)
            if 1 in show_sigmas:
                ff = self._get_fit_quantiles_(quants=[0.16, 0.84], y_unit=y_unit)
                dict_fig["ax"].fill_between(lbda, ff[0.16], ff[0.84], alpha=0.3/nsigmas, color="C0", lw=0, zorder=2)
        
        return dict_fig



    #-------------------#
    #   Properties      #
    #-------------------#
    @property
    def data_orig(self):
        """ Origin measurement data """
        return self._properties["data_orig"]
    
    @property
    def data_rand(self):
        """ Table of the random draws from the given data. """
        if self._properties["data_rand"] is None:
            self._properties["data_rand"] = {}
        return self._properties["data_rand"]

