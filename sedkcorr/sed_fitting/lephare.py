
import numpy as np
import pandas
import os
import subprocess
import pkg_resources

from propobject import BaseObject


class LePhareSEDFitter( BaseObject ):
    """
    This class is a python "wrapper" (sort of) to run LePhare SED fitting.
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

    def set_data(self, data=None, input_param_file=None, output_param_file=None, results_path=None, **kwargs):
        """
        Set up the file paths about the data, the config files (input and output) and the results path.
        
        Parameters
        ----------
        data : [string or pandas.DataFrame or dict]
            Path of the data file or a DataFrame/dict, both of them in a format readable by LePhare fitter.
            Only accept the data for one by one object.
        
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
        self._properties["input_param_file"] = pkg_resources.resource_filename(__name__, "config/") + "/lephare_zphot_input.para" \
                                               if input_param_file is None else os.path.abspath(input_param_file)
        self._properties["output_param_file"] = pkg_resources.resource_filename(__name__, "config/") + "/lephare_zphot_output.para" \
                                               if output_param_file is None else os.path.abspath(output_param_file)
        
        if data is not None:
            if type(data) == dict:
                try:
                    data = pandas.DataFrame(data)
                except(ValueError):
                    for k, v in data.items():
                        data[k] = [v]
                    data = pandas.DataFrame(data)
            if type(data) is pandas.DataFrame:
                data_path = pkg_resources.resource_filename(__name__, "config/")+"/data.csv"
                data.to_csv(data_path, sep=" ", header=False)
            elif type(data) is str:
                data_path = data
            else:
                raise TypeError("data must be a DataFrame or a string")
            self.change_param("CAT_IN", os.path.abspath(data_path))

        self.change_param("PARA_OUT", self.output_param_file, False)

        if results_path is None:
            results_path = pkg_resources.resource_filename(__name__, "results/")+"/data.out"
        
        self.change_param("CAT_OUT", os.path.abspath(results_path))
        self._side_properties["results_path"] = os.path.abspath("/".join(results_path.split("/")[:-1]) + "/")

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

    def describe_params(self, which_config="input", which_param=None):
        """
        This method prints the details (parameter name, value, comment state) of parameters.
        
        Parameters
        ----------
        which_config : [string]
            Either 'input' or 'output'.
            Only useful if we want to print every parameter details.
        
        which_param : [string or list(string) or None]
            If None, every parameter details are printed.
            If one or list of parameters, each is automatically found in the input/output file (so that 'which_config' is useless).
        
        
        Returns
        -------
        Void
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

    def run_filter(self, input_param_file=None, update=False):
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
        
        
        Returns
        -------
        Void
        """
        if input_param_file is not None:
            self._properties["input_param_file"] = os.path.abspath(input_param_file)
                
        filt_k, filt_v, _ = self._get_param_details_("FILTER_FILE")
        
        if not os.path.isfile(self.PATH_LEPHAREWORK+"/filt/"+filt_v+".filt") or update:
            cmd = "{}/source/filter -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
            subprocess.run(cmd.split())
    
    def run_sedtolib(self, input_param_file=None, update=False):
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
        
        
        Returns
        -------
        Void
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
        Execute "$LEPHAREDIR/source/mag_star -c [...].para" in the shell.
        Exception : if the "$LEPAHAREWORK/lib_mag/[...].bin" file already exists, the command is not executed, unless 'update' is True.
        
        Options
        -------
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        update : [bool]
            Set to True if you want to execute the command, even if the "$LEPAHAREWORK/lib_mag/[...].bin" file already exists.
            
            
        Returns
        -------
        Void
        """
        if input_param_file is not None:
            self._properties["input_param_file"] = os.path.abspath(input_param_file)
        
        lib_s_k, lib_s_v, _ = self._get_param_details_("STAR_LIB_OUT")
        
        if not os.path.isfile(self.PATH_LEPHAREWORK+"/lib_mag/"+lib_s_v+".bin") or update:
            cmd = "{}/source/mag_star -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
            subprocess.run(cmd.split())
    
    def run_mag_gal(self, input_param_file=None, update=False):
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
            
            
        Returns
        -------
        Void
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
        First change current directory to the results path.
        Then execute "$LEPHAREDIR/source/zphota -c [...].para" in the shell.
        
        Options
        -------
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        results_path : [string or None]  !!!!!!!!!!!!!!!!!!!!!
            If you want to set a new results path, give the new path here.
            Default is 'None', which is the path set during the class construction or an execution of 'set_data'.
            
            
        Returns
        -------
        Void
        """
        if input_param_file is not None:
            self._properties["input_param_file"] = os.path.abspath(input_param_file)
        if results_path is not None:
            self._side_properties["results_path"] = os.path.abspath(results_path)
        
        os.chdir(self.results_path)
        cmd = "{}/source/zphota -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
        subprocess.run(cmd.split())
    
    def _get_nb_filt_(self):
        """
        Return the number of filters included in the configuration file ("FILTER_LIST").
        
        
        Returns
        -------
        int
        """
        _, filt_list, _ = self._get_param_details_("FILTER_LIST")
        return len(filt_list.split(","))

    def _get_data_(self):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        nb_filt = self._get_nb_filt_()
        data = pandas.read_csv(self._get_param_details_("CAT_IN")[1], sep=" ")
        names = []
        for ii in np.arange(nb_filt):
            names.append("mag_{}".format(ii))
            names.append("mag_{}.err".format(ii))
        names += ["CONTEXT", "Z-SPEC"]
        names += ["STRING_{}".format(ii) for ii in np.arange(data.shape[1]-(2*nb_filt+3))]
        data = pandas.read_csv(self._get_param_details_("CAT_IN")[1], sep=" ", names=names)
        return data
    
    def fit_sed_errors(self, nb_fits=100):
        """
        
        
        Parameters
        ----------
        
        
        Options
        -------
        
        
        
        Returns
        -------
        
        """
        nb_filt = self._get_nb_filt_()
        data = self._get_data_()
        rand_mag = np.random.normal(loc=[data["mag_{}".format(ii)] for ii for np.arange(nb_filt)],
                                    scale=[data["mag_{}.err".format(ii)] for ii for np.arange(nb_filt)],
                                    size=(nb_filt, nb_fits))
        for ii in np.arange(nb_fits):
            for jj in np.arange(nb_filt):
                data["mag_{}".format(ii)] = rand_mag[jj, ii]
            data_path = pkg_resources.resource_filename(__name__, "config/")+"/data_buf.csv"
            data.to_csv(data_path, sep=" ", header=False)
            self.change_param("CAT_IN", os.path.abspath(data_path), False)
            results_path = pkg_resources.resource_filename(__name__, "results/")+"/data_buf.out"
            self.change_param("CAT_OUT", os.path.abspath(results_path))
            self._side_properties["results_path"] = os.path.abspath("/".join(results_path.split("/")[:-1]) + "/")
            self.run_zphota()
            
            sed_filename = "Id000000000.spec"
            idx_start = 10
            data_sed =  pandas.read_csv(os.path.expanduser(sed_dir+sed_filename),
                                        skiprows=idx_start, sep="  ", engine="python", nrows=nrows)
            while data_sed.shape[1] != 2:
                idx_start += 1
                data_buf =  pandas.read_csv(os.path.expanduser(sed_dir+sed_filename),
                                            skiprows=idx_start, sep="  ", engine="python", nrows=nrows)
            data_buf =  pandas.read_csv(os.path.expanduser(sed_dir+sed_filename),
                                        skiprows=idx_start, names=["lbda", "mag"], sep="  ",
                                        engine="python", nrows=nrows)
            
            if ii == 0:
                data_sed = data_buf
                data_sed.rename(columns={"mag":"mag0"}, inplace=True)
            else:
                data_sed["mag{}".format(ii)] = data_buf["mag"]
        
        data_sed.set_index("lbda", inplace=True)
        data_sed["mean"] = np.mean(data_sed, axis=1)
        data_sed["scale"] = np.std(data_sed, axis=1)
        return
    
    

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
    def results_path(self):
        """ Path of the result path. """
        return self._side_properties["results_path"]


