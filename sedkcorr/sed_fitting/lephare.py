
import numpy as np
import pandas
import os
import glob
import subprocess
import pkg_resources
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)

from propobject import BaseObject
from astrobject import photometry
from astrobject.instruments import sdss

from ..utils import tools





class BC03Installer( ):
    """ Test if the proper installation is made and enables to run it 

    USAGE:
    simply do: BC03Installer.install()

    """
    LEPHAREDIR = os.getenv("LEPHAREDIR",default=None)  # os.getenv('SEDM_USER',default="auto")
    
    # ----------- #
    # Initialize  #
    # ----------- #
    def __init__(self):
        """ """
    
    @classmethod
    def install(cls, compiler="gfortran", reinstall=False, verbose=False):
        """ """
        this = cls()
        if not reinstall and this.has_ised_files():
            if verbose:
                print("No need to install ised files. set reinstall to True to force reinstallation")
            return None
        
        this.run_ised_installer(compiler="gfortran", verbose=verbose)
        return
    
    # =========== #
    #  Methods    #
    # =========== #
    def get_list_sedmodel(self, which="*"):
        """ """
        return [l for l in os.listdir(self.BC03_CHAB) if l.startswith("bc2003")
               and ((which in ["*", "all", "both"]) or 
                    ((which.lower() == "ascii") and l.endswith("ASCII")) or
                    ((which.lower() == "ised") and l.endswith("ised"))
                   )]

    def build_bc03_installer(self, compiler="gfortran"):
        """ """
        os.system("%s  -O5 -Wall -ffixed-line-length-132 %s -o %s"%(compiler,
                                                                    self.bc03_installer+".f",
                                                                    self.bc03_installer))
        
    def run_ised_installer(self, compiler="gfortran", verbose=False):
        """ """
        if not self.has_bc03_installer():
            if verbose:
                print("No bc03_installer. Creating one, using the '%s' compiler"%compiler)
            self.build_bc03_installer(compiler)
            if not self.has_bc03_installer():
                raise AttributeError("build_bc03_installer() ran but still no has_bc03_installer()... sorry")
        
        if verbose:
            print("Creating the .ised files from the .ised_ASCII ones")
                  
        for file in self.get_list_sedmodel('ASCII'):
            report = os.system("%s %s"%(self.bc03_installer, self.BC03_CHAB+"/"+file))
     
    # ----------- #
    # Has tests   #
    # ----------- #
    def has_lephare_env(self):
        """ test that the global LEPHAREDIR is defined"""
        return self.LEPHAREDIR is not None

    def has_bc03_data(self):
        """ """
        if not self.has_lephare_env():
            raise IOError("$LEPHAREDIR is not defined ")
        
        return os.path.isdir( self.BC03_CHAB )
    
    def has_bc03_installer(self):
        """ """
        if not self.has_bc03_data():
            raise IOError("BC03_CHAB lephare gal library not downloaded. See http://www.cfht.hawaii.edu/~arnouts/LEPHARE/install.html")
            
        return os.path.isfile( self.bc03_installer)
    
    def has_ised_files(self):
        """ """
        return len(self.get_list_sedmodel("ised"))>0
    
    # =========== #
    #  Properties #
    # =========== #
    @property
    def BC03_CHAB(self):
        """ """
        return self.LEPHAREDIR+ "/sed/GAL/BC03_CHAB"
    
    @property
    def bc03_installer(self):
        """ """
        return self.BC03_CHAB+ "/bin_ised"
    

    






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
            
        filters : [list(string) or None]
            List of filters of the given measurements, the 'FILTER_LIST' parameter in the configuration file will be changed to the corresponding list of LePhare file names.
            The filter syntax must be like 'project.band' (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).
            If None, the filter list will be based on the configuration file.
        
        input_param_file : [string or None]
            Path of the input parameter file.
            If the given file doesn't exist yet, it will copied from the package one to the given file name.
            If 'None', the default file is imported from the package ('/config').
        
        output_param_file : [string or None]
            Path of the output parameter file.
            If the given file doesn't exist yet, it will copied from the package one to the given file name.
            If 'None', the default file is imported from the package ('/config').
        
        results_path : [string or None]
            Path for the results of the SED fitter.
            If 'None', the default folder is located in the package ('/results').
        
        flux_unit : [string]
            If 'data' is in flux, you can precise here the unit :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
            As LePhare needs the flux to be in the "Hz" (see above) unit, it will be converted if necessary.
        
        data_filename : [string]
            If your input data are in a dict/pandas.DataFrame, this option is used to change the file created to save these data for LePhare use.
            You can either enter a full path, or a simple file name without any "/".
            The latter will create the file directly in the package (less easy to find).
        
        
        Returns
        -------
        Void
        """
        if data is not None:
            self.set_data(data, **kwargs)

    def set_data(self, data=None, filters=None, input_param_file=None, output_param_file=None,
                 results_path=None, flux_unit="Hz", data_filename="data.csv", **kwargs):
        """
        Set up the file paths about the data, the config files (input and output) and the results path.
        
        Parameters
        ----------
        data : [string or pandas.DataFrame or dict]
            Path of the data file or a DataFrame/dict, both of them in a format readable by LePhare fitter.
            
        filters : [list(string) or None]
            List of filters of the given measurements, the 'FILTER_LIST' parameter in the configuration file will be changed to the corresponding list of LePhare file names.
            The filter syntax must be like "project.band" (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).
            If None, the filter list will be based on the configuration file.
        
        flux_unit : [string]
            If 'data' is in flux, you can precise here the unit :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
            As LePhare needs the flux to be in the "Hz" (see above) unit, it will be converted if necessary.
        
        Options
        -------
        input_param_file : [string or None]
            Path of the input parameter file.
            If the given file doesn't exist yet, it will copied from the package one to the given file name.
            If 'None', the default file is imported from the package ('/config/lephare_zphot_input.para').
        
        output_param_file : [string or None]
            Path of the output parameter file.
            If the given file doesn't exist yet, it will copied from the package one to the given file name.
            If 'None', the default file is imported from the package ('/config/lephare_zphot_output.para').
        
        results_path : [string or None]
            Path for the results of the SED fitter.
            If 'None', the default file is located in the package ('/results/data.out').
        
        data_filename : [string]
            If your input data are in a dict/pandas.DataFrame, this option is used to change the file created to save these data for LePhare use.
            You can either enter a full path, or a simple file name without any "/".
            The latter will create the file directly in the package (less easy to find).
        
        
        Returns
        -------
        Void
        """
        if input_param_file is not None:
            input_param_file = os.path.abspath(input_param_file)
            if not os.path.isfile(input_param_file):
                subprocess.run("cp {} {}".format(self.config_path+"/lephare_zphot_input.para", input_param_file), shell=True)
        else:
            input_param_file = self.config_path+"/lephare_zphot_input.para"
        
        if output_param_file is not None:
            output_param_file = os.path.abspath(output_param_file)
            if not os.path.isfile(output_param_file):
                subprocess.run("cp {} {}".format(self.config_path+"/lephare_zphot_ouput.para", output_param_file), shell=True)
        else:
            output_param_file = self.config_path+"/lephare_zphot_output.para"
        self._properties["input_param_file"] = input_param_file
        self._properties["output_param_file"] = output_param_file
        
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

        if filters is not None:
            print("""WARNING!!! Don't forget to check the "ERR_SCALE" parameter from the configuration file.\n""",
                  """You can check it with '.describe_params("ERR_SCALE")'. It must has the same size than the given 'filters'.""")
            self._set_filt_list_(filters)

        self.change_param("PARA_OUT", self.output_param_file, False)
        self._set_results_path_(results_path)
            
    def _set_input_type_(self):
        """
        Automatically changes the input configuration file 'INP_TYPE' parameter either the input data or magnitudes or flux.
        
        
        Returns
        -------
        Void
        """
        if self._get_param_details_("INP_TYPE")[1] == "M" and self.data_meas.iloc[0][self.data_meas.keys()[0]] < 1:
            self.change_param("INP_TYPE", "F", False)
        elif self._get_param_details_("INP_TYPE")[1] == "F" and self.data_meas.iloc[0][self.data_meas.keys()[0]] > 1:
            self.change_param("INP_TYPE", "M", False)

    def _convert_flux_(self, flux_unit):
        """
        This method automatically converts the input flux in erg.cm**-2.s**-1.Hz**-1, which necessary for LePhare.
        
        Parameters
        ----------
        flux_unit : [string]
            Input flux unit :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        
        Returns
        -------
        Void
        """
        data_buf = self.data_meas
        if self._get_param_details_("INP_TYPE")[1] == "F" and flux_unit != "Hz":
            for _filt in self.filt_list:
                data_buf["flux_"+_filt] = tools.convert_flux_unit(data_buf["flux_"+_filt], tools.FILTER_BANDS[_filt]["lbda"], flux_unit, "Hz")
                data_buf["flux_"+_filt+".err"] = tools.convert_flux_unit(data_buf["flux_"+_filt+".err"], tools.FILTER_BANDS[_filt]["lbda"], flux_unit, "Hz")
        data_buf.to_csv(self._get_param_details_("CAT_IN")[1], sep=" ", header=False)

    def _set_filt_list_(self, filters=["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"]):
        """
        This method changes the 'FILTER_LIST' parameter in the configuration file to correspond with the input 'filters' list.
        
        Parameters
        ----------
        filters : [string or list(string)]
            List of filters of the given measurements.
            The filter syntax must be like "project.band" (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).
        
        
        Returns
        -------
        Void
        """
        filters = filters if type(filters)==list else [filters]
        try:
            lp_filt_list = [tools.FILTER_BANDS[_filt]["lephare_name"] for _filt in filters]
        except(KeyError):
            raise KeyError("'filters' must be a list containing the filters with the syntax 'project.band' (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).")
        if lp_filt_list != self._get_param_details_("FILTER_LIST")[1]:
            self.change_param("FILTER_LIST", lp_filt_list, False)
        self._rename_libs_(filters=filters)
        
    def _rename_libs_(self, filters=None):
        """
        Rename the "LIB" file names in the configuration file, based on the given filters.
        
        Parameters
        ----------
        filters : [string or list(string)]
            Filter(s) used to fit the SED.
            The syntax must be "project.band" (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).
            
        
        Returns
        -------
        Void
        """
        if filters is not None:
            project_list = self._get_dict_project_(filters=filters).keys()
            gal_lib = self._get_param_details_("GAL_SED")[1].split("/")[-1].split("_")[0]
            proj_name = "_".join(project_list)
            change_params={"FILTER_FILE":"{}.filt".format(proj_name),
                           "STAR_LIB":"LIB_STAR_{}".format(proj_name),
                           "QSO_LIB":"LIB_QSO_{}".format(proj_name),
                           "GAL_LIB":"LIB_{}_{}".format(gal_lib, proj_name),
                           "STAR_LIB_OUT":"STAR_{}".format(proj_name),
                           "QSO_LIB_OUT":"QSO_{}".format(proj_name),
                           "GAL_LIB_OUT":"{}_{}".format(gal_lib, proj_name)}
            for k, v in change_params.items():
                self.change_param(k, v, False)
        else:
            raise ValueError("'filters' is None.")

    def _get_dict_project_(self, filters=None):
        """
        Return a dictionary containing, for each project, the list of corresponding bands.
        
        Parameters
        ----------
        filters : [string or list(string)]
            Filter(s) to sort in a dictionary.
            The syntax must be "project.band".
        
        
        Returns
        -------
        dict
        """
        dict_project = {}
        for _filt in filters if type(filters)==list else [filters]:
            _project, _band = _filt.split(".")
            try:
                dict_project[_project].append(_band)
            except(KeyError):
                dict_project[_project] = [_band]
        return dict_project

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
        elif len(results_path.split("/")[-1].split(".")) < 2:
            raise ValueError("You must input a compatible file name. Your input is '{}'.".format(results_path))
        
        self.change_param("CAT_OUT", os.path.abspath(results_path))
        self._side_properties["results_path"] = os.path.abspath("/".join(results_path.split("/")[:-1])) + "/"

    def _get_results_path_(self):
        """
        Return the full path of the LePhare results file.
        
        
        Returns
        -------
        string
        """
        return self._get_param_details_("CAT_OUT")[1]

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
        if param == "FILTER_LIST" and new_param_value[0] in tools.FILTER_BANDS.keys():
            self._set_filt_list_(new_param_value)
            return
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
        
        if param in ["{}_LIB".format(elt) for elt in ["GAL", "QSO", "STAR"]]:
            self.change_param(param+"_IN", new_param_value, False)
        elif param in ["{}_LIB_OUT".format(elt) for elt in ["GAL", "QSO", "STAR"]]:
            self.change_param("ZPHOTLIB", [self._get_param_details_("{}_LIB_OUT".format(elt))[1] for elt in ["GAL", "QSO", "STAR"]], False)
        elif param == "CAT_OUT":
            self._side_properties["results_path"] = os.path.abspath("/".join(new_param_value.split("/")[:-1]) + "/")
        elif param == "GAL_SED":
            self._rename_libs_(filters=self.filt_list)

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
            If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, 
            and as values a list containing the new parameter value as first element and the force comment option as second.
        
        
        Returns
        -------
        Void
        """
        self._init_changes_(input_param_file=input_param_file, results_path=None, change_params=change_params)
                
        filt_k, filt_v, _ = self._get_param_details_("FILTER_FILE")
        
        if not os.path.isfile(self.PATH_LEPHAREWORK+"/filt/"+filt_v) or update:
            cmd = "{}/source/filter -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
            try:
                subprocess.run(cmd.split())
            except:
                raise ValueError("LePhareError : unable to run 'filter'.")
    
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
                try:
                    subprocess.run(cmd.split())
                except:
                    raise ValueError("LePhareError : unable to run 'sedtolib' for '{}'.".format(elt))
    
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
            try:
                subprocess.run(cmd.split())
            except:
                raise ValueError("LePhareError : unable to run 'mag_star'.")
    
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
                try:
                    subprocess.run(cmd.split())
                except:
                    raise ValueError("LePhareError : unable to run 'mag_gal' for '{}'.".format(elt))
    
    def run_zphota(self, filters=None, input_param_file=None, output_param_file=None, results_path=None, change_params=None,
                   savefile=None, del_spec=True, **kwargs):
        """
        First change current directory to the results path.
        Then execute "$LEPHAREDIR/source/zphota -c [...].para" in the shell.
        
        Options
        -------
        filters : [list(string) or dict or None]
            If a list of filters is given, the context will be changed for every line to the corresponding one.
            If you want to change only a few indexes context, you can give a dictionary containing a list of the indexes under the key "id" and the list of filters under the key "filters".
            If None, nothing changes from the 'data_meas'.
            The filter syntax must be like 'project.band' (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).
        
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        output_param_file : [string or None]
            If you want to set a new output parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        results_path : [string or None]
            If you want to set a new results path, give the new path here.
            Default is 'None', which is the path set during the class construction or an execution of 'set_data'.
            
        change_params : [dict or None]
            If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, and as values a list containing the new parameter value as first element and the force comment option as second.
        
        savefile : [string or None]
            If not None, the 'data_sed' will be saved in the given file path.
        
        del_spec : [bool]
            If True, delete the spectra files after saving data in 'data_sed'.
            Default is True.
            
            
        Returns
        -------
        Void
        """
        if output_param_file is not None:
            self._properties["output_param_file"] = os.path.abspath(output_param_file)
            self.change_param("PARA_OUT", self.output_param_file, False)
        
        self._init_changes_(input_param_file=input_param_file, results_path=results_path, change_params=change_params)
        
        if filters is not None:
            if type(filters) == list:
                self.set_filters_context(id=None, filters=filters)
            elif type(filters) == dict:
                self.set_filters_context(**filters)
        
        os.chdir(self.results_path)
        cmd = "{}/source/zphota -c {}".format(self.PATH_LEPHAREDIR, self.input_param_file)
        try:
            subprocess.run(cmd.split())
        except:
            raise ValueError("LePhareError : unable to run 'zphota'.")
        
        self.set_data_sed()
        self.set_data_res()
        if del_spec:
            for _spec in glob.glob(self.results_path+"*.spec"):
                os.remove(_spec)
        if savefile is not None:
            self.write(savefile, None)

    def run_init(self, input_param_file=None, update=False, change_params=None, **extras):
        """
        Run shell commands to initialize LePhare fitting.
        
        Options
        -------
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        update : [bool]
            Set to True if you want to update the initialization on filters, sed libraries, etc.
        
        change_params : [dict or None]
            If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, and as values a list containing the new parameter value as first element and the force comment option as second.
        
        
        Returns
        -------
        Void
        """
        self._init_changes_(input_param_file=input_param_file, change_params=change_params)
    
        self.run_filter(input_param_file=None, update=update, change_params=None)
        self.run_sedtolib(input_param_file=None, update=update, change_params=None)
        self.run_mag_star(input_param_file=None, update=update, change_params=None)
        self.run_mag_gal(input_param_file=None, update=update, change_params=None)

    def run_fit(self, filters=None, input_param_file=None, output_param_file=None, results_path=None, update=False, change_params=None, savefile=None, **kwargs):
        """
        Run shell commands to execute LePhare fitting.
        
        Options
        -------
        filters : [list(string) or dict or None]
            If a list of filters is given, the context will be changed for every line to the corresponding one.
            If you want to change only a few indexes context, you can give a dictionary containing a list of the indexes under the key "id" and the list of filters under the key "filters".
            If None, nothing changes from the 'data_meas'.
            The filter syntax must be like 'project.band' (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).
        
        input_param_file : [string or None]
            If you want to set a new input parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
            
        output_param_file : [string or None]
            If you want to set a new output parameter file, give the new path here.
            Default is 'None', which is the file set during the class construction or an execution of 'set_data'.
        
        results_path : [string or None]
            If you want to set a new results path, give the new path here.
            Default is 'None', which is the path set during the class construction or an execution of 'set_data'.
            
        update : [bool]
            Set to True if you want to update the initialization on filters, sed libraries, etc.
            
        change_params : [dict or None]
            If you want to change any configuration parameters, put a dictionary with parameters you want to change as keys, and as values a list containing the new parameter value as first element and the force comment option as second.
        
        savefile : [string or None]
            If not None, the 'data_sed' will be saved in the given file path.
        
        **kwargs
        del_spec : [bool]
            If True, delete the spectra files after saving data in 'data_sed'.
            Default is True.
        
        
        Returns
        -------
        Void
        """
        self._init_changes_(input_param_file=input_param_file, results_path=results_path, change_params=change_params)
        
        self.run_init(input_param_file=None, update=update, change_params=None)
        self.run_zphota(filters=filters, input_param_file=None, output_param_file=output_param_file, results_path=None,
                        change_params=None, savefile=savefile, **kwargs)
    
    def _get_filt_list_(self):
        """
        Return the list of filters corresponding to the LePhare filter files included in the configuration file ("FILTER_LIST").
        
        
        Returns
        -------
        list(string)
        """
        filt_list = []
        lp_filt_list = self._get_param_details_("FILTER_LIST")[1]
        for _filt in lp_filt_list.split(","):
            for key, value in tools.FILTER_BANDS.items():
                if value["lephare_name"] == _filt:
                    filt_list.append(key)
                    break
                else:
                    continue
                raise ValueError("{} is an unknown lephare filter file name.".format(_filt))
        return filt_list #[self.lephare_filt_to_filt(tools.FILTER_BANDS, _filt) for _filt in filt_list.split(",")]
    
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
        
        Parameters
        ----------
        data : [pandas.DataFrame]
            Input data in a LePhare format.
        
        
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

    def set_filters_context(self, id=None, filters=["sdss.u", "sdss.g", "sdss.r", "sdss.i", "sdss.z"]):
        """
        Automatically changes the 'CONTEXT' value in 'data_meas' given the list filters to fit on.
        
        Parameters
        ----------
        id : [int or list(int) or None]
            Index(es) of the line(s) you want to change the context by the given filters.
            If None, the context is changed for every line.
        
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
        cat_fmt = self._get_param_details_("CAT_FMT")[1]
        names = []
        for _filt in self.filt_list:
            names.append("{}_{}".format(prefix, _filt))
            if cat_fmt == "MEME":
                names.append("{}_{}.err".format(prefix, _filt))
        if cat_fmt == "MMEE":
            for _filt in self.filt_list:
                names.append("{}_{}.err".format(prefix, _filt))
        names += ["CONTEXT", "Z-SPEC"]
        names += ["STRING_{}".format(ii) for ii in np.arange(data.shape[1]-(2*len(self.filt_list)+3))]
        data = pandas.read_csv(self._get_param_details_("CAT_IN")[1], sep=" ", names=names)
        return data
    
    def _get_header_spec_file_(self):
        """
        Return the line index to start reading the specrtrum data in LePhare .spec files.
        
        
        Returns
        -------
        int
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
    
    def _read_spec_(self, spec_filename):
        """
        Return a DataFrame of the fitted SED (wavelength, magnitude).
        
        Parameters
        ----------
        sed_filename : [string]
            File path of the LePhare spectrum file.
        
        
        Returns
        -------
        pandas.DataFrame
        """
        data_sed =  pandas.read_csv(os.path.expanduser(spec_filename),
                                    skiprows=self.header_spec_file, names=["lbda", "mag"], sep="  ",
                                    engine="python", nrows=1000)
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

    def set_data_sed(self, rm_fails=False):
        """
        Set the LePhare fit results in a dictionary containing the fitted spectrum tables.
        
        Options
        -------
        rm_fails : [bool]
            If True, remove the fit fails of the SEDs.
        
        
        Returns
        -------
        Void
        """
        jj = 0
        for ii in np.arange(len(self.data_meas)):
            if rm_fails and ii in self._get_idx_failed_():
                jj += 1
                continue
            self.data_sed[ii-jj] = self._read_spec_(self._get_sed_filename_(ii))
    
    def set_data_res(self, rm_fails=False):
        """
        Set the LePhare outfile as a pandas table attribute.
        
        Options
        -------
        rm_fails : [bool]
            If True, remove the fit fails of the SEDs and the results file.
        
        
        Returns
        -------
        Void
        """
        self._derived_properties["data_res"] = self.lephare_output_file_reader(filename=self._get_param_details_("CAT_OUT")[1], filters=self.filt_list)
        if rm_fails:
            self.data_res.drop(index=self._get_idx_failed_(), inplace=True)
    
    def _get_idx_failed_(self):
        """
        Return a list of the indexes for which the fit has failed.
        
        
        Returns
        -------
        list(int)
        """
        return [ii for ii, v in self.data_res.iterrows() if v["Z_BEST"] == -99]
    
    def _rm_fails_(self):
        """
        Remove the failed cases from the results file.
        
        
        Returns
        -------
        Void
        """
        

    def show(self, ax=None, id_sed=0, y_unit="AA", plot_sed=True, plot_phot=True, xlim=(None, None), ylim=(None, None),
             xscale="linear", yscale="linear", savefile=None, **kwargs):
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
        ax : [matplotlib.axes or None]
            Already existing axes you want to add stuff in.
            Default is None.
        
        plot_sed : [bool]
            If True, plot the fitted SED.
            Default is True.
        
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
        mask = np.array([True] * len(self.data_sed[id_sed]))
        if xlim != (None, None):
            mask = (xlim[0] < np.asarray(self.data_sed[id_sed]["lbda"])) if xlim[0] is not None else 1
            mask *= (np.asarray(self.data_sed[id_sed]["lbda"]) < xlim[1]) if xlim[1] is not None else 1
        
        x_sed = np.asarray(self.data_sed[id_sed]["lbda"])[mask]
        y_sed = np.asarray(self.data_sed[id_sed]["mag"])[mask]
        if y_unit in ["Hz", "AA", "mgy"]:
            y_sed, _ = tools.mag_to_flux(y_sed, None, band=x_sed, flux_unit=y_unit, opt_mAB0=False)
        
        #SED
        opt_sed = {"ls":"-", "marker":"", "color":"0.4"}
        if plot_sed:
            ax.plot(x_sed, y_sed, label="_nolegend_", **{**opt_sed, **kwargs})
        
        #Photometry
        if plot_phot:
            prefix = "mag" if self._get_param_details_("INP_TYPE")[1] == "M" else "flux"
            data_meas = self.data_meas if id_sed != -1 else self.data_orig
            id_sed = 0 if id_sed is None else id_sed
            for _filt in self._get_context_filters_(self.data_meas.iloc[id_sed]):
                x_phot = tools.FILTER_BANDS[_filt]["lbda"]
                y_phot = float(data_meas.iloc[id_sed]["{}_{}".format(prefix, _filt)])
                y_phot_err = float(data_meas.iloc[id_sed]["{}_{}.err".format(prefix, _filt)])
                if y_unit in ["Hz", "AA", "mgy"] and prefix == "mag":
                    y_phot, y_phot_err = tools.mag_to_flux(y_phot, y_phot_err, band=_filt, flux_unit=y_unit, opt_mAB0=True)
                elif y_unit == "mag" and prefix == "flux":
                    y_phot, y_phot_err = tools.flux_to_mag(y_phot, y_phot_err, band=_filt, flux_unit="Hz", opt_mAB0=True)
                elif y_unit in ["Hz", "AA", "mgy"] and prefix == "flux":
                    y_phot, y_phot_err = tools.convert_flux_unit((y_phot, y_phot_err), tools.FILTER_BANDS[_filt]["lbda"], "Hz", y_unit)
                ax.errorbar(x_phot, y_phot, yerr=y_phot_err, ls="", marker="o", color=tools.FILTER_BANDS[_filt]["color"], label=_filt)
        
        #Writings
        ax.set_xlabel(r"$\lambda$ [\AA]", fontsize="large")
        ylabel = "mgy" if y_unit == "mgy" else \
                 r"${{f}}_{{\nu}}$ $[erg.{{s}}^{{-1}}.{{cm}}^{{-2}}.{Hz}^{{-1}}]$" if y_unit == "Hz" else \
                 r"${{f}}_{{\lambda}}$ $[erg.{{s}}^{{-1}}.{{cm}}^{{-2}}.{\AA}^{{-1}}]$" if y_unit == "AA" else \
                 "mag"
        ax.set_ylabel(ylabel, fontsize="large")
        
        if plot_phot:
            ax.legend(loc="upper right", ncol=1)

        #Fig view
        ax.set_xlim(xlim)
        if ylim == (None, None):
            if ax.get_ylim()[0] <= 0:
                ylim = (np.min(y_sed), ylim[1])
        ax.set_ylim(ylim)

        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        
        #Save
        if savefile is not None:
            fig.savefig(savefile)
        
        return {"ax":ax, "fig":fig} if id_sed != -1 else ({"ax":ax, "fig":fig}, y_sed)
    
    def get_data_sed(self, id_sed=0):
        """
        Return a dataframe of SED spectr.um.a.
        
        Parameters
        ----------
        id_sed : [int or None]
            Index of the 'data_sed' dictionary you want to save in a file.
            If None, the whole dictionary will be saved, each column name containing the id of the income.
        
        
        Returns
        -------
        pandas.DataFrame
        """
        if type(id_sed) == int:
            return self.data_sed[id_sed]
        elif id_sed is None:
            data_out = self.data_sed[0].copy()
            data_out.set_index("lbda", inplace=True)
            data_out.rename(columns={k:k+"_id_0" for k in data_out.keys()}, inplace=True)
            for k1, v1 in self.data_sed.items():
                if k1 == 0 or k1 == -1:
                    continue
                pd_buf = v1.set_index("lbda").copy()
                for k2, v2 in pd_buf.items():
                    data_out[k2+"_id_{}".format(k1)] = v2
            data_out.reset_index(inplace=True)
            return data_out
        raise ValueError("'id_sed' must be either an integer index contained in the 'data_sed' dictionary, or None if you want to save everything.")

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
        data_out = self.get_data_sed(id_sed=id_sed)
        data_out.to_csv(savefile, sep=" ", index=False)




    
    
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
    def lephare_output_file_reader(filename, filters=None):
        """
        Read a LePhare output file, returning a pandas.DataFrame.

        Parameters
        ----------
        filename : [string]
            Path of the LePhare file.
        
        Options
        -------
        filters : [list(string) or None]
            List of filters you want to use to name the corresponding columns.
            If None, use the LePhare filter names, included in the output file.


        Returns
        -------
        pandas.DataFrame
        """
        with open(filename, "r") as f1:
            buf_file = f1.readlines()
        buf_delimiter = buf_file[0]
        
        #Filters
        for ii, line in enumerate(buf_file):
            if "FILTER_FILE" in line:
                ii_filters = ii + 1
                break
        filter_list = buf_file[ii_filters].split()[2:]
        if filters is not None:
            if len(filters) != len(filter_list):
                raise ValueError("You must input {} filter names. You've input {} of them.".format(len(filter_list), len(filters)))
            else:
                filter_list = filters

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
    
    @staticmethod
    def set_input_data(data, sn_name, sn_z, filters):
        """
        Format input data (with fluxes) in a LePhare compatible format.
        
        Parameters
        ----------
        data : [dict or pandas.DataFrame]
            Input data.
        
        sn_name : [string]
            SNeIa's name.
        
        sn_z : [float]
            SNeIa's redshift.
        
        filters : [list(string)]
            List of filters to fit in LePhare.
        
        
        Returns
        -------
        pandas.DataFrame
        """
        data_out = {}
        for _filt in [filters] if len(np.array([filters]).shape) == 1 else filters:
            data_out[_filt+"_flux"] = [data[_filt]["flux"]]
            data_out[_filt+"_flux.err"] = [data[_filt]["flux.err_up"]]
        data_out["CONTEXT"] = [np.sum([2**ii for ii, _ in enumerate(filters)])]
        data_out["Z-SPEC"] = [sn_z]
        data_out["STRING"] = [sn_name]
        return pandas.DataFrame(data_out)
    
    
    
    
    
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
        if self._derived_properties["data_res"] is None:
            self.set_data_res()
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
            
        filters : [list(string) or None]
            List of filters of the given measurements, the 'FILTER_LIST' parameter in the configuration file will be changed to the corresponding list of LePhare file names.
            The filter syntax must be like 'project.band' (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).
            If None, the filter list will be based on the configuration file.
        
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
        
        data_filename : [string]
            This option is used to change the file name created to save the Monte Carlo data for LePhare use.
            You can either enter a full path, or a simple file name without any "/".
            The latter will create the file directly in the package (less easy to find).
        
        
        Returns
        -------
        Void
        """
        if data is not None:
            self.set_data(data, **kwargs)

    def set_data(self, data=None, filters=None, input_param_file=None, output_param_file=None, results_path=None,
                 flux_unit="Hz", data_filename="data_rand.csv", nb_draw=100, **kwargs):
        """
        Set up the file paths about the data, the config files (input and output) and the results path.
        
        Parameters
        ----------
        data : [string or pandas.DataFrame or dict]
            Path of the data file or a DataFrame/dict, both of them in a format readable by LePhare fitter.
            Only accept the data for one by one object.
        
        filters : [list(string) or None]
            List of filters of the given measurements, the 'FILTER_LIST' parameter in the configuration file will be changed to the corresponding list of LePhare file names.
            The filter syntax must be like 'project.band' (ex: 'sdss.r', 'galex.FUV', 'ps1.g', ...).
            If None, the filter list will be based on the configuration file.
        
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
        
        Options
        -------
        data_filename : [string]
            This option is used to change the file name created to save the Monte Carlo data for LePhare use.
            You can either enter a full path, or a simple file name without any "/".
            The latter will create the file directly in the package (less easy to find).
        
        
        Returns
        -------
        Void
        """
        _ = super(LePhareRand, self).set_data(data=data, filters=filters, input_param_file=input_param_file, output_param_file=output_param_file,
                                              results_path=results_path, flux_unit=flux_unit, data_filename=data_filename, **kwargs)
        self.set_rand_data(nb_draw=nb_draw, data_filename=data_filename)
    
    def set_rand_data(self, nb_draw=100, data_filename="data_rand.csv"):
        """
        Initialize the Monte Carlo random data to fit through LePhare in order to get SED errors.
        
        Parameters
        ----------
        nb_draw : [int]
            Number of draw for the Monte Carlo fitting errors on the SED.
        
        Options
        -------
        data_filename : [string]
            This option is used to change the file name created to save the Monte Carlo data for LePhare use.
            You can either enter a full path, or a simple file name without any "/".
            The latter will create the file directly in the package (less easy to find).
        
        
        Returns
        -------
        Void
        """
        self._properties["data_orig"] = self.data_meas
        for _filt in self.filt_list:
            photopoint = photometry.PhotoPoint(flux=self.data_meas["flux_{}".format(_filt)],
                                               var=self.data_meas["flux_{}.err".format(_filt)]**2,
                                               lbda=tools.FILTER_BANDS[_filt]["lbda"])
            photopoint.draw_photosamplers(nsamplers=nb_draw, negative_fluxmag=40)
            self.data_rand["flux_{}".format(_filt)] = photopoint.photosamplers.samplers
            self.data_rand["flux_{}.err".format(_filt)] = np.array([self.data_meas["flux_{}.err".format(_filt)].values[0]]*nb_draw)
        self.data_rand["CONTEXT"] = np.array([self.data_meas["CONTEXT"].values[0]]*nb_draw)
        self.data_rand["Z-SPEC"] = np.array([self.data_meas["Z-SPEC"].values[0]]*nb_draw)
        for _col_name in self.data_meas.columns:
            if "STRING" in _col_name:
                self.data_rand[_col_name] = np.array([self.data_meas[_col_name].values[0]]*nb_draw)
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
            - "Hz" : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" [default] : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
        
        
        Returns
        -------
        Void
        """
        lbda = self.data_sed[0]["lbda"]
        mags = np.quantile([v["mag"] for k, v in self.data_sed.items() if (k != -1 and len(v) != 0)], quants, axis=0)
        return {k:m if y_unit=="mag" else tools.mag_to_flux(m, None, band=lbda, flux_unit=y_unit, opt_mAB0=False)[0]
                for k, m in zip(quants, mags)}
    
    def set_data_sed(self):
        """
        Set the LePhare fit results in a dictionary.
        Add the Monte Carlo results with a dictionary containing the wavelength, the spectra median and one sigma error below and above the median.
        
        
        Returns
        -------
        Void
        """
        _ = super(LePhareRand, self).set_data_sed(rm_fails=True)
        quants = self._get_fit_quantiles_(quants=[0.16, 0.5, 0.84], y_unit="mag")
        pd_none = {"lbda":self.data_sed[0]["lbda"], "mag":quants[0.5], "mag.err_low":quants[0.5]-quants[0.16], "mag.err_up":quants[0.84]-quants[0.5]}
        self.data_sed[-1] = pandas.DataFrame(pd_none)
        self.set_data_res(rm_fails=True)

    def show(self, ax=None, id_sed=None, y_unit="AA", plot_sed=True, plot_phot=True, xlim=(None, None), ylim=(None, None),
             xscale="linear", yscale="linear", savefile=None, show_sigmas=[1, 2], **kwargs):
        """
        Plot method.
        Return dict("fig", "ax").
        
        Parameters
        ----------
        id_sed : [int or None]
            Index of the SED you want to plot, corresponding to index in the 'data_meas' table.
            If None, plot the Monte Carlo fitted median of the spectra, and sigmas if specified (see 'show_sigmas').
            Default is None.
        
        y_unit : [string]
            Choice to plot "mag" or flux with :
            - "Hz" [default] : erg . cm**-2 . s**-1 . Hz**-1
            - "AA" : erg . cm**-2 . s**-1 . AA**-1 (AA = Angstrom)
            - "mgy" : mgy (mgy = maggies)
            Default is "AA".
        
        Options
        -------
        ax : [matplotlib.axes or None]
            Already existing axes you want to add stuff in.
            Default is None.
        
        plot_sed : [bool]
            If True, plot the fitted SED.
            Default is True.
        
        plot_phot : [bool]
            If True, plot the photometry points with errors, either in flux or magnitude.
            Default is True.
        
        xlim : [tuple[float or None]]
            Set the limits on the x axis.
            If (None, None), the figure has free x axis limits.
            Default is (None, None).
        
        ylim : [tuple[float or None]]
            Set the limits on the y axis.
            If (None, None), the figure has free y axis limits.
            Default is (None, None).
        
        xscale : [string]
            Scale of the x axis : "linear", "log", ...
            Default is "linear".
        
        yscale : [string]
            Scale of the y axis : "linear", "log", ...
            Default is "linear".
        
        savefile : [string or None]
            If None, the figure won't be saved.
            To save it, input a path directory + filename.
        
        show_sigmas : [int or list(int) or None]
            Show 1 and/or 2 (or None) sigmas error of the SED.
            Default is [1, 2]
        
        **kwargs : [dict]
            pyplot.plot options to apply on the SED spectrum.
        
        
        Returns
        -------
        dict
        """
        kwargs = {**{"zorder":3}, **kwargs}
        dict_fig, y_sed = super(LePhareRand, self).show(ax=ax, id_sed=id_sed, y_unit=y_unit, plot_sed=plot_sed, plot_phot=plot_phot,
                                                        xlim=xlim, ylim=ylim, xscale="linear", yscale="linear", savefile=None, **kwargs)
        
        if id_sed is None:
            lbda = self.data_sed[0]["lbda"]
            mask = np.array([True] * len(lbda))
            if xlim != (None, None):
                mask = (xlim[0] < np.asarray(lbda)) if xlim[0] is not None else 1
                mask *= (np.asarray(lbda) < xlim[1]) if xlim[1] is not None else 1
            nsigmas = len(np.atleast_1d(show_sigmas))
            if plot_sed:
                if 2 in show_sigmas:
                    ff = self._get_fit_quantiles_(quants=[0.05, 0.95], y_unit=y_unit)
                    dict_fig["ax"].fill_between(lbda[mask], np.asarray(ff[0.05])[mask], np.asarray(ff[0.95])[mask], alpha=0.3/nsigmas, color="C0", lw=0, zorder=1)
                if 1 in show_sigmas:
                    ff = self._get_fit_quantiles_(quants=[0.16, 0.84], y_unit=y_unit)
                    dict_fig["ax"].fill_between(lbda[mask], np.asarray(ff[0.16])[mask], np.asarray(ff[0.84])[mask], alpha=0.3/nsigmas, color="C0", lw=0, zorder=2)
    
        dict_fig["ax"].set_xlim(xlim)
        if ylim == (None, None):
            if dict_fig["ax"].get_ylim()[0] <= 0:
                ylim = (np.min(y_sed), ylim[1])
        dict_fig["ax"].set_ylim(ylim)

        dict_fig["ax"].set_xscale(xscale)
        dict_fig["ax"].set_yscale(yscale)
        
        if savefile is not None:
            dict_fig["fig"].savefig(savefile)
        
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

