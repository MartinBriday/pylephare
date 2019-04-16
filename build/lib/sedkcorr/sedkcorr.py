
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from . import sed_lephare
from propobject import BaseObject





class KCorrection( BaseObject ):
    """
    Kcorrection based on SED fits.
    """
    
    PROPERTIES         = ["meas_data", "sed_dir", "sed_data"]
    SIDE_PROPERTIES    = ["nb_sn"]
    DERIVED_PROPERTIES = ["kcorr_data"]
    
    def set_data(self, meas_dir=None, meas_data=None, col_syntax=["mag_band", "mag_band_err"], sed_dir=None):
        """
        Load measurements data and set the SED fits directory.
        
        Parameters
        ----------
        meas_dir : [string or None]
            Path of the measurements data.
            You can not use both meas_dir and meas_data.
        
        meas_data : [pandas.DataFrame or None]
            Measurements data input in a DataFrame.
            You can not use both meas_data and meas_dir.
        
        col_syntax : [list[string]]
            Syntax of measurements and errors column names in the meas_data table.
            Replace the filter band in the column names with the word "band" (eg: ["mag_band", "mag_band_err"]).
        
        sed_dir : [string]
            Path to the SED files folder.
            Needed if sed_index used.
        
        
        Returns
        -------
        Void
        """
        #Measurements table
        if meas_dir is not None and meas_data is None:
            self._properties["meas_data"] = pd.read_table(os.path.expanduser(meas_dir), sep="  ")
        elif meas_dir is None and meas_data is not None:
            self._properties["meas_data"] = meas_data
        else:
            raise ValueError("You must input one (and only one) data set option : the filename directory or a DataFrame.")
        
        #Spectra directory
        self._properties["sed_dir"] = os.path.expanduser(sed_dir)
        
        #SED data
        self.set_sed_data(col_syntax=col_syntax)
    
    def set_sed_data(self, col_syntax=["mag_band", "mag_band_err"]):
        """
        Create an SED object for every SNeIa from the meas_data table.
        
        Parameters
        ----------
        col_syntax : [list[string]]
            Syntax of measurements and errors column names in the meas_data table.
            Replace the filter band in the column names with the word "band" (eg: ["mag_band", "mag_band_err"]).
        
        
        Returns
        -------
        Void
        """
        for idx in self.sed_data:
            self.sed_data[idx].set_sed_data(sed_index=idx, sed_dir=self.sed_dir)
            self.sed_data[idx].set_meas_data(meas_data=self.meas_data.loc[idx], col_syntax=col_syntax)
        
    def k_correction(self, list_bands=["FUV", "NUV", "u", "g", "r", "i", "z"], save_data=False, save_filename=None, save_dir=None):
        """
        Apply k corrections to measurements.
        
        Parameters
        ----------
        list_bands : [list[string]]
            List of the filter bands to apply the k correction.
        
        Options
        -------
        save_data : [bool]
            If True, the data will be saved in the given directory.
        
        save_filename : [string]
            Name of the file to save the data.
        
        save_dir : [string or None]
            Directory where the data will be saved.
            If None, the default value is the Desktop (Unix users...).
        
        
        Returns
        -------
        Void
        """
        kcorr_data = {}
        for band in list_bands:
            kcorr_data["mag_"+band] = [self.sed_data[idx].kcorr_data[band]["mag"] for idx in range(self.nb_sn)]
            kcorr_data["mag_"+band+"_err"] = [self.sed_data[idx].kcorr_data[band]["mag.err"] for idx in range(self.nb_sn)]
        self._derived_properties["kcorr_data"] = pd.DataFrame(kcorr_data)
        self.kcorr_data["SN_name"] = self.meas_data["SN_name"]
        self.kcorr_data.set_index("SN_name", inplace=True)
        
        if save_data:
            save_dir = "~/Desktop/" if save_dir is None else save_dir
            directory_out = os.path.expanduser(save_dir)
            if save_filename is None or type(save_filename)!=str:
                raise ValueError("If you want to save the data, you must input a file name.")
            self.kcorr_data.to_csv(directory_out+save_filename, sep=" ")
        
        
    
    
    # ================ #
    #  Properties      #
    # ================ #
    @property
    def meas_data(self):
        """ DataFrame of the measurements """
        return self._properties["meas_data"]
    
    @property
    def sed_dir(self):
        """ Path of the SED files """
        return self._properties["sed_dir"]
    
    @property
    def sed_data(self):
        """ Dictionnary of SED objects for every SNeIa """
        if self._properties["sed_data"] is None:
            self._properties["sed_data"] = {idx:sed_lephare.SED() for idx in range(self.nb_sn)}
        return self._properties["sed_data"]
    
    @property
    def nb_sn(self):
        """ Number of SNeIa """
        if self._side_properties["nb_sn"] is None:
            self._side_properties["nb_sn"] = len(self.meas_data)
        return self._side_properties["nb_sn"]
    
    @property
    def kcorr_data(self):
        """ K corrected data """
        return self._derived_properties["kcorr_data"]













