
import pandas
import numpy as np
import matplotlib.pyplot as plt

from . import sed_lephare
from propobject import BaseObject





class KCorrection( BaseObject ):
    """
    Kcorrection based on SED fits.
    """
    
    PROPERTIES         = ["meas_data", "sed_dir", "sed_data"]
    SIDE_PROPERTIES    = ["nb_sn"]
    DERIVED_PROPERTIES = ["kcorr_data"]
    
    def set_data(self, meas_dir=None, meas_data=None, sed_dir=None):
        """
        
        """
        #Measurements table
        if meas_dir is not None and meas_data is None:
            self._properties["meas_data"] = pd.read_table(filename, sep="  ")
        elif meas_dir is None and meas_data is not None:
            self._properties["meas_data"] = meas_data
        else:
            raise ValueError("You must input one (and only one) data set option : the filename directory or a DataFrame.")
        
        #Spectra directory
        self._properties["sed_dir"] = sed_dir
    
    def set_sed_data(self, z=None, col_syntax=["mag_band", "mag_band_err"]):
        """
        
        """
        for idx in self.sed_data:
            self.sed_data[idx].set_sed_data(sed_index=idx, sed_dir=self.sed_dir)
            self.sed_data[idx].set_meas_data(meas_data=self.meas_data.loc[idx], z=z, col_syntax=col_syntax)
        
    def k_correction(self, list_bands=["FUV", "NUV", "u", "g", "r", "i", "z"]):
        """
        
        """
        kcorr_data = {}
        for band in list_bands:
            kcorr_data["mag_"+band] = [self.sed_data[idx].kcorr_data["mag"] for idx in range(self.nb_sn)]
            kcorr_data["mag_"+band+"_err"] = [self.sed_data[idx].kcorr_data["mag.err"] for idx in range(self.nb_sn)]
        self._derived_properties["kcorr_data"] = pd.DataFrame(kcorr_data)
        
        
    
    
    # ================ #
    #  Properties      #
    # ================ #
    @property
    def meas_data(self):
        """  """
        return self._properties["meas_data"]
    
    @property
    def sed_dir(self):
        """  """
        return self._properties["sed_dir"]
    
    @property
    def sed_data(self):
        """  """
        if self._properties["sed_data"] is None:
            self._properties["sed_data"] = {idx:sed_lephare.SED() for idx in range(self.nb_sn)}
        return self._properties["sed_data"]
    
    @property
    def nb_sn(self):
        """  """
        if _side_properties["nb_sn"] is None:
            _side_properties["nb_sn"] = len(self.data)
        return _side_properties["nb_sn"]
    
    @property
    def kcorr_data(self):
        """  """
        return self._derived_properties["kcorr_data"]













