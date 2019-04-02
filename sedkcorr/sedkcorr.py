
import pandas
import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from sncosmo import bandpasses

from pyifu import spectroscopy
from propobject import BaseObject



# ----------------- #
#                   #
#  Functions        #
#                   #
# ------------------#
def lbda_z0(lbda, z):
    return lbda/(1+z)

def flux_z0(flux, z):
    return flux*(1+z)**2

def mag_to_flux(mag):
    return np.power(10, (mag + 48.585)/(-2.5))

def mag_to_flux_err(mag, mag_err, mag0_err=0.005):
    return (mag_to_flux(mag)*np.log(10)/(-2.5))*np.sqrt(mag_err**2 + mag0_err**2)


LBD_INFO = {"FUV":1539, "NUV":2316, "u":3562, "g":4719, "r":6185, "i":7500, "z":8961}

COLOR_INFO = {"FUV":"xkcd:purple", "NUV":"xkcd:violet", "u":"xkcd:blue", "g":"xkcd:green", "r":"xkcd:red", "i":"xkcd:cyan", "z":"xkcd:magenta"}



class SED( BaseObject ):
    """
    LePhare SED fits
    """
    
    PROPERTIES         = ["sed_data"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
    def set_data(self, filename, nrows=1057):
        """
        
        """
        self._properties["sed_data"] = pd.read_table(filename, skiprows=20, names=["lbd", "mag"], sep="  ", engine="python", nrows=nrows)
        self.sed_data["flux"] = mag_to_flux(spec_data["mag"])
    
    def set_filter_bandpass(self, band, 
    
    
    
    
    # ================ #
    #  Properties      #
    # ================ #
    @property
    def sed_data(self):
        """  """
        return self._property["sed_data"]
    



class KCorrection( BaseObject ):
    """
    
    """
    
    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
    














