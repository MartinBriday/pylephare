
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






class KCorrection( BaseObject ):
    """
    LePhare SED fits
    """
    
    PROPERTIES         = ["sed_data", "filter_bandpass", "meas_data", "z"]
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = ["sed_shifted", "kcorr_data"]
    
    FILTER_BANDS = ["FUV", "NUV", "u", "g", "r", "i", "z"]

    LBDA_INFO = {"FUV":1539, "NUV":2316, "u":3562, "g":4719, "r":6185, "i":7500, "z":8961}

    COLOR_INFO = {"FUV":"xkcd:purple", "NUV":"xkcd:violet", "u":"xkcd:blue", "g":"xkcd:green", "r":"xkcd:red", "i":"xkcd:cyan", "z":"xkcd:magenta"}
    
    def set_sed_data(self, sed_filename=None, sed_data=None, nrows=1057):
        """
        
        """
        if sed_filename is not None and sed_data is None:
            self._properties["sed_data"] = pd.read_table(filename, skiprows=20, names=["lbda", "mag"], sep="  ", engine="python", nrows=nrows)
        if sed_filename is None and sed_data is not None:
            self._properties["sed_data"] = sed_pandas
        else:
            raise ValueError("You must input one (and only one) data set option : the filename directory or a DataFrame with columns = ['lbda', 'mag'].")
        self.sed_data["flux"] = mag_to_flux(self.sed_data["mag"])
    
    def set_filter_bandpass(self, opt_bands=None, from_sncosmo=False):
        """
        
        """
        opt_bands = self.FILTER_BANDS if opt_bands is None else opt_bands if type(opt_bands)==list else [opt_bands]
        list_sdss_bands = ["u", "g", "r", "i", "z"]
        if from_sncosmo:
            self.filter_bandpass = {band:pd.DataFrame({"lbda":bandpasses.get_bandpass("sdss"+band).wave, 
                                                       "trans":bandpasses.get_bandpass("sdss"+band).trans}) 
                                    for band in opt_bands if band in list_sdss_bands}
        else:
            self.filter_bandpass = {band:pd.read_table('filter_bandpass/SLOAN_SDSS.'+band+'.dat', sep=' ', names=["lbda", "trans"])
                                    for band in opt_bands if band in list_sdss_bands}
        if "FUV" in opt_bands:
            self.filter_bandpass["FUV"] = pd.read_table('filter_bandpass/GALEX_GALEX.FUV.dat', sep=' ', names=["lbda", "trans"])
        if "NUV" in opt_bands:
            self.filter_bandpass["NUV"] = pd.read_table('filter_bandpass/GALEX_GALEX.NUV.dat', sep=' ', names=["lbda", "trans"])
    
    def set_meas_data(self, meas_data=None, z=None, col_syntax=["mag_band", "mag_band_err"]):
        """
        
        """
        z = z if z is not None else self.meas_data["Z-SPEC"]
        col_syntax = ["band", "band.err"] if col_syntax is None else col_syntax
        
        self._properties["meas_data"] = {band:{"mag":meas_data[col_syntax[0].replace("band",band)], 
                                               "mag.err":meas_data[col_syntax[1].replace("band",band)], 
                                               "flux":mag_to_flux(meas_data[col_syntax[0].replace("band",band)]), 
                                               "flux.err":mag_to_flux_err(meas_data[col_syntax[0].replace("band",band)], 
                                                                          meas_data[col_syntax[1].replace("band",band)])}
                                         for band in self.FILTER_BANDS}
        
        self._properties["z"] = z
    
    def shift_sed(self):
        """
        
        """
        self._derived_properties["sed_shifted"] = pd.DataFrame({"lbda":lbda_z0(self.sed_data["lbda"], self.z), 
                                                                "flux":flux_z0(self.sed_data["flux"], self.z)})
    
    def k_correction(self):
        """
        
        """
        self._derived_properties["kcorr_data"] = {band:{"mag.err":self.meas_data[band]["mag.err"], 
                                                        "flux":spectroscopy.synthesize_photometry(self.sed_shifted["lbda"], self.sed_shifted["flux"], 
                                                                                                  self.filter_bandpass[band]["lbda"], self.filter_bandpass[band]["trans"]), 
                                                        "flux.err":self.meas_data[band]["flux.err"]}
                                                  for band in self.FILTER_BANDS}
        for band in self.FILTER_BANDS:
            self.kcorr_data[band]["mag"] = -9999999. #TO BE DONE
    
    def show(self, ):
        """
        
        """
        
    
    
    
    
    
    
    # ================ #
    #  Properties      #
    # ================ #
    @property
    def sed_data(self):
        """  """
        return self._property["sed_data"]
    
    @property
    def filter_bandpass(self):
        """  """
        if self._properties["filter_bandpass"] is None:
            self._properties["filter_bandpass"] = {}
        return self._properties["filter_bandpass"]
    
    @property
    def meas_data(self):
        """  """
        if self._properties["meas_data"] is None:
            self._properties["meas_data"] = {}
        return self._properties["meas_data"]
    
    @property
    def z(self):
        """  """
        return self._properties["z"]
    
    @property
    def sed_shifted(self):
        """  """
        if self._derived_properties["sed_shifted"] is None:
            self.shift_sed()
        return self._derived_properties["sed_shifted"]
    
    @property
    def kcorr_data(self):
        """  """
        if self._derived_properties["kcorr_data"] is None:
            self.k_correction()
        return self._derived_properties["kcorr_data"]
    
    














