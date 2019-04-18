
import numpy as np

from prospect.io.read_results import results_from, get_sps
from prospect.io.read_results import traceplot, subcorner


from ..sed_fitting import prospector
from . import basesed




class SED_prospector( basesed.SED ):
    """
    
    """
    
    PROPERTIES         = ["p_res", "p_obs", "p_mod"]
    SIDE_PROPERTIES    = ["p_run_params", "p_sps"]
    DERIVED_PROPERTIES = []
    
    def __init__(self, filename=None, **kwargs):
        """
        
        """
        if filename is not None:
            self.read_fit_results(filename)
            self.set_data(**kwargs)

    def read_fit_results(self, filename=None):
        """
        
        """
        res, obs, mod = results_from(filename, dangerous=False)
        self._properties["p_res"] = res
        self._properties["p_obs"] = obs
        self._properties["p_mod"] = prospector.load_model(**self.p_run_params)
    
    def set_data_sed(self, **extras):
        """
        
        """
        imax = np.argmax(self.p_res["lnprobability"])
        theta_max = self.p_res["chain"][imax, :].copy()
        data_sed = {"lbda":self.get_sed_wavelength,
                    "flux":self.get_sed_flux(theta_max)}
        _ = super(SED_prospector, self).set_data_sed( pandas.DataFrame(data_sed))
    
    def context_filters(self, context):
        """
        Return a list of the concerned filter bands relative to the given context.
        
        Parameters
        ----------
        context : [int]
        LePhare type context, it defines the used filter bands for the SED fitting.
        
        
        Returns
        -------
        list(string)
        """
        idx = []
        for ii in range(len(basesed.LIST_BANDS)-1,-1,-1):
            if (context - 2**ii) >= 0:
                context = context - 2**ii
                idx.append(ii)
        return [band for band in basesed.LIST_BANDS if basesed.FILTER_BANDS[band]["context_id"] in idx]
    
    def get_data_meas(self, data_meas=None, z=None, col_syntax=["mag_band", "mag_band_err"], list_bands=None, **extras):
        """
        Set the host redshift and the measured magnitudes for every filter bands used in SED fitting.
        
        Parameters
        ----------
        data_meas : [table like]
        Table like (eg : DataFrame line) of the measurements.
        
        z : [float or None]
        Redshift of the SNeIa host.
        If None, the redshift is supposed to be in the data_meas table under the name "Z-SPEC".
        Default is None.
        
        col_syntax : [list[string]]
        Syntax of measurements and errors column names in the data_meas table.
        Replace the filter band in the column names with the word "band" (eg: ["mag_band", "mag_band_err"]).
        
        Options
        -------
        list_bands : [list[string] or None]
        List of the filter bands used in SED fitting.
        If None, the LePhare context set list_bands. The context is supposed to be in the data_meas table under the name "CONTEXT".
        
        
        Returns
        -------
        Void
        """
        
        ############# Add an option to read phot data from prospector results ###############
        
        self._side_properties["list_bands"] = self.context_filters(data_meas["CONTEXT"]) if list_bands is None else list_bands
        z = z if z is not None else data_meas["Z-SPEC"]
        
        data_meas = {band:{"mag":data_meas[col_syntax[0].replace("band",band)],
                           "mag.err":data_meas[col_syntax[1].replace("band",band)]}
                     for band in self.list_bands}
        
        return data_meas, z
    
    def get_sed_flux(self, theta):
        """
        
        """
        mspec, mphot, mextra = self.p_mod.mean_model(theta, self.p_obs, sps=self.p_sps)
        return mspec
    
    def get_sed_wavelength(self):
        """
        
        """
        if self.p_obs["wavelength"] is None:
            # *restframe* spectral wavelengths, since obs["wavelength"] is None
            a = 1.0 + self.p_run_params.get('zred', 0.0)
            wspec = self.p_sps.wavelengths.copy()
            wspec *= a #redshift them
        else:
            wspec = self.p_obs["wavelength"]
        return wspec
    
    def k_correction(self):
        """
        Recover the integrated flux from every filter bands from the shifted SED.
        Then convert them into magnitudes.
        
        
        Returns
        -------
        Void
        """
        _ = super(SED_LePhare, self).k_correction()
        for band in self.list_bands:
            self.data_kcorr[band]["mag.err"] = self.data_meas[band]["mag.err"]
            self.data_kcorr[band]["flux.err"] = self.data_meas[band]["flux.err"]
        
        for band in basesed.FILTER_BANDS:
            if band not in self.list_bands:
                self.data_kcorr[band]["flux.err"] = 0.
                self.data_kcorr[band]["mag.err"] = 0.

    
    





    #-------------------#
    #   Properties      #
    #-------------------#
    @property
    def p_res(self):
        """  """
        return self._properties["p_res"]

    @property
    def p_obs(self):
        """  """
        return self._properties["p_obs"]

    @property
    def p_mod(self):
        """  """
        return self._properties["p_mod"]

    @property
    def p_run_params(self):
        """  """
        if self._side_properties["p_run_params"] is None:
            self._side_properties["p_run_params"] = self.prospector_res["run_params"]
        return self._side_properties["p_run_params"]
    
    @property
    def p_sps(self):
        """  """
        if self._side_properties["p_sps"] is None:
            self._side_properties["p_sps"] = prospector.load_sps(**self.p_run_params)
        return self._side_properties["p_sps"]
