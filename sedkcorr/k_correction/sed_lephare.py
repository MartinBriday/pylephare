
import pandas
import os

from . import basesed




class SED_LePhare( basesed.SED ):
    """
    LePhare SED fits
    """
    
    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
    def __init__(self, data_sed=None, data_meas=None, z=None, **kwargs):
        """
        
        """
        if data_sed is not None and data_meas is not None and z is not None:
            self.set_data(data_sed, data_meas, z, **kwargs)
    
    def set_data(self, data_sed=None, data_meas=None, z=None, **kwargs):
        """
        
        """
        data_sed = self.get_sed_data(sed_data=data_sed, **kwargs)
        data_meas, z = self.get_meas_data(meas_data=data_meas, z=z, **kwargs)
    
        _ = super(SED_LePhare, self).set_data(data_sed, data_meas, z, **kwargs)
    
    
    def get_sed_data(self, sed_index=None, sed_data=None, sed_dir=None, nrows=1057, **kwargs):
        """
        Load the SED spectrum.
        
        Parameters
        ----------
        sed_index : [int or None]
            Object SED file index.
            A SED directory is needed to be usable (sed_dir)
            You can not use both sed_index and sed_data.
        
        sed_data : [pandas.DataFrame or None]
            Object SED data input in a DataFrame.
            You can not use both sed_data and sed_index.
        
        sed_dir : [string]
            Path to the SED files folder.
            Needed if sed_index used.
        
        
        Returns
        -------
        Void
        """
        if sed_index is not None and sed_data is None:
            if sed_dir is None:
                raise ValueError("You did not input the sprectrum directory.")
            sed_index = str(sed_index)
            while len(sed_index)<9:
                sed_index = "0"+sed_index
            sed_filename = "Id"+sed_index+".spec"
            return pandas.read_table(os.path.expanduser(sed_dir+sed_filename), skiprows=20, names=["lbda", "mag"], sep="  ", engine="python", nrows=nrows)
        elif sed_index is None and sed_data is not None:
            return sed_data
        else:
            raise ValueError("You must input one (and only one) data set option : the spectrum index or a DataFrame with columns = ['lbda', 'mag'].")
    
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
        return [band for band in basesed.LIST_BANDS if FILTER_BANDS[band]["context_id"] in idx]
    
    def get_meas_data(self, meas_data=None, z=None, col_syntax=["mag_band", "mag_band_err"], list_bands=None, **kwargs):
        """
        Set the host redshift and the measured magnitudes for every filter bands used in SED fitting.
        
        Parameters
        ----------
        meas_data : [table like]
            Table like (eg : DataFrame line) of the measurements.
        
        z : [float or None]
            Redshift of the SNeIa host.
            If None, the redshift is supposed to be in the meas_data table under the name "Z-SPEC".
            Default is None.
        
        col_syntax : [list[string]]
            Syntax of measurements and errors column names in the meas_data table.
            Replace the filter band in the column names with the word "band" (eg: ["mag_band", "mag_band_err"]).
        
        Options
        -------
        list_bands : [list[string] or None]
            List of the filter bands used in SED fitting.
            If None, the LePhare context set list_bands. The context is supposed to be in the meas_data table under the name "CONTEXT".
        
        
        Returns
        -------
        Void
        """
        self._side_properties["list_bands"] = self.context_filters(meas_data["CONTEXT"]) if list_bands is None else list_bands
        z = z if z is not None else meas_data["Z-SPEC"]
        
        data_meas = {band:{"mag":meas_data[col_syntax[0].replace("band",band)],
                           "mag.err":meas_data[col_syntax[1].replace("band",band)]}
                     for band in self.list_bands}

        return data_meas, z
    
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

    
    


