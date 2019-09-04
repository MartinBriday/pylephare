
import numpy as np
import pandas
import os

from . import basesed


def lephare_output_file_reader(filename, filter_list):
    """
    Read a LePhare output file, returning a pandas.DataFrame.
    
    Parameters
    ----------
    filename : [string]
        Path of the LePhare file.
    
    filter_list : [list(string)]
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


class SED_LePhare( basesed.SED ):
    """
    This class read LePhare SED fits and apply k-correction on it.
    Child of the class 'basesed.SED'.
    """
    
    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []
    
    def set_data_sed(self, data_sed=None, sed_index=None, sed_dir=None, nrows=1057, **kwargs):
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
        
        nrows : [int]
            Number of rows to account for the spectrum.
            Generally, LePhare output three series of data (GAL, QSO, STAR) in a row.
            This number allow you to only keep GAL serie.
        
        
        Returns
        -------
        Void
        """
        if sed_index is not None and data_sed is None:
            if sed_dir is None:
                raise ValueError("You did not input the sprectrum directory.")
            sed_index = str(sed_index)
            while len(sed_index)<9:
                sed_index = "0"+sed_index
            sed_filename = "Id"+sed_index+".spec"
            idx_start = 10
            data_sed =  pandas.read_csv(os.path.expanduser(sed_dir+sed_filename),
                                        skiprows=idx_start, sep="  ", engine="python", nrows=nrows)
            while data_sed.shape[1] != 2:
                idx_start += 1
                data_sed =  pandas.read_csv(os.path.expanduser(sed_dir+sed_filename),
                                            skiprows=idx_start, sep="  ", engine="python", nrows=nrows)
            data_sed =  pandas.read_csv(os.path.expanduser(sed_dir+sed_filename),
                                        skiprows=idx_start, names=["lbda", "mag"], sep="  ",
                                        engine="python", nrows=nrows)
        elif sed_index is None and data_sed is not None:
            data_sed = data_sed
        else:
            raise ValueError("You must input one (and only one) data set option : the spectrum index or a DataFrame with columns = ['lbda', 'mag'].")
        
        data_sed["mag.err"] = np.zeros(len(data_sed))
        _ = super(SED_LePhare, self).set_data_sed( data_sed )
    
    def context_filters(self, context):
        """
        Return a list of the concerned filter bands relative to the given context.
        
        Parameters
        ----------
        context : [int]
            LePhare type context, it defines the used filter bands for the SED fitting.
            Knowing the bands in 'self.list_bands', the number we want to set is : sum(2**[band_nb]).
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
        idx = []
        for ii in range(len(basesed.LIST_BANDS)-1,-1,-1):
            if (context - 2**ii) >= 0:
                context = context - 2**ii
                idx.append(ii)
        return [band for band in basesed.LIST_BANDS if basesed.FILTER_BANDS[band]["context_id"] in idx]
    
    def set_data_meas(self, data_meas=None, z=None, col_syntax=["mag_band", "mag_band_err"], list_bands=None, **kwargs):
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
        self._side_properties["list_bands"] = self.context_filters(int(data_meas["CONTEXT"])) if list_bands is None else list_bands
        z = z if z is not None else data_meas["Z-SPEC"]
        
        data_meas = {band:{"mag":data_meas[col_syntax[0].replace("band",band)],
                           "mag.err":data_meas[col_syntax[1].replace("band",band)]}
                     for band in self.list_bands}
        _ = super(SED_LePhare, self).set_data_meas(data_meas=data_meas, z=z)
    
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

    
    


