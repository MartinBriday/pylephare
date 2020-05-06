

from .io import get_filter_bandpass

class _FilterHolder_( object ) :
    """ Handle filters (attributes, bandpasses, labels, numbers, ...) """
    
    def set_filters(self, filters):
        """
        Set the filters as attributes and reset the bandpasses.
        
        Parameters
        ----------
        filters : [list(string)]
            List of filters used for the SED fitting.
            Must be known filters by LePhare with the format instrument.band, for instance sdss.u or ps1.z.
        
        
        Returns
        -------
        Void
        """
        self._filters = filters
        self._filter_bandpass = None

    # ============== #
    #  Properties    #
    # ============== #
    @property
    def filters(self):
        """ List of filters available in the input data catalog """
        if not hasattr(self,"_filters"):
            self._filters = None
        return self._filters

    def has_filters(self):
        """ Test that filters are loaded """
        return self.filters is not None and len(self.filters)>0

    @property
    def _filter_labels(self):
        """ Filter label, used for file making """
        return "_".join([l.replace(".","") for l in self.filters])

    @property
    def filter_bandpasses(self):
        """ Dictionary containing the bandpasses for each used filter """
        if not hasattr(self, "_filter_bandpass") or self._filter_bandpass is None:
            if not self.has_filters():
                raise AttributeError("No filter loaded")
            self._filter_bandpass = {filter_: get_filter_bandpass(filter_) for filter_ in self.filters}
            
        return self._filter_bandpass
            
    @property
    def nfilters(self):
        """ Number of filters """
        if not self.has_filters():
            raise AttributeError("No filter set yet.")
        return len(self.filters)
