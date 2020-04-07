

from .io import get_filter_bandpass

class _FilterHolder_( object ) :
    
    def set_filters(self, filters):
        """ """
        self._filters = filters
        self._filter_bandpass = None

    # ============== #
    #  Properties    #
    # ============== #
    @property
    def filters(self):
        """ """
        if not hasattr(self,"_filters"):
            self._filters = None
        return self._filters

    def has_filters(self):
        """ """
        return self.filters is not None and len(self.filters)>0

    @property
    def _filter_labels(self):
        """ """
        return "_".join([l.replace(".","") for l in self.filters])

    @property
    def filter_bandpasses(self):
        """ """
        if not hasattr(self, "_filter_bandpass") or self._filter_bandpass is None:
            if not self.has_filters():
                raise AttributeError("No filter loaded")
            self._filter_bandpass = {filter_: get_filter_bandpass(filter_) for filter_ in self.filters}
            
        return self._filter_bandpass
            
    @property
    def nfilters(self):
        """ """
        if not self.has_filters():
            raise AttributeError("No filter set yet.")
        return len(self.filters)
