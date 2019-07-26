
import numpy as np
import pandas
import pkg_resources

from propobject import BaseObject


class LePhareSEDFitter( BaseObject ):
    """
    
    """

    PROPERTIES         = []
    SIDE_PROPERTIES    = []
    DERIVED_PROPERTIES = []

    def __init__(self, **kwargs):
        """
        
        """
        self.set_data(**kwargs)

    def set_data(self, data_path=None, ):
        """
        
        """
        return
    
    def _write_para_(self, line, new_param_value):
        """
        
        """
        if type(line) == int:
            idx_line = line
        elif type(line) == str:
            with open(pkg_resources.resource_filename(__name__, "bin/lephare_zphot_input.para"), "r") as file:
                file_buf = [line for line in file]
            for ii, line in enumerate(file_buf):
                if line.split(" ")[0] == line:
                    idx_line = ii
                    break
        else:
            raise ValueError("'line' must be either integer (line index) or string (first word of the line)")
        
        line_splited = file_buf[idx_line].split()
        line_splited[1] = new_param_value
        file_buf[idx_line] = " ".join(line_splited) + "\n"

        with open(pkg_resources.resource_filename(__name__, "bin/lephare_zphot_input.para"), "w") as file:
            for line in file_buf:
                file.write(line)
    
    def _set_data_path_(self, data_path):
        """
        
        """
        self._write_para_("CAT_IN", data_path)

    #-------------------#
    #   Properties      #
    #-------------------#
    @property



