from os import path

import pandas


class ARGDrugTable:
    DEFAULT_DATA_DIR = path.join(path.dirname(__file__), 'data')
    DEFAULT_INFO_FILE = path.join(DEFAULT_DATA_DIR, 'info.ini')

    def __init__(self, file=None, info_file=DEFAULT_INFO_FILE):
        self._info_file = info_file
        self._file = file

        if file is not None:
            self._data = pandas.read_csv(file, sep="\t")

    def get_resistance_table_info(self):
        """
        Gets information about the antimcirobial resistance gene drug table versions.
        :return: A list of key/value for the ResFinder and PointFinder versions.
        """
        database_info = pandas.read_csv(self._info_file, sep="=", index_col=False, header=None)
        return database_info.as_matrix().tolist()
