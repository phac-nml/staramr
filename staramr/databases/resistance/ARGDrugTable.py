import configparser
from collections import OrderedDict
from os import path

import pandas as pd

"""
Class which provides access to gene/drug mappings stored in tabular files. 
"""


class ARGDrugTable:
    DEFAULT_DATA_DIR = path.join(path.dirname(__file__), 'data')
    DEFAULT_INFO_FILE = path.join(DEFAULT_DATA_DIR, 'info.ini')
    DTYPES = {}

    def __init__(self, file=None, info_file=DEFAULT_INFO_FILE):
        """
        Creates a new ARGDrugTable with the given file, and info file (storing versions of each gene/drug table).
        :param file: The file containing the gene/drug mappings.
        :param info_file: The info file containing version information for the gene/drug mapping files.
        """
        self._info_file = info_file
        self._file = file

        if file is not None:
            self._data = pd.read_csv(file, sep='\t', dtype=self.DTYPES)

    def get_resistance_table_info(self):
        """
        Gets information about the antimcirobial resistance gene drug table versions.
        :return: A dictionary of the database gene drug table versions.
        """
        config = configparser.ConfigParser()
        config.read(self._info_file)
        return OrderedDict(config['Versions'])

    def _drug_string_to_correct_separators(self, drug):
        """
        Converts a drug string (separated by commas) to use correct separators/spacing.
        :param drug: The drug string.
        :return: The drug string with correct separators/spacing.
        """
        return ', '.join(drug.split(','))
