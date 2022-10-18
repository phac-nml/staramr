from os import path

import pandas as pd

"""
A Class used to parse out a list of genes to exclude from the results.
"""


class ExcludeGenesList:
    DEFAULT_EXCLUDE_FILE = path.join(path.dirname(__file__), 'data', 'genes_to_exclude.tsv')

    def __init__(self, file=DEFAULT_EXCLUDE_FILE):
        self._data = pd.read_csv(file, sep='\t')

    def tolist(self):
        """
        Converts the exclude genes data to a list.
        :return: A list with genes to exclude.
        """
        return self._data['gene_id'].tolist()

    @classmethod
    def get_default_exclude_file(cls):
        """
        Get the default file containing the list of genes to exclude.
        :return: The default file containing the list of genes to exclude.
        """
        return cls.DEFAULT_EXCLUDE_FILE
