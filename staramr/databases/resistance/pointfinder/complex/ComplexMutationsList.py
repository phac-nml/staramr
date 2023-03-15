from os import path

import pandas as pd

"""
A Class used to parse out a list of point mutation complexes that together confer resistance.
"""


class ComplexMutationsList:
    DEFAULT_COMPLEX_FILE = path.join(path.dirname(__file__), 'data', 'complex_phenotypes.tsv')

    def __init__(self, file=DEFAULT_COMPLEX_FILE):
        self._data = pd.read_csv(file, sep='\t')

    @classmethod
    def get_default_complex_file(cls):
        """
        Get the default file containing the list of complexes.
        :return: The default file containing the list of complexes.
        """
        return cls.DEFAULT_COMPLEX_FILE
