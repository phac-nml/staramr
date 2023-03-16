from os import path

import pandas as pd
import re

"""
A Class used to parse out a list of point mutation complexes that together confer resistance.
"""


class ComplexMutations:
    DEFAULT_COMPLEX_FILE = path.join(path.dirname(__file__), 'data', 'complex_mutations.tsv')

    def __init__(self, file=DEFAULT_COMPLEX_FILE):
        self._data = pd.read_csv(file, sep='\t')

    @classmethod
    def get_default_mutation_file(cls):
        """
        Get the default file containing the list of complex mutations.
        :return: The default file containing the list of complex mutations.
        """
        return cls.DEFAULT_COMPLEX_FILE
    
    def get_matches(self, mutation_codes):
        matches = []
        for row in self._data.itertuples():
            positions = re.split(', *', row.positions)
            mandatory = re.split(', *', row.mandatory)
            #print("start printing mandatory")
            print(mandatory)
            #print("stop printing mandatory")
            #print(row)
            if set(mandatory).issubset(set(mutation_codes)):
                print("YES MATCH")
                matches.append((set(positions).intersection(set(mutation_codes)), row.phenotype))
                print(matches)


