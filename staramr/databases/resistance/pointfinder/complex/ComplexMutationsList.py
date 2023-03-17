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
    
    def get_matches(self, results_table, hit):

        matches = []

        # We can't use "Pointfinder Position", because the gene is missing from it, and we can't use
        # "Gene", because those co-ords may contain a position correction for indels, meaning they
        # won't match the position expected in the database files.
        mutation_codes = list(results_table["Isolate ID"] + "(" + results_table["Pointfinder Position"] + ")")

        for row in self._data.itertuples():
            positions = re.split(', *', row.positions)
            mandatory = re.split(', *', row.mandatory)

            if set(mandatory).issubset(set(mutation_codes)):
                intersection = list(set(positions).intersection(set(mutation_codes)))
                mutation_positions = [re.sub("[^0-9]", "", i) for i in list(intersection)]

                result = [hit.get_genome_id(),
                        ", ".join(list(intersection)), # Technically these are also Pointfinder co-ords
                        row.phenotype,
                        "complex",
                        ", ".join(mutation_positions),
                        "complex", # Creating a mutation string would be confusing for this.
                        hit.get_pid(),
                        hit.get_plength(),
                        str(hit.get_hsp_length()) + "/" + str(hit.get_amr_gene_length()),
                        hit.get_genome_contig_id(),
                        hit.get_genome_contig_start(),
                        hit.get_genome_contig_end(),
                        ", ".join(list(intersection))
                        ]

                matches.append(result)
        
        return matches
