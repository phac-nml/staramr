from os import path

import pandas as pd
import re

"""
A class used to parse out a list of point mutation complexes that together confer resistance.
"""


class ComplexMutations:
    DEFAULT_COMPLEX_FILE = path.join(path.dirname(__file__), 'data', 'complex_mutations.tsv')

    def __init__(self, file=DEFAULT_COMPLEX_FILE):
        self._data = pd.read_csv(file, sep='\t')

        # Convert strings to integer lists for positions:
        self._data["positions"] = self._data["positions"].apply(lambda str: re.split(', *', str))
        self._data["mandatory"] = self._data["mandatory"].apply(lambda str: re.split(', *', str))

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
        mutation_codes = list(str(hit.get_amr_gene_name()) + " (" + results_table["Pointfinder Position"] + ")")

        for row in self._data.itertuples():

            if set(row.mandatory).issubset(set(mutation_codes)):
                intersection = list(set(row.positions).intersection(set(mutation_codes)))

                # Sorting, because it makes output reproducable.
                # .sort() is in-place.
                intersection.sort()

                # Similar operation as above for positions, except we need to sort them
                # as integers, rather than strings:      
                mutation_positions = [re.sub("[^0-9]", "", i) for i in intersection]
                mutation_positions = list(map(int, mutation_positions))
                mutation_positions.sort()  # in-place
                mutation_positions = list(map(str, mutation_positions))

                result = [hit.get_genome_id(),
                        ", ".join(intersection), # Technically these are also Pointfinder co-ords
                        row.phenotype,
                        pd.NA,  # CGE-predicted phenotype
                        "complex",  # Type
                        ", ".join(mutation_positions),
                        "complex", # Creating a mutation string would be confusing for this.
                        hit.get_pid(),
                        hit.get_plength(),
                        str(hit.get_hsp_length()) + "/" + str(hit.get_amr_gene_length()),
                        hit.get_genome_contig_id(),
                        hit.get_genome_contig_start(),
                        hit.get_genome_contig_end(),
                        ", ".join(intersection),
                        pd.NA, # The CGE notes.
                        pd.NA,
                        pd.NA,
                        pd.NA]

                matches.append(result)

        return matches
