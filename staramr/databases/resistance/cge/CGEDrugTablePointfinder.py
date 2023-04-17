import logging
import os
from pathlib import Path

from staramr.databases.resistance.ARGDrugTable import ARGDrugTable

logger = logging.getLogger("CGEDrugTablePointfinder")

"""
A class used to load up and search a file containing gene/drug mappings for CGE Pointfinder results.
"""


class CGEDrugTablePointfinder(ARGDrugTable):
    DTYPES = {'#Gene_accession': str, 'TypeGene': str, 'Mutation ID': str, 'Codon_pos': str,
              'Ref_nuc': str, "Ref_codon": str, "Res_codon": str,
              'Class': str, 'Phenotype': str, 'PMID': str,
              'Mechanism of resistance': str, 'Notes': str, 'Required mutation': str,
              '% abundance mutation required': str}

    def __init__(self, file):
        """
        Builds a new CGEDrugTablePointfinder from the given file.
        :param file: The file containing the gene/drug mappings.
        """
        super().__init__(file=file)

    def get_drug(self, gene, position):
        """
        Gets the drug given the organism, gene, and position of point mutation.
        :param gene: The gene.
        :param position: The position of the point mutation (may be codon position or nucleotide position depending on gene).
        :return: The drug this mutation causes resistance to, or None if no such drug.
        """
        table = self._data
        print("table")
        print(table)

        # The "Mutation ID" column looks like: GLY-15
        # Where the first 3 letters are the gene in upper case, followed by a dash,
        # followed by the position of the mutation.
        mutation_id = gene.upper() + "-" + str(position)

        drug = table[(table['Mutation ID'] == mutation_id)]['Phenotype']

        if (drug.empty):
            logger.warning("No drug found for gene=%s, position=%s", gene, position)
            return None
        else:
            return drug.iloc[0]