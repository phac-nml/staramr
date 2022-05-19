import logging
from os import path

from staramr.databases.resistance.ARGDrugTable import ARGDrugTable

logger = logging.getLogger("ARGDrugTablePointfinder")

"""
A Class used to load up and search a file containing gene/drug mappings for PointFinder results.
"""


class ARGDrugTablePointfinder(ARGDrugTable):
    DEFAULT_FILE = path.join(ARGDrugTable.DEFAULT_DATA_DIR, 'ARG_drug_key_pointfinder.tsv')
    DTYPES = {'Organism': str, 'Gene': str, 'Codon Pos.': int, 'Drug': str}

    def __init__(self, file=DEFAULT_FILE):
        """
        Builds a new ARGDrugTablePointfinder from the given file.
        :param file: The file containing the gene/drug mappings.
        """
        super().__init__(file=file)

    def get_drug(self, organism, gene, position):
        """
        Gets the drug given the organism, gene, and position of point mutation.
        :param organism: The organism.
        :param gene: The gene.
        :param position: The position of the point mutation (may be codon position or nucleotide position depending on gene).
        :return: The drug this mutation causes resistance to, or None if no such drug.
        """
        table = self._data

        drug = table[(table['Organism'] == organism) & (table['Gene'] == gene) & (
                table['Codon Pos.'] == position)]['Drug']
        if (drug.empty):
            logger.warning("No drug found for organism=%s, gene=%s, position=%s", organism, gene, position)
            return None
        else:
            return self._drug_string_to_correct_separators(drug.iloc[0])
