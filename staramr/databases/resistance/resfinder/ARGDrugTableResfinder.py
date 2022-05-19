import logging
from os import path

from staramr.databases.resistance.ARGDrugTable import ARGDrugTable

logger = logging.getLogger("ARGDrugTableResfinder")

"""
A Class used to load up and search a file containing gene/drug mappings for ResFinder results.
"""


class ARGDrugTableResfinder(ARGDrugTable):
    DEFAULT_FILE = path.join(ARGDrugTable.DEFAULT_DATA_DIR, 'ARG_drug_key_resfinder.tsv')
    DTYPES = {'Class': str, 'Gene': str, 'Accession': str, 'Drug': str}

    def __init__(self, file=DEFAULT_FILE):
        """
        Builds a new ARGDrugTableResfinder from the given file.
        :param file: The file containing the gene/drug mappings.
        """
        super().__init__(file=file)

    def get_drug(self, drug_class, gene_plus_variant, accession):
        """
        Gets the drug given the drug class, gene (plus variant of gene encoded in ResFinder database) and accession.
        :param drug_class: The drug class.
        :param gene_plus_variant: The gene plus variant (e.g., {gene}_{variant} = {blaIMP-58}_{1}).
        :param accession: The accession in the resfinder database (e.g., KU647281).
        :return: The particular drug, or None if no matching drug was found.
        """
        table = self._data

        drug = table[(table['Class'] == drug_class) & (table['Gene'] == gene_plus_variant) & (
                table['Accession'] == accession)]['Drug']
        if (drug.empty):
            logger.warning("No drug found for drug_class=%s, gene=%s, accession=%s", drug_class, gene_plus_variant,
                           accession)
            return None
        else:
            return self._drug_string_to_correct_separators(drug.iloc[0])
