import logging
from os import path

from staramr.databases.resistance.ARGDrugTable import ARGDrugTable

logger = logging.getLogger("ARGDrugTableResfinder")


class ARGDrugTableResfinder(ARGDrugTable):
    DEFAULT_FILE = path.join(ARGDrugTable.DEFAULT_DATA_DIR, 'ARG_drug_key_resfinder.tsv')

    def __init__(self, file=DEFAULT_FILE):
        super().__init__(file=file)

    def get_drug(self, drug_class, gene_plus_variant, accession):
        table = self._data

        drug = table[(table['Class'] == drug_class) & (table['Gene'] == gene_plus_variant) & (
                table['Accession'] == accession)]['Drug']
        if (drug.empty):
            logger.warning(
                "No drug found for drug_class=" + drug_class + ", gene=" + gene_plus_variant + ", accession=" + accession)
            return None
        else:
            return drug.iloc[0]
