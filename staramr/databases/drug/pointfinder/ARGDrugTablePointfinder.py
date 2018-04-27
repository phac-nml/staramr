import logging
from os import path

from staramr.databases.drug.ARGDrugTable import ARGDrugTable

logger = logging.getLogger("ARGDrugTablePointfinder")


class ARGDrugTablePointfinder(ARGDrugTable):

    DEFAULT_FILE = path.join(ARGDrugTable.DEFAULT_DATA_DIR, 'ARG_drug_key_pointfinder.tsv')

    def __init__(self, file = DEFAULT_FILE):
        super().__init__(file)

    def get_drug(self, organism, gene, codon_position):
        table = self.data

        drug = table[(table['Organism'] == organism) & (table['Gene'] == gene) & (
                table['Codon Pos.'] == codon_position)]['Drug']
        if (drug.empty):
            logger.warning(
                "No drug found for organism=" + organism + ", gene=" + gene + ", codon_position=" + codon_position)
            return None
        else:
            return drug.iloc[0]
