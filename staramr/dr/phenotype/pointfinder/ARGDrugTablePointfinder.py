import logging

from staramr.dr.phenotype.ARGDrugTable import ARGDrugTable

logger = logging.getLogger("ARGDrugTablePointfinder")


class ARGDrugTablePointfinder(ARGDrugTable):

    def __init__(self, file):
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
