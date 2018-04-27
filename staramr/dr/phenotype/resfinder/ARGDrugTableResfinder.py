import logging

from staramr.dr.phenotype.ARGDrugTable import ARGDrugTable

logger = logging.getLogger("ARGDrugTableResfinder")


class ARGDrugTableResfinder(ARGDrugTable):

    def __init__(self, file):
        super().__init__(file)

    def get_drug(self, drug_class, gene_plus_variant, accession):
        table = self.data

        drug = table[(table['Class'] == drug_class) & (table['Gene'] == gene_plus_variant) & (
                table['Accession'] == accession)]['Drug']
        if (drug.empty):
            logger.warning(
                "No drug found for drug_class=" + drug_class + ", gene=" + gene_plus_variant + ", accession=" + accession)
            return None
        else:
            return drug.iloc[0]
