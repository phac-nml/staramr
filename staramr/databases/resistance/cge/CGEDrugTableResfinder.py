import logging
import os
from pathlib import Path

from staramr.databases.resistance.ARGDrugTable import ARGDrugTable
from staramr.results.AMRDetectionSummaryResistance import AMRDetectionSummaryResistance

logger = logging.getLogger("CGEDrugTableResfinder")

"""
A Class used to load up and search a file containing gene/drug mappings for CGE ResFinder results.
"""


class CGEDrugTableResfinder(ARGDrugTable):
    DTYPES = {'Gene_accession no.': str, 'Class': str, 'Phenotype': str, 'PMID': str,
              'Mechanism of resistance': str, "Notes": str, "Required_gene": str}

    def __init__(self, file):
        """
        Builds a new CGEDrugTableResfinder from the given file.
        :param file: The file containing the gene/drug mappings.
        """
        super().__init__(file=file)

        self._data['Class'] = self._data['Class'].str.lower()

    def get_drug(self, drug_class, gene_plus_variant, accession):
        """
        Gets the drug given the drug class, gene (plus variant of gene encoded in ResFinder database) and accession.
        :param drug_class: The drug class.
        :param gene_plus_variant: The gene plus variant (e.g., {gene}_{variant} = {blaIMP-58}_{1}).
        :param accession: The accession in the resfinder database (e.g., KU647281).
        :return: The particular drug, or None if no matching drug was found.
        """
        table = self._data

        gene_accession = str(gene_plus_variant) + "_" + str(accession)
        drug = table[(table['Gene_accession no.'] == gene_accession)]['Phenotype']

        if (drug.empty):
            logger.warning("No drug found for drug_class=%s, gene=%s, accession=%s", drug_class, gene_plus_variant,
                           accession)
            return None
        elif len(drug) > 1:
            logger.warning("Multiple entries found for drug_class=%s, gene=%s, accession=%s", drug_class,
                           gene_plus_variant, accession)
            return AMRDetectionSummaryResistance.aggregate_phenotype(drug)
        else:
            return drug.iloc[0]
        
    def get_cge_notes(self, gene_plus_variant, accession):
        """
        Gets the notes for a gene (plus variant of gene encoded in ResFinder database) and accession.
        :param gene_plus_variant: The gene plus variant (e.g., {gene}_{variant} = {blaIMP-58}_{1}).
        :param accession: The accession in the resfinder database (e.g., KU647281).
        :return: A string containtain the Notes, if they exist, or the empty string ("") if
                 there are no Notes.
        """
        table = self._data

        gene_accession = str(gene_plus_variant) + "_" + str(accession)
        matches = table[(table['Gene_accession no.'] == gene_accession)]
        matches = matches.fillna("")

        # There's a chance of having multiple matches:
        notes = ';'.join(matches["Notes"])

        return notes
