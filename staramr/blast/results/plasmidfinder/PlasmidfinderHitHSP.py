import logging
import re
from typing import List

import pandas as pd

from staramr.blast.results.AMRHitHSP import AMRHitHSP

logger = logging.getLogger('PlasmidfinderHitHSP')

"""
A Class storing a PlasmidFinder-specific BLAST hit/HSP.
"""


class PlasmidfinderHitHSP(AMRHitHSP):

    def __init__(self, file: str, blast_record: pd.Series) -> None:
        """
        Builds a new PlasmidfinderHitHSP.
        :param file: The input file.
        :param blast_record: The Bio.Blast.Record this hit came from.
        :param hit: The particular Bio.Blast.Record.Alignment.
        :param hsp: The particular Bio.Blast.Record.HSP.
        """
        super().__init__(file, blast_record)

        logger.debug("record=%s", self._blast_record)

        splitList = re.split('_', self.get_amr_gene_id())
        re_search = list(filter(None, splitList))  # type: List[str]

        if not re_search:
            raise Exception("Could not split up seq name for [" + self.get_amr_gene_id() + "]")

        length = len(re_search)

        # Add empty string if gene variant is missing
        if length == 2:
            re_search.insert(1, '')
            length += 1

        if 3 <= length <= 5:
            self._gene = re_search[0]
            self._gene_variant = re_search[1]

            if length == 3:
                self._accession = re_search[2]
            elif length == 4:
                self._accession = re_search[3]
            elif length == 5:
                self._accession = re_search[3] + "_" + re_search[4]
            else:
                raise Exception("Length not valid for variable assignment, for [" + self.get_amr_gene_id() + "]")
        else:
            raise Exception("Index length not within range, for [" + self.get_amr_gene_id() + "]")

    def get_amr_gene_name(self) -> str:
        """
        Gets the gene name for the PlasmidFinder hit.
        :return: The gene name.
        """
        return self._gene

    def get_amr_gene_name_with_variant(self) -> str:
        """
        Gets the gene name + variant number for the PlasmidFinder hit.
        :return: The gene name + variant number.
        """
        return self.get_amr_gene_name() + '_' + self._gene_variant

    def get_amr_gene_variant_accession(self) -> str:
        """
        Gets the gene name + variant number + accession for the PlasmidFinder hit.
        :return: The gene name + variant number + accession.
        """
        return self._gene + '_' + self._gene_variant + '_' + self._accession

    def get_amr_gene_variant(self) -> str:
        """
        Gets the variant number for the PlasmidFinder hit.
        :return: The variant number.
        """
        return self._gene_variant

    def get_amr_gene_accession(self) -> str:
        """
        Gets the accession for the PlasmidFinder hit.
        :return: The accession.
        """
        return self._accession
