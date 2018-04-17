import logging
import re

from staramr.blast.results.AMRHitHSP import AMRHitHSP

logger = logging.getLogger('ResfinderHitHSP')

"""
A Class storing a ResFinder-specific BLAST hit/HSP.
"""


class ResfinderHitHSP(AMRHitHSP):

    def __init__(self, file, blast_record, hit, hsp):
        """
        Builds a new ResfinderHitHSP.
        :param file: The input file.
        :param blast_record: The Bio.Blast.Record this hit came from.
        :param hit: The particular Bio.Blast.Record.Alignment.
        :param hsp: The particular Bio.Blast.Record.HSP.
        """
        super().__init__(file, blast_record, hit, hsp)

        re_search = re.search(r'([^_]+)_([^_]+)_([^_\s]+)$', hit.hit_id)
        if not re_search:
            raise Exception("Could not split up seq name for [" + hit.hit_id + "]")
        self._gene = re_search.group(1)
        self._gene_variant = re_search.group(2)
        self._accession = re_search.group(3)

        logger.debug("hit_id=" + str(hit.hit_id))

    def get_gene(self):
        """
        Gets the gene name for the ResFinder hit.
        :return: The gene name.
        """
        return self._gene

    def get_gene_with_variant(self):
        """
        Gets the gene name + variant number for the ResFinder hit.
        :return: The gene name + variant number.
        """
        return self.get_gene() + '_' + self._gene_variant

    def get_accession(self):
        """
        Gets the accession for the ResFinder hit.
        :return: The accession.
        """
        return self._accession

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            return self.__dict__ == other.__dict__
