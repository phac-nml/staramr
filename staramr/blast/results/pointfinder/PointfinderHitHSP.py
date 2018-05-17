import logging

from staramr.blast.results.AMRHitHSP import AMRHitHSP
from staramr.blast.results.pointfinder.codon.CodonMutationPosition import CodonMutationPosition

logger = logging.getLogger('PointfinderHitHSP')

"""
Class storing a PointFinder BLAST hit/hsp.
"""


class PointfinderHitHSP(AMRHitHSP):

    def __init__(self, file, blast_record):
        """
        Creates a new PointfinderHitHSP.
        :param file: The input file.
        :param blast_record: The Bio.Blast.Record this hit came from.
        """
        super().__init__(file, blast_record)

    def get_amr_gene_name(self):
        """
        Gets the particular gene name for the PointFinder hit.
        :return: The gene name.
        """
        return self._blast_record['sseqid']

    def _get_match_positions(self):
        return [i for i, (x,y) in enumerate(zip(self._blast_record['sseq'], self._blast_record['qseq'])) if x != y]

    def _get_mutation_positions(self, start, database_strand):
        return [CodonMutationPosition(i, self._blast_record['sseq'], self._blast_record['qseq'], start, database_strand) for i
                in self._get_match_positions()]

    def get_mutations(self):
        """
        Gets a list of NucleotideMutationPosition for the individual mutations.
        :return: A list of NucleotideMutationPosition.
        """
        start = self._blast_record['sstart']
        database_strand = self.get_amr_database_strand()
        return self._get_mutation_positions(start, database_strand)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        else:
            return self.__dict__ == other.__dict__
