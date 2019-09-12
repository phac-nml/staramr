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
        return self._blast_record['qseqid']

    def _get_match_positions(self):
        amr_seq = self.get_amr_gene_seq()
        genome_seq = self.get_genome_contig_hsp_seq()

        return [i for i, (x, y) in enumerate(zip(amr_seq, genome_seq)) if x != y]

    def _get_mutation_positions(self, start):
        mutation_positions_filtered = []
        codon_starts = []

        amr_seq = self.get_amr_gene_seq()
        genome_seq = self.get_genome_contig_hsp_seq()

        # Only return mutation position objects with unique codon start positions
        mutation_positions = [CodonMutationPosition(i, amr_seq, genome_seq, start) for i in self._get_match_positions()]

        for m in mutation_positions:
            if m._codon_start not in codon_starts:
                codon_starts.append(m._codon_start)
                mutation_positions_filtered.append(m)

        # @formatter:off
        return mutation_positions_filtered
        # @formatter:on

    def get_mutations(self):
        """
        Gets a list of NucleotideMutationPosition for the individual mutations.
        :return: A list of NucleotideMutationPosition.
        """
        return self._get_mutation_positions(self.get_amr_gene_start())
