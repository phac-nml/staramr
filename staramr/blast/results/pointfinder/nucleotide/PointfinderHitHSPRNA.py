from staramr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP
from staramr.blast.results.pointfinder.nucleotide.NucleotideMutationPosition import NucleotideMutationPosition


class PointfinderHitHSPRNA(PointfinderHitHSP):

    def __init__(self, file, blast_record):
        """
        Creates a new PointfinderHitHSPRNA.
        :param file: The input file.
        :param blast_record: The Bio.Blast.Record this hit came from.
        """
        super().__init__(file, blast_record)

    def _get_mutation_positions(self, start):
        amr_seq = self.get_amr_gene_seq()
        genome_seq = self.get_genome_contig_hsp_seq()

        # @formatter:off
        return [NucleotideMutationPosition(i, amr_seq, genome_seq, start) for i in self._get_match_positions()]
        # @formatter:on
