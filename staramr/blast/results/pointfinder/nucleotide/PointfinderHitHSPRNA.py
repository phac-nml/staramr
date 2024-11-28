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

    def get_amr_gene_name(self):
        """
        Gets the particular gene name for the PointfinderHitHSPRNA hit.
        :return: The gene name.
        """

        # This is explicitly overriding the parent class's function
        # because we want to return this gene name without modification.
        return self._blast_record['qseqid']

    def _get_mutation_positions(self, start):
        mutation_positions = []

        amr_seq = self.get_amr_gene_seq()
        genome_seq = self.get_genome_contig_hsp_seq()

        amr_pos = 0

        for i in range(len(amr_seq)):

            # Insertion: "-" in the reference:
            if amr_seq[i] == "-":
                # left side
                offset = i - amr_pos  # accounting for string index and reference index possibly being different
                mutation = NucleotideMutationPosition(i - 1, amr_seq, genome_seq, start, offset=offset)
                mutation_positions.append(mutation)
            # Mismatch or Deletion:
            elif (amr_seq[i] != genome_seq[i]):
                offset = i - amr_pos
                mutation = NucleotideMutationPosition(i, amr_seq, genome_seq, start, offset=offset)
                mutation_positions.append(mutation)
                amr_pos += 1
            # Match:
            else:
                amr_pos += 1

        return mutation_positions
