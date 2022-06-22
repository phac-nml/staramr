from staramr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP
from staramr.blast.results.pointfinder.codon.CodonMutationPosition import CodonMutationPosition
from staramr.blast.results.pointfinder.nucleotide.NucleotideMutationPosition import NucleotideMutationPosition


class PointfinderHitHSPPromoter(PointfinderHitHSP):

    def __init__(self, file, blast_record, database_name):
        """
        Creates a new PointfinderHitHSPPromoter.
        :param file: The input file.
        :param blast_record: The Bio.Blast.Record this hit came from.
        :param database_name: The name of the BLAST database.
        """

        self._parse_database_name(database_name)
        super().__init__(file, blast_record)

    def _get_mutation_positions(self, start):
        # nucleotide space for promoter part
        # codon space for gene part

        amr_seq = self.get_amr_gene_seq()
        genome_seq = self.get_genome_contig_hsp_seq()

        match_positions = self._get_match_positions()

        # Get all the nucleotide match positions up to the offset:
        nucleotide_match_positions = list(filter(lambda x: x < self.offset, match_positions))
        print("nucleotide_match_positions = " + str(nucleotide_match_positions))

        # @formatter:off
        nucleotide_mutations = [NucleotideMutationPosition(i, amr_seq, genome_seq, start, self.offset + 1) for i in nucleotide_match_positions]
        # TODO: the +1 is weird, but seemingly needed
        # @formatter:on

        print("nucleotide_mutations = " + str(nucleotide_mutations))


        # Get all the codon match positions after the offset:
        codon_match_positions = list(filter(lambda x: x >= self.offset, match_positions))
        print("codon_match_positions = " + str(codon_match_positions))

        mutation_positions_filtered = []
        codon_starts = []

        # Only return mutation position objects with unique codon start positions
        mutation_positions = [CodonMutationPosition(i, amr_seq, genome_seq, start, self.offset) for i in codon_match_positions]

        for m in mutation_positions:
            if m._codon_start not in codon_starts:
                codon_starts.append(m._codon_start)
                mutation_positions_filtered.append(m)
        
        print("mutation_positions_filtered = " + str(mutation_positions_filtered))
        print("combined list = " + str(nucleotide_mutations + mutation_positions_filtered))

        return nucleotide_mutations + mutation_positions_filtered

    def _parse_database_name(self, database_name):
        """
        Parses the name of the database in order to obtain the promoter offset.
        The database name is expected to have the following format:

        [GENENAME]_promoter_size_[SIZE]bp

        example:

        embA_promoter_size_115bp
        """

        tokens = database_name.split("_")  # split the name into tokens
        size_string = tokens[len(tokens) - 1]  # get the last token
        size = int(size_string.replace('bp', ''))  # remove the 'bp' and convert to an int

        self.offset = size # + 1  # TODO: REALLY double-check this; I think BLAST starts on index 1, so we need to adjust for that
        print(self.offset)
