from staramr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP
from staramr.blast.results.pointfinder.codon.CodonMutationPosition import CodonMutationPosition
from staramr.blast.results.pointfinder.nucleotide.NucleotideMutationPosition import NucleotideMutationPosition


class PointfinderHitHSPPromoter(PointfinderHitHSP):
    """
    This class represents a Pointfinder BLAST hit and is responsible for parsing the related BLAST hit.
    """

    def __init__(self, file, blast_record, database_name):
        """
        Creates a new PointfinderHitHSPPromoter from a promoter-related BLAST hit.
        :param file: The input file.
        :param blast_record: The Bio.Blast.Record this hit came from.
        :param database_name: The name of the BLAST database. This name must be named in the same way the
            Pointfinder promoters are named (ex: "embA_promoter_size_115bp").
        """

        self._parse_database_name(database_name)
        super().__init__(file, blast_record)

    def _get_mutation_positions(self, start):
        nucleotide_mutations = []
        codon_mutations = []

        amr_seq = self.get_amr_gene_seq()
        genome_seq = self.get_genome_contig_hsp_seq()

        amr_pos = 0

        # Get all the nucleotide match positions up to the offset:
        for i in range(self.offset):

            # Insertion: "-" in the reference:
            if amr_seq[i] == "-":
                # left side
                offset = i - amr_pos + self.offset + 1  # accounting for string index and reference index possibly being different
                mutation = NucleotideMutationPosition(i, amr_seq, genome_seq, start, offset=offset + 1)
                nucleotide_mutations.append(mutation)
            # Mismatch or Deletion:
            elif (amr_seq[i] != genome_seq[i]):
                offset = i - amr_pos + self.offset
                mutation = NucleotideMutationPosition(i, amr_seq, genome_seq, start, offset=offset + 1)
                nucleotide_mutations.append(mutation)
                amr_pos += 1
            # Match:
            else:
                amr_pos += 1

        # Merge adjacent nucleotide insertions:
        nucleotide_mutations_merged = []
        while(len(nucleotide_mutations) > 0):
            current = nucleotide_mutations.pop(0)

            # If the merged list has something in it
            # and the current mutation is an insertion
            # and the last item in the merged list is an insertion
            # and the last item in the merged list has the same position:
            if(len(nucleotide_mutations_merged) > 0
            and current.get_database_amr_gene_mutation() == 'ins'
            and nucleotide_mutations_merged[-1].get_database_amr_gene_mutation() == 'ins'
            and nucleotide_mutations_merged[-1].get_mutation_position() == current.get_mutation_position()):
                # Add the nucleotide to the existing:
                nucleotide_mutations_merged[-1]._database_amr_gene_mutation += current._database_amr_gene_mutation
                nucleotide_mutations_merged[-1]._input_genome_mutation += current._input_genome_mutation

            else:
                nucleotide_mutations_merged.append(current)

        # Get all the codon match positions after the offset:
        for i in range(self.offset, len(amr_seq)):

            # Insertion: "-" in the reference:
            if amr_seq[i] == "-":
                # left side
                offset = i - amr_pos + self.offset # accounting for string index and reference index possibly being different
                mutation = CodonMutationPosition(i - 1, amr_seq, genome_seq, start, offset=offset)
                codon_mutations.append(mutation)
            # Mismatch or Deletion:
            elif (amr_seq[i] != genome_seq[i]):
                offset = i - amr_pos + self.offset
                mutation = CodonMutationPosition(i, amr_seq, genome_seq, start, offset=offset)
                codon_mutations.append(mutation)
                amr_pos += 1
            # Match:
            else:
                amr_pos += 1

        codon_mutations_filtered = []
        codon_starts = []

        # Only return codon mutation position objects with unique codon start positions
        for m in codon_mutations:
            if m._codon_start not in codon_starts:
                codon_starts.append(m._codon_start)
                codon_mutations_filtered.append(m)

        # Combine lists and return all positions:
        combined_mutations = nucleotide_mutations_merged + codon_mutations_filtered

        return combined_mutations

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

        self.offset = size
