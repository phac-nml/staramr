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
        name = self._blast_record['qseqid']

        # CGE has been changing FASTA record headers to include accession
        # numbers, which need to be removed. See PointfinderHitHSP.get_amr_gene_name()
        # for more information. Naming schemes are also inconsistent:
        # pointfinder/campylobacter/23S.fsa -> 23S_1_LR134511.1
        # pointfinder/neisseria_gonorrhoeae/23S-rRNA-a1.fsa -> 23S-rRNA-a1_1_AE004969.1
        if name.startswith("16S_rrs"): name = name.split("_")[0] + "_" + name.split("_")[1]
        elif name.startswith("16S-rrs"): name = name.split("_")[0].replace("-", "_", 1) # Ex: 16S-rrsD_1_CP049983.1
        elif name.startswith("23S"): name = "23S"
        else: name = name.split("_")[0]

        return name

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
