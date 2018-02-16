import pandas
import logging

from amr.blast.results.BlastResultsParser import BlastResultsParser
from amr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP

logger = logging.getLogger('BlastResultsParserPointfinder')

class BlastResultsParserPointfinder(BlastResultsParser):

    def __init__(self, file_blast_map, blast_database, pid_threshold, plength_threshold):
        super().__init__(file_blast_map, blast_database, pid_threshold, plength_threshold)

    def _create_hit(self, file, blast_record, alignment, hsp):
        return PointfinderHitHSP(file, blast_record, alignment, hsp)

    def _create_data_frame(self, results):
        return pandas.DataFrame(results,
                                columns=('FILE', 'GENE', 'CODON_POSITION', 'NUCLEOTIDE',
                                         '%IDENTITY', '%OVERLAP', 'DB_SEQ_LENGTH/QUERY_HSP'))

    def _append_results_to(self, hit, results):
        nucleotide_mutation_positions = hit.get_nucleotide_mutation_positions()
        codon_mutation_positions = hit.get_codon_mutation_positions_at(nucleotide_mutation_positions)
        codons = hit.get_nucleotide_codons_at(nucleotide_mutation_positions)
        gene=hit.get_gene()

        resistance_codon_mutation_positions = self._blast_database.get_resistance_codon_mutation_positions(gene, codon_mutation_positions)
        resistance_codons = self._blast_database.get_resistance_codons(gene, codon_mutation_positions)

        if len(resistance_codon_mutation_positions) == 0:
            logger.debug("No mutations for ["+hit.get_hit_id()+"]")
        elif len(resistance_codon_mutation_positions) == 1:
            results.append([hit.get_file(),
                            hit.get_hit_id(),
                            resistance_codon_mutation_positions[0],
                            resistance_codons[0],
                            hit.get_pid(),
                            hit.get_plength(),
                            str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length())
                            ])
        else:
            raise Exception("Error, multiple resistance mutations for ["+hit.get_hit_id()+"], mutations "+str(resistance_codon_mutation_positions))
