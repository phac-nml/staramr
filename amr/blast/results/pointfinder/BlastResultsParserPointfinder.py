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
        database_nucleotide_mutation_positions = hit.get_database_nucleotide_mutation_positions()
        database_codon_mutation_positions = hit.get_codon_mutation_positions_at(database_nucleotide_mutation_positions)
        gene=hit.get_gene()

        database_resistance_codon_mutation_positions = self._blast_database.get_resistance_codon_mutation_positions(gene, database_codon_mutation_positions)
        database_resistance_codons_reference = hit.get_database_nucleotide_codons_at(database_nucleotide_mutation_positions)
        database_resistance_codons = self._blast_database.get_resistance_codons(gene, database_codon_mutation_positions)

        logger.info("database_resistance_codons_reference="+str(database_resistance_codons_reference))
        logger.info("database_resistance_codons="+str(database_resistance_codons))

        if (database_resistance_codons_reference != database_resistance_codons):
            raise Exception("Error, did not extract identical codons from both blast database file and resfinder info file. blast="+
                            str(database_resistance_codons_reference) + ", resfinder info="+str(database_resistance_codons))
        elif len(database_resistance_codon_mutation_positions) == 0:
            logger.debug("No mutations for ["+hit.get_hit_id()+"]")
        elif len(database_resistance_codon_mutation_positions) == 1:
            results.append([hit.get_file(),
                            hit.get_hit_id(),
                            database_resistance_codon_mutation_positions[0],
                            database_resistance_codons[0],
                            hit.get_pid(),
                            hit.get_plength(),
                            str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length())
                            ])
        else:
            raise Exception("Error, multiple resistance mutations for ["+hit.get_hit_id()+"], mutations "+str(database_resistance_codon_mutation_positions))
