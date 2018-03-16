import logging

import pandas

from staramr.blast.results.BlastResultsParser import BlastResultsParser
from staramr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP
from staramr.blast.results.pointfinder.nucleotide.PointfinderHitHSPRNA import PointfinderHitHSPRNA

logger = logging.getLogger('BlastResultsParserPointfinder')

"""
A Class for parsing BLAST results specific to PointFinder.
"""
logger = logging.getLogger('BlastResultsParserPointfinder')


class BlastResultsParserPointfinder(BlastResultsParser):

    def __init__(self, file_blast_map, blast_database, pid_threshold, plength_threshold):
        """
        Creates a new BlastResultsParserPointfinder.
        :param file_blast_map: A map/dictionary linking input files to BLAST results files.
        :param blast_database: The particular staramr.blast.AbstractBlastDatabase to use.
        :param pid_threshold: A percent identity threshold for BLAST results.
        :param plength_threshold: A percent length threshold for results.
        """
        super().__init__(file_blast_map, blast_database, pid_threshold, plength_threshold)

    def _create_hit(self, file, database_name, blast_record, alignment, hsp):
        logger.debug("database_name="+database_name)
        if database_name == '16S_rrsD':
            return PointfinderHitHSPRNA(file, blast_record, alignment, hsp)
        else:
            return PointfinderHitHSP(file, blast_record, alignment, hsp)

    def _create_data_frame(self, results):
        df = pandas.DataFrame(results,
                              columns=('FILE', 'GENE', 'TYPE', 'POSITION', 'MUTATION',
                                       '%IDENTITY', '%OVERLAP', 'DB_SEQ_LENGTH/QUERY_HSP'))
        return df.set_index('FILE')

    def _do_append(self, hit, db_mutation, results):
        results.append([hit.get_file(),
                        hit.get_hit_id() + " (" + db_mutation.get_mutation_string_short() + ")",
                        db_mutation.get_type(),
                        db_mutation.get_mutation_position(),
                        db_mutation.get_mutation_string(),
                        hit.get_pid(),
                        hit.get_plength(),
                        str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length())
                        ])

    def _append_results_to(self, hit, database_name, results):
        database_mutations = hit.get_mutations()

        gene = hit.get_gene()

        logger.debug("")
        logger.debug("gene=" + str(gene))
        logger.debug("sbjct_start=" + str(hit.hsp.sbjct_start))
        logger.debug("sbjct_end=" + str(hit.hsp.sbjct_end))
        logger.debug("sbjct_frame=" + str(hit.get_database_frame()))
        logger.debug("query_frame=" + str(hit.get_query_frame()))
        logger.debug("query_start=" + str(hit.hsp.query_start))
        logger.debug("query_end=" + str(hit.hsp.query_end))
        for x in database_mutations:
            logger.debug("database_mutations: position=" + str(x.get_mutation_position()) + ", mutation=" + x.get_mutation_string())

        if database_name == '16S_rrsD':
            database_resistance_mutations = self._blast_database.get_resistance_nucleotides(gene, database_mutations)
        else:
            database_resistance_mutations = self._blast_database.get_resistance_codons(gene, database_mutations)
        logger.debug("database_resistance_mutations=" + str(database_resistance_mutations))

        logger.debug("gaps=" + str(hit.hsp.gaps))

        if len(database_resistance_mutations) == 0:
            logger.debug("No mutations for [" + hit.get_hit_id() + "]")
        elif len(database_resistance_mutations) == 1:
            db_mutation = database_resistance_mutations[0]
            self._do_append(hit, db_mutation, results)
        else:
            raise Exception("Error, multiple resistance mutations for [" + hit.get_hit_id() + "], mutations " + str(
                database_resistance_mutations))
