import logging

import pandas

from staramr.blast.results.BlastResultsParser import BlastResultsParser
from staramr.blast.results.pointfinder.PointfinderHitHSP import PointfinderHitHSP

logger = logging.getLogger('BlastResultsParserPointfinder')

"""
A Class for parsing BLAST results specific to PointFinder.
"""


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

    def _create_hit(self, file, blast_record, alignment, hsp):
        return PointfinderHitHSP(file, blast_record, alignment, hsp)

    def _create_data_frame(self, results):
        df = pandas.DataFrame(results,
                              columns=('FILE', 'GENE', 'RESFINDER_PHENOTYPE', 'CODON_POSITION', 'NUCLEOTIDE',
                                       'AMINO_ACID', '%IDENTITY', '%OVERLAP', 'DB_SEQ_LENGTH/QUERY_HSP'))
        return df.set_index('FILE')

    def _do_append(self, hit, db_codon, results):
        results.append([hit.get_file(),
                        hit.get_hit_id() + " (" + db_codon.get_database_amino_acid() + str(
                            db_codon.get_codon_start()) + db_codon.get_query_amino_acid() + ")",
                        self._blast_database.get_phenotype(hit.get_hit_id(), db_codon),
                        db_codon.get_codon_start(),
                        db_codon.get_database_codon() + ' -> ' + db_codon.get_query_codon(),
                        db_codon.get_database_amino_acid() + ' -> ' + db_codon.get_query_amino_acid(),
                        hit.get_pid(),
                        hit.get_plength(),
                        str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length())
                        ])

    def _append_results_to(self, hit, database_name, results):
        database_nucleotide_mutations = hit.get_nucleotide_mutations()

        gene = hit.get_gene()

        logger.debug("")
        logger.debug("gene=" + str(gene))
        logger.debug("sbjct_start=" + str(hit.hsp.sbjct_start))
        logger.debug("sbjct_end=" + str(hit.hsp.sbjct_end))
        logger.debug("sbjct_frame=" + str(hit.get_database_frame()))
        logger.debug("query_frame=" + str(hit.get_query_frame()))
        logger.debug("query_start=" + str(hit.hsp.query_start))
        logger.debug("query_end=" + str(hit.hsp.query_end))
        for x in database_nucleotide_mutations:
            logger.debug("database_nucleotide_mutation_positions codon=" + str(
                x.get_database_codon()) + ", aa=" + x.get_query_amino_acid())

        database_resistance_codons = self._blast_database.get_resistance_codons(gene, database_nucleotide_mutations)
        logger.debug("database_resistance_codons=" + str(database_resistance_codons))

        logger.debug("gaps=" + str(hit.hsp.gaps))

        if len(database_resistance_codons) == 0:
            logger.debug("No mutations for [" + hit.get_hit_id() + "]")
        elif len(database_resistance_codons) == 1:
            db_codon = database_resistance_codons[0]
            self._do_append(hit, db_codon, results)
        else:
            raise Exception("Error, multiple resistance mutations for [" + hit.get_hit_id() + "], mutations " + str(
                database_resistance_codons))
