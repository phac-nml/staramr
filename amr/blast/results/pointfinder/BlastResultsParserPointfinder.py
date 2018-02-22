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
        df = pandas.DataFrame(results,
                                columns=('FILE', 'GENE', 'RESFINDER_PHENOTYPE', 'CODON_POSITION', 'NUCLEOTIDE',
                                         'AMINO_ACID', '%IDENTITY', '%OVERLAP', 'DB_SEQ_LENGTH/QUERY_HSP'))
        return df.set_index('FILE')

    def _append_results_to(self, hit, results):
        database_nucleotide_mutations = hit.get_nucleotide_mutations()
        database_codon_start_positions = [x.get_codon_start() for x in database_nucleotide_mutations]

        gene=hit.get_gene()

        logger.info("")
        logger.info("gene="+str(gene))
        logger.info("sbjct_start="+str(hit.hsp.sbjct_start))
        logger.info("sbjct_end="+str(hit.hsp.sbjct_end))
        logger.info("sbjct_frame="+str(hit.get_database_frame()))
        logger.info("query_frame=" + str(hit.get_query_frame()))
        logger.info("query_start="+str(hit.hsp.query_start))
        logger.info("query_end="+str(hit.hsp.query_end))
        logger.info("database_nucleotide_mutation_positions="+str(database_nucleotide_mutations))

        database_resistance_codons = self._blast_database.get_resistance_codons(gene, database_nucleotide_mutations)
        logger.info("database_resistance_codons="+str(database_resistance_codons))

        logger.info("gaps="+str(hit.hsp.gaps))

        if len(database_resistance_codons) == 0:
            logger.debug("No mutations for ["+hit.get_hit_id()+"]")
        elif len(database_resistance_codons) == 1:
            db_codon = database_resistance_codons[0]
            results.append([hit.get_file(),
                            hit.get_hit_id(),
                            self._blast_database.get_phenotype(hit.get_hit_id(), db_codon),
                            db_codon.get_codon_start(),
                            db_codon.get_database_codon() + ' -> ' + db_codon.get_query_codon(),
                            db_codon.get_database_amino_acid() + ' -> ' + db_codon.get_query_amino_acid(),
                            hit.get_pid(),
                            hit.get_plength(),
                            str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length())
                            ])
        else:
            raise Exception("Error, multiple resistance mutations for ["+hit.get_hit_id()+"], mutations "+str(database_resistance_codons))
