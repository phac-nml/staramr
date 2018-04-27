import logging

import pandas

from staramr.blast.results.pointfinder.BlastResultsParserPointfinder import BlastResultsParserPointfinder

logger = logging.getLogger('BlastResultsParserPointfinderPlus')


class BlastResultsParserPointfinderPlus(BlastResultsParserPointfinder):

    def __init__(self, file_blast_map, arg_drug_table, blast_database, pid_threshold, plength_threshold, report_all=False, output_dir=None):
        super().__init__(file_blast_map, blast_database, pid_threshold, plength_threshold, report_all, output_dir=output_dir)
        self._arg_drug_table = arg_drug_table

    def _create_data_frame(self, results):
        df = pandas.DataFrame(results,
                              columns=('Isolate ID', 'Gene', 'Predicted Phenotype', 'Type', 'Position', 'Mutation',
                                       '%Identity', '%Overlap', 'HSP Length/Total Length', 'Contig', 'Start', 'End'))
        return df.set_index('Isolate ID')

    def _do_append(self, hit, db_mutation, results):
        drug = self._arg_drug_table.get_drug(self._blast_database.get_organism(), hit.get_hit_id(),
                                             db_mutation.get_mutation_position())
        results.append([hit.get_isolate_id(),
                        hit.get_hit_id() + " (" + db_mutation.get_mutation_string_short() + ")",
                        drug,
                        db_mutation.get_type(),
                        db_mutation.get_mutation_position(),
                        db_mutation.get_mutation_string(),
                        hit.get_pid(),
                        hit.get_plength(),
                        str(hit.get_hsp_alignment_length()) + "/" + str(hit.get_alignment_length()),
                        hit.get_contig(),
                        hit.get_contig_start(),
                        hit.get_contig_end()
                        ])
