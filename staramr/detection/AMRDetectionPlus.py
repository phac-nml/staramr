from staramr.detection.AMRDetection import AMRDetection
from staramr.blast.results.pointfinder.BlastResultsParserPointfinderPlus import BlastResultsParserPointfinderPlus
from staramr.blast.results.resfinder.BlastResultsParserResfinderPlus import BlastResultsParserResfinderPlus
from staramr.results.AMRDetectionSummaryPlus import AMRDetectionSummaryPlus


class AMRDetectionPlus(AMRDetection):

    def __init__(self, resfinder_database, arg_drug_table_resfinder, amr_detection_handler, arg_drug_table_pointfinder,
                 pointfinder_database=None, include_negative_results=False, output_dir=None):
        super().__init__(resfinder_database, amr_detection_handler, pointfinder_database, include_negative_results, output_dir=output_dir)
        self._arg_drug_table_resfinder = arg_drug_table_resfinder
        self._arg_drug_table_pointfinder = arg_drug_table_pointfinder

    def _create_resfinder_dataframe(self, resfinder_blast_map, pid_threshold, plength_threshold, report_all):
        resfinder_parser = BlastResultsParserResfinderPlus(resfinder_blast_map, self._arg_drug_table_resfinder,
                                                           self._resfinder_database, pid_threshold,
                                                           plength_threshold, report_all, output_dir=self._output_dir)
        return resfinder_parser.parse_results()

    def _create_pointfinder_dataframe(self, pointfinder_blast_map, pid_threshold, plength_threshold, report_all):
        pointfinder_parser = BlastResultsParserPointfinderPlus(pointfinder_blast_map, self._arg_drug_table_pointfinder,
                                                               self._pointfinder_database,
                                                               pid_threshold, plength_threshold, report_all, output_dir=self._output_dir)
        return pointfinder_parser.parse_results()

    def _create_amr_summary(self, files, resfinder_dataframe, pointfinder_dataframe):
        amr_detection_summary = AMRDetectionSummaryPlus(files, resfinder_dataframe, pointfinder_dataframe)
        return amr_detection_summary.create_summary(self._include_negative_results)
