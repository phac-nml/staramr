from staramr.blast.results.pointfinder.BlastResultsParserPointfinder import BlastResultsParserPointfinder
from staramr.blast.results.resfinder.BlastResultsParserResfinder import BlastResultsParserResfinder
from staramr.results.AMRDetectionSummary import AMRDetectionSummary


class AMRDetection:

    def __init__(self, resfinder_database, amr_detection_handler, pointfinder_database=None, include_negative_results=False):
        self._resfinder_database = resfinder_database
        self._amr_detection_handler = amr_detection_handler
        self._pointfinder_database = pointfinder_database
        self._include_negative_results = include_negative_results

        if pointfinder_database is None:
            self._has_pointfinder = False
        else:
            self._has_pointfinder = True

    def run_amr_detection(self, files, pid_threshold, plength_threshold):
        self._amr_detection_handler.run_blasts(files)

        resfinder_blast_map = self._amr_detection_handler.get_resfinder_outputs()
        resfinder_parser = BlastResultsParserResfinder(resfinder_blast_map, self._resfinder_database, pid_threshold,
                                                       plength_threshold)
        self._resfinder_dataframe = resfinder_parser.parse_results()

        if self._has_pointfinder:
            pointfinder_blast_map = self._amr_detection_handler.get_pointfinder_outputs()
            pointfinder_parser = BlastResultsParserPointfinder(pointfinder_blast_map, self._pointfinder_database,
                                                               pid_threshold, plength_threshold)
            self._pointfinder_dataframe = pointfinder_parser.parse_results()
        else:
            self._pointfinder_dataframe = None

        amr_detection_summary = AMRDetectionSummary(files, self._resfinder_dataframe, self._pointfinder_dataframe)
        self._summary_dataframe = amr_detection_summary.create_summary(self._include_negative_results)

    def get_resfinder_results(self):
        return self._resfinder_dataframe

    def get_pointfinder_results(self):
        return self._pointfinder_dataframe

    def get_summary_results(self):
        return self._summary_dataframe
