from staramr.blast.results.pointfinder.BlastResultsParserPointfinder import BlastResultsParserPointfinder
from staramr.blast.results.resfinder.BlastResultsParserResfinder import BlastResultsParserResfinder
from staramr.results.AMRDetectionSummary import AMRDetectionSummary

"""
A Class to handle scanning files for AMR genes.
"""


class AMRDetection:

    def __init__(self, resfinder_database, amr_detection_handler, pointfinder_database=None,
                 include_negative_results=False, output_dir=None):
        """
        Builds a new AMRDetection object.
        :param resfinder_database: The staramr.blast.resfinder.ResfinderBlastDatabase for the particular ResFinder database.
        :param amr_detection_handler: The staramr.blast.BlastHandler to use for scheduling BLAST jobs.
        :param pointfinder_database: The staramr.blast.pointfinder.PointfinderBlastDatabase to use for the particular PointFinder database.
        :param include_negative_results:  If True, include files lacking AMR genes in the resulting summary table.
        :param output_dir: The directory where output fasta files are to be written into (None for no output fasta files).
        """
        self._resfinder_database = resfinder_database
        self._amr_detection_handler = amr_detection_handler
        self._pointfinder_database = pointfinder_database
        self._include_negative_results = include_negative_results

        if pointfinder_database is None:
            self._has_pointfinder = False
        else:
            self._has_pointfinder = True

        self._output_dir = output_dir

    def _create_amr_summary(self, files, resfinder_dataframe, pointfinder_dataframe):
        amr_detection_summary = AMRDetectionSummary(files, resfinder_dataframe,
                                                    pointfinder_dataframe)
        return amr_detection_summary.create_summary(self._include_negative_results)

    def _create_resfinder_dataframe(self, resfinder_blast_map, pid_threshold, plength_threshold, report_all):
        resfinder_parser = BlastResultsParserResfinder(resfinder_blast_map, self._resfinder_database, pid_threshold,
                                                       plength_threshold, report_all, output_dir=self._output_dir)
        return resfinder_parser.parse_results()

    def _create_pointfinder_dataframe(self, pointfinder_blast_map, pid_threshold, plength_threshold, report_all):
        pointfinder_parser = BlastResultsParserPointfinder(pointfinder_blast_map, self._pointfinder_database,
                                                           pid_threshold, plength_threshold, report_all,
                                                           output_dir=self._output_dir)
        return pointfinder_parser.parse_results()

    def run_amr_detection(self, files, pid_threshold, plength_threshold_resfinder, plength_threshold_pointfinder,
                          report_all=False):
        """
        Scans the passed files for AMR genes.
        :param files: The files to scan.
        :param pid_threshold: The percent identity threshold for BLAST results.
        :param plength_threshold_resfinder: The percent length overlap for BLAST results (resfinder).
        :param plength_threshold_pointfinder: The percent length overlap for BLAST results (pointfinder).
        :param report_all: Whether or not to report all blast hits.
        :return: None
        """
        self._amr_detection_handler.run_blasts(files)

        resfinder_blast_map = self._amr_detection_handler.get_resfinder_outputs()
        self._resfinder_dataframe = self._create_resfinder_dataframe(resfinder_blast_map, pid_threshold,
                                                                     plength_threshold_resfinder, report_all)

        if self._has_pointfinder:
            pointfinder_blast_map = self._amr_detection_handler.get_pointfinder_outputs()
            self._pointfinder_dataframe = self._create_pointfinder_dataframe(pointfinder_blast_map, pid_threshold,
                                                                             plength_threshold_pointfinder, report_all)
        else:
            self._pointfinder_dataframe = None

        self._summary_dataframe = self._create_amr_summary(files, self._resfinder_dataframe,
                                                           self._pointfinder_dataframe)

    def get_resfinder_results(self):
        """
        Gets a pandas.DataFrame for the ResFinder results.
        :return: A pandas.DataFrame for the ResFinder results.
        """
        return self._resfinder_dataframe

    def get_pointfinder_results(self):
        """
        Gets a pandas.DataFrame for the PointFinder results.
        :return: A pandas.DataFrame for the PointFinder results.
        """
        return self._pointfinder_dataframe

    def get_summary_results(self):
        """
        Gets a pandas.DataFrame for a summary table of the results.
        :return: A pandas.DataFrame for a summary table of the results.
        """
        return self._summary_dataframe
