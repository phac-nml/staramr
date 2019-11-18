import copy
import logging
import os
import pandas as pd
from os import path
import re
from collections import Counter
from typing import List, Dict, Optional

from Bio import SeqIO
from pandas import DataFrame

from staramr.blast.plasmidfinder.PlasmidfinderBlastDatabase import PlasmidfinderBlastDatabase
from staramr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from staramr.blast.results.BlastResultsParser import BlastResultsParser
from staramr.blast.results.plasmidfinder.BlastResultsParserPlasmidfinder import BlastResultsParserPlasmidfinder
from staramr.blast.results.pointfinder.BlastResultsParserPointfinder import BlastResultsParserPointfinder
from staramr.blast.results.resfinder.BlastResultsParserResfinder import BlastResultsParserResfinder
from staramr.results.AMRDetectionSummary import AMRDetectionSummary
from staramr.results.QualityModule import QualityModule

logger = logging.getLogger("AMRDetection")

"""
A Class to handle scanning files for AMR genes.
"""


class AMRDetection:

    def __init__(self, resfinder_database: ResfinderBlastDatabase, amr_detection_handler,
                 pointfinder_database: PointfinderBlastDatabase = None,
                 include_negative_results: bool = False, output_dir: str = None, genes_to_exclude: list = [],
                 plasmidfinder_database: PlasmidfinderBlastDatabase = None) -> None:
        """
        Builds a new AMRDetection object.
        :param resfinder_database: The staramr.blast.resfinder.ResfinderBlastDatabase for the particular ResFinder database.
        :param amr_detection_handler: The staramr.blast.JobHandler to use for scheduling BLAST jobs.
        :param pointfinder_database: The staramr.blast.pointfinder.PointfinderBlastDatabase to use for the particular PointFinder database.
        :param plasmidfinder_database: The staramr.blast.plasmidfinder.PlasmidfinderBlastDatabase for the particular PlasmidFinder database.
        :param include_negative_results:  If True, include files lacking AMR genes in the resulting summary table.
        :param output_dir: The directory where output fasta files are to be written into (None for no output fasta files).
        :param genes_to_exclude: A list of gene IDs to exclude from the results.
        """
        self._resfinder_database = resfinder_database
        self._amr_detection_handler = amr_detection_handler
        self._pointfinder_database = pointfinder_database
        self._plasmidfinder_database = plasmidfinder_database
        self._include_negative_results = include_negative_results

        if pointfinder_database is None:
            self._has_pointfinder = False
        else:
            self._has_pointfinder = True

        self._output_dir = output_dir

        self._genes_to_exclude = genes_to_exclude

    def _create_amr_summary(self, files: List[str], resfinder_dataframe: DataFrame,quality_module_dataframe: DataFrame,
                            pointfinder_dataframe: Optional[BlastResultsParserPointfinder],
                            plasmidfinder_dataframe: DataFrame, mlst_dataframe: DataFrame) -> DataFrame:
        amr_detection_summary = AMRDetectionSummary(files, resfinder_dataframe,quality_module_dataframe,
                                                    pointfinder_dataframe, plasmidfinder_dataframe, mlst_dataframe)
        return amr_detection_summary.create_summary(self._include_negative_results)

    def _create_detailed_amr_summary(self, files: List[str], resfinder_dataframe: DataFrame,quality_module_dataframe: DataFrame,
                                     pointfinder_dataframe: Optional[BlastResultsParserPointfinder],
                                     plasmidfinder_dataframe: DataFrame, mlst_dataframe: DataFrame) -> DataFrame:
        amr_detection_summary = AMRDetectionSummary(files, resfinder_dataframe,quality_module_dataframe,
                                                    pointfinder_dataframe, plasmidfinder_dataframe, mlst_dataframe)
        return amr_detection_summary.create_detailed_summary(self._include_negative_results)

    def _create_resfinder_dataframe(self, resfinder_blast_map: Dict, pid_threshold: float, plength_threshold: int,
                                    report_all: bool) -> DataFrame:
        resfinder_parser = BlastResultsParserResfinder(resfinder_blast_map, self._resfinder_database, pid_threshold,
                                                       plength_threshold, report_all, output_dir=self._output_dir,
                                                       genes_to_exclude=self._genes_to_exclude)
        return resfinder_parser.parse_results()

    def _create_pointfinder_dataframe(self, pointfinder_blast_map: Dict, pid_threshold: float, plength_threshold: int,
                                      report_all: bool) -> DataFrame:
        pointfinder_parser = BlastResultsParserPointfinder(pointfinder_blast_map, self._pointfinder_database,
                                                           pid_threshold, plength_threshold, report_all,
                                                           output_dir=self._output_dir,
                                                           genes_to_exclude=self._genes_to_exclude)
        return pointfinder_parser.parse_results()

    def _create_plasmidfinder_dataframe(self, plasmidfinder_blast_map: Dict[str, BlastResultsParser],
                                        pid_threshold: float, plength_threshold: int,
                                        report_all: bool) -> DataFrame:
        plasmidfinder_parser = BlastResultsParserPlasmidfinder(plasmidfinder_blast_map, self._plasmidfinder_database,
                                                               pid_threshold,
                                                               plength_threshold, report_all,
                                                               output_dir=self._output_dir,
                                                               genes_to_exclude=self._genes_to_exclude)
        return plasmidfinder_parser.parse_results()

    def create_quality_module_dataframe(self,files,genome_size_lower_bound,genome_size_upper_bound,minimum_N50_value,
                                        minimum_contig_length,unacceptable_num_contigs) ->DataFrame:
        quality_module = QualityModule(files,genome_size_lower_bound,genome_size_upper_bound,minimum_N50_value,
                                        minimum_contig_length,unacceptable_num_contigs)

        return quality_module.create_quality_module_dataframe()

    def _generate_empty_columns(self, row: list, max_cols: int, cur_cols: int) -> list:
        if(cur_cols < max_cols):
            for i in range(max_cols-cur_cols):
                row.append('-')

        return row

    def _create_mlst_dataframe(self, mlst_data: str) -> DataFrame:

        columns = ['Isolate ID', 'Scheme', 'Sequence Type']
        curr_data = []
        max_columns = 0
        extension = None

        mlst_split = mlst_data.splitlines()

        # Parse and format the current row
        for row in mlst_split:
            array_format = re.split('\t', row)
            num_columns = len(array_format)

            # We want the file name without the extension
            array_format[0] = path.basename(path.splitext(array_format[0])[0])

            if max_columns < num_columns:
                max_columns = num_columns

            curr_data.append(array_format)

        # Go through each row and append additional columns for the dataframes
        curr_data = list(map(lambda x: self._generate_empty_columns(x, max_columns, len(x)), curr_data))

        # Append Locus Column names if any
        locus_columns = max_columns - len(columns)
        if locus_columns > 0:
            for x in range(0, locus_columns):
                columns.append(("Locus {}").format(x+1))

        mlst_dataframe = pd.DataFrame(curr_data, columns=columns)
        mlst_dataframe = mlst_dataframe.set_index('Isolate ID')

        return mlst_dataframe

    def run_amr_detection(self,files, pid_threshold, plength_threshold_resfinder, plength_threshold_pointfinder,
                          plength_threshold_plasmidfinder, genome_size_lower_bound,genome_size_upper_bound,
                          minimum_N50_value,minimum_contig_length,unacceptable_num_contigs,
                          report_all=False, ignore_invalid_files=False, mlst_scheme=None) -> None:
        """
        Scans the passed files for AMR genes.
        :param files: The files to scan.
        :param pid_threshold: The percent identity threshold for BLAST results.
        :param plength_threshold_resfinder: The percent length overlap for BLAST results (resfinder).
        :param plength_threshold_pointfinder: The percent length overlap for BLAST results (pointfinder).
        :param plength_threshold_plasmidfinder: The percent length overlap for BLAST results (plasmidfinder).
        :param genome_size_lower_bound: The lower bound for the genome size as defined by the user for quality metrics
        :param genome_size_upper_bound: The upper bound for the genome size as defined by the user for quality metrics
        :param minimum_N50_value: The minimum N50 value as defined by the user for quality metrics
        :param minimum_contig_length: The minimum contig length as defined by the user for quality metrics
        :param unacceptable_num_contigs: The number of contigs in a file, equal to or above our minimum contig length, for which to raise a flag as defined by the user for quality metrics
        :param report_all: Whether or not to report all blast hits.
        :param ignore_invalid_files: Skips the invalid input files if set.
        :param mlst_scheme: Specifys scheme name MLST uses if set.
        :return: None
        """

        files_copy = copy.deepcopy(files)
        files = self._validate_files(files_copy, ignore_invalid_files)

        self._quality_module_dataframe=self.create_quality_module_dataframe(files,genome_size_lower_bound,genome_size_upper_bound,minimum_N50_value,minimum_contig_length,unacceptable_num_contigs)

        self._amr_detection_handler.run_blasts_mlst(files, mlst_scheme)

        resfinder_blast_map = self._amr_detection_handler.get_resfinder_outputs()
        self._resfinder_dataframe = self._create_resfinder_dataframe(resfinder_blast_map, pid_threshold,
                                                                     plength_threshold_resfinder, report_all)

        plasmidfinder_blast_map = self._amr_detection_handler.get_plasmidfinder_outputs()
        self._plasmidfinder_dataframe = self._create_plasmidfinder_dataframe(plasmidfinder_blast_map, pid_threshold,
                                                                             plength_threshold_plasmidfinder,
                                                                             report_all)

        mlst_data = self._amr_detection_handler.get_mlst_outputs()
        self._mlst_dataframe = self._create_mlst_dataframe(mlst_data)

        self._pointfinder_dataframe = None
        if self._has_pointfinder:
            pointfinder_blast_map = self._amr_detection_handler.get_pointfinder_outputs()
            self._pointfinder_dataframe = self._create_pointfinder_dataframe(pointfinder_blast_map, pid_threshold,
                                                                             plength_threshold_pointfinder, report_all)

        self._summary_dataframe = self._create_amr_summary(files, self._resfinder_dataframe,self._quality_module_dataframe,
                                                           self._pointfinder_dataframe, self._plasmidfinder_dataframe, self._mlst_dataframe)

        self._detailed_summary_dataframe = self._create_detailed_amr_summary(files, self._resfinder_dataframe,self._quality_module_dataframe,
                                                                             self._pointfinder_dataframe,
                                                                             self._plasmidfinder_dataframe,
                                                                             self._mlst_dataframe)

    def _validate_files(self, files: List[str], ignore_invalid_files: bool) -> List[str]:
        total_files = len(files)
        removeable_files = []

        for file in files:
            # Check if the file is not a directory
            if os.path.isdir(file):
                if ignore_invalid_files:
                    logger.warning('--ignore-invalid-files is set, skipping directory {}'.format(file))
                    removeable_files.append(file)
                else:
                    raise Exception(
                        'Directory {} is invalid, please use --ignore-invalid-files to skip over this directory'.format(
                            file))
            else:
                # Will raise an error if the input returns an empty generator on non-FASTA files, returns a boolean
                validInput = any(SeqIO.parse(file, "fasta"))

                if not validInput:
                    if ignore_invalid_files:
                        logger.warning('--ignore-invalid-files is set, skipping file {}'.format(file))
                        removeable_files.append(file)
                    else:
                        raise Exception(
                            'File {} is invalid, please use --ignore-invalid-files to skip over invalid input files'.format(
                                file))
                else:
                    # Check if there are any duplicate sequence id's in the valid files
                    record = []
                    # Store all the sequence id's in a list
                    for sequence in SeqIO.parse(file, "fasta"):
                        record.append(sequence.id)

                    duplicates = []

                    # Each sequence contains a tuple (sequence id, frequency)
                    for sequence in (Counter(record)).items():
                        if sequence[1] > 1:
                            # We want the sequence id's that are duplicates
                            duplicates.append(sequence[0])

                    # Raise an error if there's any duplicates in the file
                    if len(duplicates) > 0:
                        raise Exception(
                            'File {} contains the following duplicate sequence IDs: {}'.format(file, duplicates))

        # Check to see if the invalid file is not the only file in the directory
        if total_files == len(removeable_files):
            raise Exception('Cannot produce output due to no valid input files')

        # Remove the skipped files
        if ignore_invalid_files:
            for file in removeable_files:
                files.remove(file)

        return files

    def get_mlst_results(self):
        """
        Gets a pd.DataFrame for the MLST results.
        :return: A pd.DataFrame for the MLST results.
        """

        return self._mlst_dataframe

    def get_resfinder_results(self):
        """
        Gets a pd.DataFrame for the ResFinder results.
        :return: A pd.DataFrame for the ResFinder results.
        """
        return self._resfinder_dataframe

    def get_pointfinder_results(self):
        """
        Gets a pd.DataFrame for the PointFinder results.
        :return: A pd.DataFrame for the PointFinder results.
        """
        return self._pointfinder_dataframe

    def get_plasmidfinder_results(self):
        """
        Gets a pd.DataFrame for the PlasmidFinder results.
        :return: A pd.DataFrame for the PlasmidFinder results.
        """
        self._plasmidfinder_dataframe = self._plasmidfinder_dataframe.rename({'Gene':'Plasmid'}, axis=1)
        return self._plasmidfinder_dataframe

    def get_summary_results(self):
        """
        Gets a pd.DataFrame for a summary table of the results.
        :return: A pd.DataFrame for a summary table of the results.
        """
        return self._summary_dataframe

    def get_detailed_summary_results(self):
        """
        Gets a pd.DataFrame for a detailed summary table of the results.
        :return: A pd.DataFrame for a detailed summary table of the results.
        """

        self._detailed_summary_dataframe = self._detailed_summary_dataframe.rename({'Gene':'Data'}, axis=1)
        return self._detailed_summary_dataframe
