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
        amr_detection_summary = AMRDetectionSummary(files, resfinder_dataframe,
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
        :param unacceptable_num_contigs: The number of contigs under our minimum length for which to raise a flag as defined by the user for quality metrics
        :param report_all: Whether or not to report all blast hits.
        :param ignore_invalid_files: Skips the invalid input files if set.
        :param mlst_scheme: Specifys scheme name MLST uses if set.
        :return: None
        """

        files_copy = copy.deepcopy(files)
        files = self._validate_files(files_copy, ignore_invalid_files)

        self._quality_module_dataframe=self._create_quality_module_dataframe(files,genome_size_lower_bound,genome_size_upper_bound,minimum_N50_value,minimum_contig_length,unacceptable_num_contigs)

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

    def _create_quality_module_dataframe(self,files,genome_size_lower_bound,genome_size_upper_bound,minimum_N50_value,minimum_contig_length,unacceptable_num_contigs):
        """
        Goes through the files and creates a dataframe consisting of the file's genome length, N50 value and the number of contigs less than the minimum length as
        specified by the quality metrics. It also consists of the feedback for whether or not the file passed the quality metrics and if it didn't feedback on why it failed
        :param genome_size_lower_bound: The lower bound for the genome size as defined by the user for quality metrics
        :param genome_size_upper_bound: The upper bound for the genome size as defined by the user for quality metrics
        :param minimum_N50_value: The minimum N50 value as defined by the user for quality metrics
        :param minimum_contig_length: The minimum contig length as defined by the user for quality metrics
        :param unacceptable_num_contigs: The number of contigs under our minimum length for which to raise a flag as defined by the user for quality metrics
        :return: A pd.dataframe containing the genome size, N50 value, number of contigs under our user defined minimum length
        as well as the results of our quality metrics (pass or fail) and the corresponding feedback
        """
        name_set=[]
        for myFile in files:
            name_set.append(path.splitext(path.basename(myFile))[0])
        files_contigs_and_genomes_lengths=self._get_files_contigs_and_genomes_lengths(files)
        files_genome_length_feedback = self._get_genome_length_feedback(files,files_contigs_and_genomes_lengths[1],genome_size_lower_bound,genome_size_upper_bound) #array where first element is the genome lengths and second is their corresponding feedback
        files_N50_value_feedback=self._get_N50_feedback(files_contigs_and_genomes_lengths[0],files_contigs_and_genomes_lengths[1],minimum_N50_value)
        file_num_contigs_under_minimum_bp_feedback= self._get_num_contigs_under_minimum_bp_feedback(files_contigs_and_genomes_lengths[0],minimum_contig_length,unacceptable_num_contigs)
        quality_module = self._get_quality_module(files_genome_length_feedback,files_N50_value_feedback[1],file_num_contigs_under_minimum_bp_feedback[1],files)
        quality_module_feedback = quality_module[0]
        quality_module_result = quality_module[1]
        quality_module_frame=pd.DataFrame([[t,u,v,w,x,y] for t,u,v,w,x,y in zip(name_set,files_contigs_and_genomes_lengths[1],files_N50_value_feedback[0],file_num_contigs_under_minimum_bp_feedback[0],quality_module_result,quality_module_feedback)],
            columns=('Isolate ID', 'Genome Length','N50 value','Number of Contigs Under '+str(minimum_contig_length)+' bp','Quality Module','Quality Module Feedback')).set_index('Isolate ID')
        return quality_module_frame

    def _get_files_contigs_and_genomes_lengths(self,files):
        #This solution was taken almost directly from https://github.com/phac-nml/sistr_cmd and was in no way specifically designed for starAMR
        """
        Goes through the files and determines their genome length as well as the length of each contig
        :param files: The files for which we wish to determine the genome length as well as the length of each contig
        :return: An array where the first element is itself an array where each element represents the corresponding 
        file and is itself an array where each element is the length for the corresponding contig inside of this file.
        The second element is itself an array where each element is the genome length for the corresponding file
        """
        files_contigs_lengths =[]
        files_genomes_lengths =[]
        feedback = []
        for filepath in files:
            genome_length = 0
            with open(filepath,'r') as g:
                contig_lengths = []
                length = 0
                for line in g:
                    line = line.strip()
                    if line == '':
                        continue
                    if line[0] == '>':
                        if length == 0:
                            continue
                        else:
                            genome_length = genome_length + length
                            contig_lengths.append(length)
                            length = 0        
                    else:
                        length = length + len(line)
            contig_lengths.append(length)
            files_contigs_lengths.append(contig_lengths)
            files_genomes_lengths.append(genome_length+length)
        feedback.append(files_contigs_lengths)
        feedback.append(files_genomes_lengths)
        return feedback

    def _get_genome_length_feedback(self,files,files_genome_lengths,lb_gsize,ub_gsize):
        """
        Goes through the files and determines whether or not they pass the quality metrics for genome length
        :param files: The files for which we wish to determine whether or not they pass the quality metrics for genome length
        :param files_genome_lengths: An array where each element is the genome length of the corresponding file
        :param lb_gsize: The lower bound for the genome size as defined by the user for quality metrics
        :param ub_gsize: The upper bound for the genome size as defined by the user for quality metrics
        :return: An array where each element corresponds to the feedback (true or false) for the corresponding file in regards to the
        genome size quality metric
        """
        files_genome_feedback=[]
        for genome_length in files_genome_lengths:
            if genome_length >= lb_gsize and genome_length <= ub_gsize:
                files_genome_feedback.append(True)
            else:
                files_genome_feedback.append(False)
        return files_genome_feedback

    def _get_N50_feedback(self,files_contigs_lengths,files_genome_lengths,minimum_N50):
        #For information on what N50 is and how it is calculated see https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics
        """
        Goes through the files and determines whether or not they pass the quality metrics for N50 value
        :param files_contigs_lengths: The lengths of the contigs for the files
        :param files_genome_lengths: An array where each element is the genome length of the corresponding file
        :param minimum_N50_value: The minimum N50 value as defined by the user for quality metrics
        :return: An array where the first element is itself an array where each element is the N50 value for 
        the corresponding file. The second element is itself an array where each element is the feedback (true or false) 
        for whether the corresponding file passes the N50 quality metrics
        """
        feedback = []
        file_index =0
        files_N50=[]
        N50_feedback = []
        for files_genome_length in files_genome_lengths:
            contig_lengths = files_contigs_lengths[file_index]
            half_length=(files_genome_lengths[file_index])/2
            contig_lengths.sort()
            contig_num=len(contig_lengths)
            contig_index = 1
            sum_of_largest_contigs=0
            while contig_index < contig_num:
                if sum_of_largest_contigs+contig_lengths[contig_num-contig_index] >=half_length:
                    break
                else: 
                    sum_of_largest_contigs=sum_of_largest_contigs+contig_lengths[contig_num-contig_index]
                    contig_index=contig_index+1
            files_N50.append(contig_lengths[contig_num-contig_index])
            file_index = file_index +1
        for file_N50_value in files_N50:
            if file_N50_value > minimum_N50:
                N50_feedback.append(True)
            else:
                N50_feedback.append(False)

        feedback.append(files_N50)
        feedback.append(N50_feedback)
        return feedback
    
    def _get_num_contigs_under_minimum_bp_feedback(self,files_contigs_lengths,minimum_contig_length,unacceptable_num_contigs_under_minimum_bp):
        """
        Goes through the files and determines whether or not they pass the quality metrics for the acceptable number of contigs under the minimum length
        :param files_contigs_lengths: The lengths of the contigs for the files
        :param minimum_contig_length: The minimum contig length as defined by the user for quality metrics
        :param unacceptable_num_contigs: The number of contigs under our minimum length for which to raise a flag as defined by the user for quality metrics
        :return: An array where the first element is itself an array where each element is the number of contigs under the minimum length for
        the corresponding file. The second element is itself an array where each element is the feedback (true or false) 
        for whether the corresponding file passes the acceptable number of contigs under the minimum length quality metric
        """
        feedback=[]
        file_num_contigs=[]
        contigs_under_minimum_bp_feedback=[]
        file_index = 0
        for files_contigs_length in files_contigs_lengths:
            num_contigs = 0
            for contig in files_contigs_lengths[file_index]:
                if contig < minimum_contig_length:
                    num_contigs = num_contigs+1
            file_num_contigs.append(num_contigs)
            file_index=file_index+1
        for file_num_contigs_under_minimum_bp in file_num_contigs:
            if file_num_contigs_under_minimum_bp >= unacceptable_num_contigs_under_minimum_bp:
                contigs_under_minimum_bp_feedback.append(False)
            else:
                contigs_under_minimum_bp_feedback.append(True)
        feedback.append(file_num_contigs)
        feedback.append(contigs_under_minimum_bp_feedback)
        return feedback

    def _get_quality_module(self,genome_length_feedback,N50_feedback,contigs_under_minimum_bp_feedback,files):
        """
        Goes through the files and for each provides detailed feedback for why they failed the quality metrics
        :param genome_length_feedback: An array where each element is the feedback (true or false) for the corresponding file in regards to the
        genome length quality metric
        :param N50_feedback: An array where each element is the feedback (true or false) for the corresponding file in regards to the
        the N50 quality metric
        :param contigs_under_minimum_bp_feedback: An array where each element is the feedback (true or false) for the corresponding file in regards to the
        the acceptable number of contigs under the minimum length quality metric
        :param files: The files for which we wish to create quality modules
        :return: An array where the first element is itself an array where each element is the detailed quality metric feedback for
        the corresponding file. The second element is itself an array where each element is the feedback (true or false) 
        for whether the corresponding file passes all of the quality metrics
        """
        file_index = 0
        feedback = []
        quality_parameter = []
        quality_parameter_feedback = []
        for file in files:
            if genome_length_feedback[file_index] == True & N50_feedback[file_index] == True & contigs_under_minimum_bp_feedback[file_index] == True:
                quality_parameter_feedback_for_file=("")
                quality_parameter.append("Passed")
            else:
                quality_parameter_feedback_for_file=""
                quality_parameter.append("Failed")
                if genome_length_feedback[file_index] == False:
                    quality_parameter_feedback_for_file = quality_parameter_feedback_for_file + "Genome length is not within the acceptable length range. "
                if N50_feedback[file_index] == False:
                    quality_parameter_feedback_for_file = quality_parameter_feedback_for_file + "N50 value is not greater than the specified minimum value. "
                if contigs_under_minimum_bp_feedback[file_index] == False:
                    quality_parameter_feedback_for_file = quality_parameter_feedback_for_file + "Number of Contigs with a length less than the minimum length exceeds the acceptable number. "

            quality_parameter_feedback.append(quality_parameter_feedback_for_file)
            file_index=file_index+1
        feedback.append(quality_parameter_feedback)
        feedback.append(quality_parameter)
        return feedback



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
