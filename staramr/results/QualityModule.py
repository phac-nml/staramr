import logging
from os import path
from typing import Set

import pandas as pd
from pandas import DataFrame

logger = logging.getLogger("QualityModule")

"""
Summarizes the quality metric feedback into a single table.
"""


class QualityModule:

    def __init__(self, files, genome_size_lower_bound,genome_size_upper_bound,minimum_N50_value,minimum_contig_length,unacceptable_num_contigs) -> None:
        """
        Constructs an object for summarizing our quality module.
        :param files: The list of genome files we have scanned against.
        :param genome_size_lower_bound: The lower bound for the genome size as defined by the user for quality metrics
        :param genome_size_upper_bound: The upper bound for the genome size as defined by the user for quality metrics
        :param minimum_N50_value: The minimum N50 value as defined by the user for quality metrics
        :param minimum_contig_length: The minimum contig length as defined by the user for quality metrics
        :param unacceptable_num_contigs: The number of contigs under our minimum length for which to raise a flag as defined by the user for quality metrics
        """
        #self._names = [path.splitext(path.basename(x))[0] for x in files]
        self._files = files
        self._genome_size_lower_bound = genome_size_lower_bound
        self._genome_size_upper_bound = genome_size_upper_bound
        self._minimum_N50_value = minimum_N50_value
        self._minimum_contig_length = minimum_contig_length
        self._unacceptable_num_contigs = unacceptable_num_contigs
        
    
    def _create_quality_module_dataframe(self):
        """
        Goes through the files and creates a dataframe consisting of the file's genome length, N50 value and the number of contigs less than the minimum length as
        specified by the quality metrics. It also consists of the feedback for whether or not the file passed the quality metrics and if it didn't feedback on why it failed
        :return: A pd.dataframe containing the genome size, N50 value, number of contigs under our user defined minimum length
        as well as the results of our quality metrics (pass or fail) and the corresponding feedback
        """
        name_set=[]
        for myFile in self._files:
            name_set.append(path.splitext(path.basename(myFile))[0])
        files_contigs_and_genomes_lengths=self._get_files_contigs_and_genomes_lengths(self._files)
        files_genome_length_feedback = self._get_genome_length_feedback(self._files,files_contigs_and_genomes_lengths[1],self._genome_size_lower_bound,self._genome_size_upper_bound) #array where first element is the genome lengths and second is their corresponding feedback
        files_N50_value_feedback=self._get_N50_feedback(files_contigs_and_genomes_lengths[0],files_contigs_and_genomes_lengths[1],self._minimum_N50_value)
        file_num_contigs_under_minimum_bp_feedback= self._get_num_contigs_under_minimum_bp_feedback(files_contigs_and_genomes_lengths[0],self._minimum_contig_length,self._unacceptable_num_contigs)
        quality_module = self._get_quality_module(files_genome_length_feedback,files_N50_value_feedback[1],file_num_contigs_under_minimum_bp_feedback[1],self._files)
        quality_module_feedback = quality_module[0]
        quality_module_result = quality_module[1]
        quality_module_frame=pd.DataFrame([[t,u,v,w,x,y] for t,u,v,w,x,y in zip(name_set,files_contigs_and_genomes_lengths[1],files_N50_value_feedback[0],file_num_contigs_under_minimum_bp_feedback[0],quality_module_result,quality_module_feedback)],
            columns=('Isolate ID', 'Genome Length','N50 value','Number of Contigs Under '+str(self._minimum_contig_length)+' bp','Quality Module','Quality Module Feedback')).set_index('Isolate ID')
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

